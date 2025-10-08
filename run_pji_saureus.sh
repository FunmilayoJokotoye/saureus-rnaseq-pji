#!/usr/bin/env bash
# --------------------------
# 0) Basic project structure
# --------------------------
PROJECT="PJI_Saureus_rnaseq"
mkdir -p "${PROJECT}"/{fastq,trimmed,fastqc,align,counts,ref,meta,plots,scripts,logs}
cd "${PROJECT}"
# ---------------------------------
# 1) (Optional) create conda env
# ---------------------------------
# Uncomment if you want this script to build the environment for you.
# mamba create -y -n pji-saureus \
#   sra-tools fastp fastqc multiqc hisat2 samtools subread \
#   r-base r-deseq2 r-ggplot2 r-pheatmap r-data.table r-tidyverse \
#   r-biocmanager r-EnhancedVolcano r-ComplexHeatmap r-scales r-viridis
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate pji-saureus

# ---------------------------------
# 2) Study metadata & sample table
# ---------------------------------
# Acute:  SRR20959680, SRR20959681, SRR20959682
# Chronic:SRR20959676, SRR20959677, SRR20959678, SRR20959679
cat > meta/samples.tsv <<'TSV'
sample_id	condition	srr
A1_acute	acute	SRR20959680
A2_acute	acute	SRR20959681
A3_acute	acute	SRR20959682
C1_chronic	chronic	SRR20959676
C2_chronic	chronic	SRR20959677
C3_chronic	chronic	SRR20959678
C4_chronic	chronic	SRR20959679
TSV

# ---------------------------------
# 3) Reference genome & annotation
# ---------------------------------
# Provide your reference FASTA + GFF/GTF here (download once and reuse).
# Examples (replace with the exact files you want to use):
#   wget -O ref/genome.fa "https://.../Saureus_ref.fa"
#   wget -O ref/annotation.gff "https://.../Saureus_ref.gff"
#
# REQUIRED paths:
REF_FASTA="ref/genome.fa"
REF_GFF="ref/annotation.gff"

if [[ ! -s "${REF_FASTA}" || ! -s "${REF_GFF}" ]]; then
  echo ">> ERROR: Missing reference files. Place your S. aureus reference at:"
  echo "          ${REF_FASTA}"
  echo "          ${REF_GFF}"
  echo "   Then re-run this script."
  exit 1
fi

# Build HISAT2 index if not present
if [[ ! -e ref/sa_index.1.ht2 ]]; then
  echo ">> Building HISAT2 index..."
  hisat2-build -p 8 "${REF_FASTA}" ref/sa_index > logs/hisat2_build.log 2>&1
fi

# ---------------------------------
# 4) Download FASTQ from SRA
# ---------------------------------
echo ">> Downloading FASTQ from SRA..."
cut -f3 meta/samples.tsv | tail -n +2 | while read -r SRR; do
  if compgen -G "fastq/${SRR}*.fastq.gz" > /dev/null; then
    echo "   - Found existing fastq for ${SRR}, skipping."
    continue
  fi
  fasterq-dump "${SRR}" -O fastq --threads 8 --progress > "logs/${SRR}_fqdump.log" 2>&1
  # gzip all produced fastqs
  find fastq -maxdepth 1 -type f -name "${SRR}*.fastq" -print0 | xargs -0 -I{} pigz -p 8 {}
done

# Determine if runs are paired-end or single-end by checking files
IS_PAIRED=0
if ls fastq/*_1.fastq.gz >/dev/null 2>&1; then
  IS_PAIRED=1
fi
echo ">> Library type auto-detect: $([[ ${IS_PAIRED} -eq 1 ]] && echo 'PAIRED-END' || echo 'SINGLE-END')"

# ---------------------------------
# 5) QC + TRIMMING
# ---------------------------------
echo ">> Running FastQC (pre-trim)..."
fastqc -t 8 -o fastqc fastq/*.fastq.gz > logs/fastqc_pre.log 2>&1

echo ">> Trimming with fastp..."
if [[ ${IS_PAIRED} -eq 1 ]]; then
  # Paired-end
  for r1 in fastq/*_1.fastq.gz; do
    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    base=$(basename "${r1/_1.fastq.gz/}")
    fastp \
      -i "${r1}" -I "${r2}" \
      -o "trimmed/${base}_1.trimmed.fastq.gz" \
      -O "trimmed/${base}_2.trimmed.fastq.gz" \
      --length_required 30 --thread 8 \
      --html "fastqc/${base}.fastp.html" \
      --json "fastqc/${base}.fastp.json" > "logs/${base}_fastp.log" 2>&1
  done
else
  # Single-end
  for fq in fastq/*.fastq.gz; do
    base=$(basename "${fq/.fastq.gz/}")
    fastp \
      -i "${fq}" -o "trimmed/${base}.trimmed.fastq.gz" \
      --length_required 30 --thread 8 \
      --html "fastqc/${base}.fastp.html" \
      --json "fastqc/${base}.fastp.json" > "logs/${base}_fastp.log" 2>&1
  done
fi

echo ">> MultiQC summary..."
multiqc -o fastqc fastqc > logs/multiqc.log 2>&1

# ---------------------------------
# 6) Alignment (HISAT2) + BAM sort
# ---------------------------------
echo ">> Aligning reads with HISAT2..."
while read -r sid cond srr; do
  [[ "${sid}" == "sample_id" ]] && continue
  if [[ ${IS_PAIRED} -eq 1 ]]; then
    r1="trimmed/${srr}_1.trimmed.fastq.gz"
    r2="trimmed/${srr}_2.trimmed.fastq.gz"
    [[ ! -s "${r1}" || ! -s "${r2}" ]] && { echo "Missing ${r1} or ${r2}"; exit 1; }
    hisat2 -p 8 -x ref/sa_index -1 "${r1}" -2 "${r2}" \
      2> "logs/${sid}_hisat2.log" \
      | samtools view -bS - \
      | samtools sort -@ 8 -o "align/${sid}.sorted.bam"
  else
    r="trimmed/${srr}.trimmed.fastq.gz"
    [[ ! -s "${r}" ]] && { echo "Missing ${r}"; exit 1; }
    hisat2 -p 8 -x ref/sa_index -U "${r}" \
      2> "logs/${sid}_hisat2.log" \
      | samtools view -bS - \
      | samtools sort -@ 8 -o "align/${sid}.sorted.bam"
  fi
  samtools index "align/${sid}.sorted.bam"
done < meta/samples.tsv

# ---------------------------------
# 7) featureCounts (gene-level)
# ---------------------------------
# Adjust -t/-g to match your GFF. Common for bacteria:
#   -t gene  -g ID        OR
#   -t gene  -g locus_tag OR
#   -t CDS   -g Parent    (less common for whole genes)
#
# Start with -t gene -g ID; change if needed after inspecting your GFF attributes.
echo ">> Running featureCounts..."
featureCounts -T 8 -s 0 \
  -a "${REF_GFF}" -g ID -t gene \
  -o counts/featureCounts_raw.txt align/*.sorted.bam > logs/featureCounts.log 2>&1

# Clean header: strip 'align/' from column names for easier matching
awk 'BEGIN{FS=OFS="\t"} NR==2{for(i=7;i<=NF;i++) gsub(/align\//,"",$i)} {print}' \
  counts/featureCounts_raw.txt > counts/featureCounts_clean.txt

# ---------------------------------
# 8) DE analysis (DESeq2) via heredoc R
# ---------------------------------
echo ">> Running DESeq2 (R)..."
cat > scripts/deseq2_saureus.R <<'RSCRIPT'
suppressPackageStartupMessages({
  library(data.table); library(DESeq2); library(ggplot2)
  library(pheatmap); library(EnhancedVolcano); library(viridis)
})

count_file <- "counts/featureCounts_clean.txt"
coldata_file <- "meta/coldata.tsv"
outdir <- "plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# load counts
dt <- fread(count_file)
anno <- dt[, 1:6]
cts  <- as.matrix(dt[, 7:ncol(dt)])
rownames(cts) <- dt$Geneid

# clean colnames
colnames(cts) <- sub("\\.sorted\\.bam$", "", colnames(cts))
colnames(cts) <- sub("\\.bam$", "", colnames(cts))

# load metadata
coldata <- fread(coldata_file)
stopifnot(all(coldata$sample_id %in% colnames(cts)))
cts <- cts[, coldata$sample_id, drop=FALSE]
coldata$condition <- factor(coldata$condition, levels=c("chronic","acute"))

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData   = as.data.frame(coldata),
                              design    = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

# VST/PCA
vsd <- vst(dds, blind=FALSE)
pca <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=3) + geom_text(nudge_y = 1, size=3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  scale_color_viridis(discrete=TRUE, end=0.8) +
  theme_bw(base_size=12)
ggsave(file.path(outdir, "PCA_condition.png"), p, width=6.5, height=5.2, dpi=300)

# Results: acute vs chronic (reference = chronic)
res <- results(dds, contrast=c("condition","acute","chronic"))
# lfcShrink for stability
if ("apeglm" %in% rownames(installed.packages())) {
  res <- lfcShrink(dds, coef="condition_acute_vs_chronic", res=res, type="apeglm")
}
res_dt <- as.data.table(res, keep.rownames="gene")
fwrite(res_dt[order(padj)], "counts/DESeq2_results_acute_vs_chronic.tsv", sep="\t")

# Volcano
keyp <- res_dt[!is.na(padj)]
vp <- EnhancedVolcano(keyp,
  lab = keyp$gene,
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,
  FCcutoff = 1.0,
  labSize = 3.0,
  pointSize = 2.0,
  drawConnectors = TRUE,
  legendPosition = 'right')
ggsave(file.path(outdir, "Volcano_acute_vs_chronic.png"), vp, width=7.5, height=6, dpi=300)

# Top-30 heatmap
top <- keyp[order(padj)][1:min(30, nrow(keyp))]$gene
mat <- assay(vsd)[top, , drop=FALSE]
mat <- t(scale(t(mat)))
ann <- data.frame(condition = coldata$condition)
rownames(ann) <- coldata$sample_id
pheatmap(mat, annotation_col = ann, show_rownames=TRUE, cluster_cols=TRUE,
         filename = file.path(outdir, "Heatmap_top30.png"), width=7, height=9)
RSCRIPT

Rscript scripts/deseq2_saureus.R > logs/deseq2.log 2>&1

# ---------------------------------
# 9) Done: print key outputs
# ---------------------------------
echo "============================================================"
echo " DONE. Key outputs:"
echo "  - Counts matrix:        counts/featureCounts_clean.txt"
echo "  - DE results (TSV):     counts/DESeq2_results_acute_vs_chronic.tsv"
echo "  - PCA plot:             plots/PCA_condition.png"
echo "  - Volcano plot:         plots/Volcano_acute_vs_chronic.png"
echo "  - Heatmap (top 30):     plots/Heatmap_top30.png"
echo
echo " If featureCounts yielded too few counts:"
echo "   * Inspect ${REF_GFF} attributes; adjust -t/-g in featureCounts (e.g., -t gene -g locus_tag)."
echo "   * Verify strandedness (-s 0/1/2)."
echo
echo " (Optional) Enrichment:"
echo "   * Build a mapping (gene -> GO/KEGG) via eggNOG/UniProt/KEGGREST, then use clusterProfiler::enricher()."
echo "============================================================"
```
