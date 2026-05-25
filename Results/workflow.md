## Bioinformatics Workflow
### Step 1 - Environment Setup
```bash
conda create -n 750
conda activate 750
conda install -c bioconda -c conda-forge \
  bwa samtools bcftools snippy bedtools subread r-base fastqc multiqc fastp -y
```

### Step 2 - Quality Control
```bash
# FastQC on all raw reads
mkdir -p ~/750/results/raw_qc_results
fastqc *.fastq -o ~/750/results/raw_qc_results/

# MultiQC summary
cd ~/750/results/raw_qc_results
multiqc .

# fastp trimming — mutant
fastp -i mutant_GeneX_R1.fastq -I mutant_GeneX_R2.fastq \
  -o ~/750/results/raw_qc_results/mutant_GeneX_R1.fastq \
  -O ~/750/results/raw_qc_results/mutant_GeneX_R2.fastq \
  --html ~/750/results/raw_qc_results/mutant_report.html

# fastp trimming — normal
fastp -i normal_GeneX_R1.fastq -I normal_GeneX_R2.fastq \
  -o ~/750/results/raw_qc_results/normal_GeneX_R1.fastq \
  -O ~/750/results/raw_qc_results/normal_GeneX_R2.fastq \
  --html ~/750/results/raw_qc_results/normal_report.html
```

### Step 3 - Alignment and Variant Calling
```bash
# Index reference
bwa index GeneX_reference.fa

# Align reads and sort
mkdir -p ~/750/results/alignment

bwa mem -t 2 GeneX_reference.fa normal_GeneX_R1.fastq normal_GeneX_R2.fastq \
  | samtools sort -@ 2 -o ~/750/results/alignment/normal_GeneX_sorted.bam

bwa mem -t 2 GeneX_reference.fa mutant_GeneX_R1.fastq mutant_GeneX_R2.fastq \
  | samtools sort -@ 2 -o ~/750/results/alignment/mutant_GeneX_sorted.bam

# Index BAM files
samtools index ~/750/results/alignment/normal_GeneX_sorted.bam
samtools index ~/750/results/alignment/mutant_GeneX_sorted.bam
```

```bash
# Variant calling
mkdir -p ~/750/results/variants

bcftools mpileup -f GeneX_reference.fa \
  ~/750/results/alignment/normal_GeneX_sorted.bam \
  ~/750/results/alignment/mutant_GeneX_sorted.bam \
  | bcftools call -mv -Ob -o ~/750/results/variants/gene_x_variants.bcf

# Convert to VCF for viewing
bcftools view ~/750/results/variants/gene_x_variants.bcf \
  > ~/750/results/variants/gene_x_variants.vcf
```

### Step 4 - CUT&RUN Binding Analysis
```bash
mkdir -p ~/750/results/binding

# Sort files
bedtools sort -i gene_x_binding_sites.bed > gene_x_binding_sites_sorted.bed
bedtools sort -i genes.gff3 > genes_sorted.gff3

# Find nearest genes to binding peaks
bedtools closest -a gene_x_binding_sites_sorted.bed -b genes_sorted.gff3 \
  > ~/750/results/binding/binding_sites_nearest_genes.txt

# Find genes within 2000 bp of binding peaks
bedtools window -a genes_sorted.gff3 -b gene_x_binding_sites_sorted.bed -w 2000 \
  > ~/750/results/binding/genes_with_proximal_binding.txt

# Extract unique bound gene IDs
grep -o "ID=GENE_[0-9]*" ~/750/results/binding/genes_with_proximal_binding.txt \
  | sort | uniq > ~/750/results/binding/bound_gene_ids.txt

wc -l ~/750/results/binding/bound_gene_ids.txt
```

### Step 5 - Differential Expression Analysis (R-studio/DESeq2)
```r
library(DESeq2)
library(tidyverse)

# Load counts data
counts_data <- read.csv("~/750/dataset/gene_counts.txt",
                        sep="\t", row.names=1)

# Create sample metadata
colData <- data.frame(
  condition = factor(c(rep("normal", 5), rep("mutant", 5)))
)
rownames(colData) <- colnames(counts_data)

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ condition)

# Filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Set reference level
dds$condition <- relevel(dds$condition, ref = "normal")

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]

# Save results
res_df <- as.data.frame(res)
write.csv(res_df,
  "~/750/results/differential_expression/DESeq2_results.csv")

# Significant DEGs
sig <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig,
  "~/750/results/differential_expression/significant_DEGs.csv")

# Plots
vsd <- vst(dds, blind=FALSE)
png("~/750/results/differential_expression/PCA_plot.png",
    width=800, height=600)
plotPCA(vsd, intgroup="condition")
dev.off()

png("~/750/results/differential_expression/MA_plot.png",
    width=800, height=600)
plotMA(res, ylim=c(-5,5))
dev.off()
```

### Step 6 - Integration: Binding vs Differential Expression
```r
# Load bound gene IDs
bound_genes <- read.table(
  "~/750/results/binding/bound_gene_ids.txt",
  header=FALSE, stringsAsFactors=FALSE)
bound_genes$V1 <- gsub("ID=", "", bound_genes$V1)

# Overlap with significant DEGs
overlap <- intersect(rownames(sig), bound_genes$V1)
cat("Genes bound AND differentially expressed:", length(overlap), "\n")

write.csv(data.frame(gene=overlap),
  "~/750/results/differential_expression/bound_and_DE_genes.csv")
```

