# LIFE750 Advanced Genomic Analysis
## Assignment Cycle 2 Scientific Report
__Assignment Overview__
Work individually to produce a scientific report that presents a __bioinformatic analysis and interpretation__ of the supplied data set. The focus of the report is to demonstrate how bioinformatic data can be analysed to generate biologically meaningful insights.
The goal is to characterise a mutated epigenomic factor (Gene X) and investigate its effect on genomic-wide gene expression using multiplegenomic data types.


## Project Description
We suspect an epigenomic factor (Gene X) is important for gene expression.
To test this, we intergrated four data types:
- A genome assembly (FASTA) and annotation (GFF)
- RUN-seq count data from 5 normal and 5 mutant individuals
- CUT&RUN binding data (BED file) showing Gene X genomic binding sites
- DNA sequencing data from a cloned Gene X plasmid (normal and mutant)


## Objectives
1. Characterise the mutations present in the mutant Gene X sequence
2. Identify genomic regions bound by Gene X using CUT&RUN data
3. Identify differentially expressed genes between normal and mutant individuals
4. Interate binding and expression data to identify direct Gene X regulatory tagets
   
## Backgroud
Epigenomic factors regulate gene expression by controlling chromatin accessibility and transcription factor binding. CUT&RUN (Cleavage Under Tagets and Release Using Nuclease) is a high-resolution method for mapping protein-DNA interactions genomic-wide (Skene, P.J. and Henikoff, S., 2017). RNA-seq combined with DESeq2 enables quantification of transcriptional differences between conditions (Love et al., 2014).

## Dataset
Download dataset from 
```bash
wget -R "index.html*" -r -nd -nH -np https://cgr.liv.ac.uk/454/acdarby/LIFE750/Life750_datasets/dataset_tbpsrisa/
```

## Dataset Contents 
- `genome.fa`: Genome assembly in FASTA format
- `gene.gff3`: Genome annotation in GFF3 format
- `Gene_reference.fa`: Gene X reference sequence
- `normal_GeneX_R1.fasta/_R2.fasta`: Paired-end DNA reads - normal Gene X
- `mutant_GeneX_R1.fasta/_R2.fasta`: Paired-end DNA reads - mutant Gene X
- `gene_counts`: RNA-seq counts
- `gene_x_binding_sites.bad`: CUT&RUN peaks for Gene X binding

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

## Repository Contents
- `README.md`
- `Dataset` folder contains: Gene X generated files
- `Results` folder contains:
    - `raw_qc_results` folder contains: mutant_report and normal_report
    - `alignment` folder contains: mutant_GeneX_sorted.bam and normal_GeneX_sorted.bam
    - `variants` folder contains: gene_x_variants.vcf
    - `binding` folder contains: binding_sites_nearest_genes.txt, gene_with_proximal_binding.txt, and bound_gene_ids.txt
    - `differential_expression` folder contains: DESeq2_results.csv, significant_DEGs.csv, bound_and_DE_genes.csv, PCA_plot.png, and MA_plot.png

## Reference
Baptista, T., Grünberg, S., Minoungou, N., Koster, M.J.E., Timmers, H.T.M., Hahn, S., Devys, D. and Tora, L. (2017) ‘Saga is a general cofactor for rna polymerase ii transcription’, Molecular Cell, 68(1), pp. 130-143.e5. Available at: https://doi.org/10.1016/j.molcel.2017.08.016.

Churchman, L.S. and Weissman, J.S. (2011) ‘Nascent transcript sequencing visualizes transcription at nucleotide resolution’, Nature, 469(7330), pp. 368–373. Available at: https://doi.org/10.1038/nature09652.

Gkikopoulos, T., Schofield, P., Singh, V., Pinskaya, M., Mellor, J., Smolle, M., Workman, J.L., Barton, G.J. and Owen-Hughes, T. (2011) ‘A role for snf2-related nucleosome-spacing enzymes in genome-wide nucleosome organization’, Science, 333(6050), pp. 1758–1760. Available at: https://doi.org/10.1126/science.1206097.

Kacsoh, B.Z., Lynch, Z.R., Mortimer, N.T. and Schlenke, T.A. (2013) ‘Fruit flies medicate offspring after seeing parasites’, Science, 339(6122), pp. 947–950. Available at: https://doi.org/10.1126/science.1229625.

Nojima, T., Gomes, T., Grosso, A.R.F., Kimura, H., Dye, M.J., Dhir, S., Carmo-Fonseca, M. and Proudfoot, N.J. (2015) ‘Mammalian net-seq reveals genome-wide nascent transcription coupled to rna processing’, Cell, 161(3), pp. 526–540. Available at: https://doi.org/10.1016/j.cell.2015.03.027.

Rhee, H.S. and Pugh, B.F. (2011) ‘Comprehensive genome-wide protein-dna interactions detected at single-nucleotide resolution’, Cell, 147(6), pp. 1408–1419. Available at: https://doi.org/10.1016/j.cell.2011.11.013.

Sanz, L.A., Hartono, S.R., Lim, Y.W., Steyaert, S., Rajpurkar, A., Ginno, P.A., Xu, X. and Chédin, F. (2016) ‘Prevalent, dynamic, and conserved r-loop structures associate with specific epigenomic signatures in mammals’, Molecular Cell, 63(1), pp. 167–178. Available at: https://doi.org/10.1016/j.molcel.2016.05.032.

Skene, P.J. and Henikoff, S. (2017) An efficient targeted nuclease strategy for high-resolution mapping of DNA binding sites, eLife. Available at: https://doi.org/10.7554/eLife.21856.

Sun, M., Schwalb, B., Pirkl, N., Maier, K.C., Schenk, A., Failmezger, H., Tresch, A. and Cramer, P. (2013) ‘Global analysis of eukaryotic mrna degradation reveals xrn1-dependent buffering of transcript levels’, Molecular Cell, 52(1), pp. 52–62. Available at: https://doi.org/10.1016/j.molcel.2013.09.010.

Wissink, E.M., Vihervaara, A., Tippens, N.D. and Lis, J.T. (2019) ‘Nascent RNA analyses: tracking transcription and its regulation’, Nature Reviews Genetics, 20(12), pp. 705–723. Available at: https://doi.org/10.1038/s41576-019-0159-6.

## Software reference
Bioinformatics tools and packages for coding

### R-studio packages
Love, M.I., Huber, W. and Anders, S. (2014) ‘Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2’, Genome Biology, 15(12), p. 550. Available at: https://doi.org/10.1186/s13059-014-0550-8.

Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T., Miller, E., Bache, S., Müller, K., Ooms, J., Robinson, D., Seidel, D., Spinu, V., Takahashi, K., Vaughan, D., Wilke, C., Woo, K. and Yutani, H. (2019) ‘Welcome to the tidyverse’, Journal of Open Source Software, 4(43), p. 1686. Available at: https://doi.org/10.21105/joss.01686.

### Bioinformatics tools
Andrews, S. (2019). Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data. [online] Babraham.ac.uk. Available at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

Chen, S., Zhou, Y., Chen, Y. and Gu, J. (2018) ‘Fastp: an ultra-fast all-in-one fastq preprocessor’, Bioinformatics, 34(17), pp. i884–i890. Available at: https://doi.org/10.1093/bioinformatics/bty560.

Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M. and Li, H. (2021) ‘Twelve years of samtools and bcftools’, GigaScience, 10(2), p. giab008. Available at: https://doi.org/10.1093/gigascience/giab008.

Ewels, P., Magnusson, M., Lundin, S. and Käller, M. (2016) ‘MultiQC: summarize analysis results for multiple tools and samples in a single report’, Bioinformatics, 32(19), pp. 3047–3048. Available at: https://doi.org/10.1093/bioinformatics/btw354.

Li, H. and Durbin, R. (2009) ‘Fast and accurate short read alignment with Burrows–Wheeler transform’, Bioinformatics, 25(14), pp. 1754–1760. Available at: https://doi.org/10.1093/bioinformatics/btp324.

Quinlan, A.R. and Hall, I.M. (2010) ‘BEDTools: a flexible suite of utilities for comparing genomic features’, Bioinformatics, 26(6), pp. 841–842. Available at: https://doi.org/10.1093/bioinformatics/btq033.
