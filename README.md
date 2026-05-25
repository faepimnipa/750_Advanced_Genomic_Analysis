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
Epigenomic factors regulate gene expression by controlling chromatin accessibility and transcription factor binding. CUT&RUN (Cleavage Under Tagets and Release Using Nuclease) is a high-resolution method for mapping protein-DNA interactions genomic-wide (
## Dataset
Download dataset from 
```bash
wget -R "index.html*" -r -nd -nH -np https://cgr.liv.ac.uk/454/acdarby/LIFE750/Life750_datasets/dataset_tbpsrisa/
```
## Bioinformatics Workflow


## Repository Contents
- `README.md`
- `Dataset` folder contains: Gene X generated files
- `Results` folder contains:
    - `raw_qc_results

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

# R-studio packages
Love, M.I., Huber, W. and Anders, S. (2014) ‘Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2’, Genome Biology, 15(12), p. 550. Available at: https://doi.org/10.1186/s13059-014-0550-8.

Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T., Miller, E., Bache, S., Müller, K., Ooms, J., Robinson, D., Seidel, D., Spinu, V., Takahashi, K., Vaughan, D., Wilke, C., Woo, K. and Yutani, H. (2019) ‘Welcome to the tidyverse’, Journal of Open Source Software, 4(43), p. 1686. Available at: https://doi.org/10.21105/joss.01686.

# Bioinformatics tools
Andrews, S. (2019). Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data. [online] Babraham.ac.uk. Available at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

Chen, S., Zhou, Y., Chen, Y. and Gu, J. (2018) ‘Fastp: an ultra-fast all-in-one fastq preprocessor’, Bioinformatics, 34(17), pp. i884–i890. Available at: https://doi.org/10.1093/bioinformatics/bty560.

Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M. and Li, H. (2021) ‘Twelve years of samtools and bcftools’, GigaScience, 10(2), p. giab008. Available at: https://doi.org/10.1093/gigascience/giab008.

Ewels, P., Magnusson, M., Lundin, S. and Käller, M. (2016) ‘MultiQC: summarize analysis results for multiple tools and samples in a single report’, Bioinformatics, 32(19), pp. 3047–3048. Available at: https://doi.org/10.1093/bioinformatics/btw354.

Li, H. and Durbin, R. (2009) ‘Fast and accurate short read alignment with Burrows–Wheeler transform’, Bioinformatics, 25(14), pp. 1754–1760. Available at: https://doi.org/10.1093/bioinformatics/btp324.

Quinlan, A.R. and Hall, I.M. (2010) ‘BEDTools: a flexible suite of utilities for comparing genomic features’, Bioinformatics, 26(6), pp. 841–842. Available at: https://doi.org/10.1093/bioinformatics/btq033.
