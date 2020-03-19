# RNAseqCNV

### A tool in R for analysis of copy number variations from RNA-seq data
The aim of this R package is to enable analysis, visualization and automatic estimation of large-scale (arm-level) copy number variations (CNVs) from RNA-seq data. Users can use either a wrapper function or a shiny app to create clear figures and automatically estimate CNVs on chromosome arm level. The app serves also as an interface to view and check the reliability of the results.

### Installation
Users must have [R](https://www.r-project.org/) installed. [Rstudio](https://rstudio.com/products/rstudio/) is optional but recommended IDE.

To download RNAseqCNV package from GitHub, it is convenient to use the package devtools, namely function: install_github. The package installation can take a few minutes.
```
# install devtools
install.packages("devtools")

# install RNAseqCNV package
devtools::install_github(repo = "honzee/RNAseqCNV")
```

### Functionality
The results are generated either by a wrapper function: RNAseqCNV_wrapper() or through a shiny app which is deployed by launchApp() function. The RNAseqCNV_wrapper() provides more flexibility in terms of function arguments. The app on the other hand enables easier browsing and checking of the results. 

```
# Examples of basic function calls:
library(RNAseqCNV)

# example of wrapper function (DO NOT RUN) 
RNAseqCNV_wrapper(config = "path/to/config", metadata = "path/to/metadata", snv_format = "vcf")

# launch the shiny app with:
launchApp()
```

#### Input
Both the wrapper and the shiny app receive the same required input. Per-gene read count and SNV minor allele frequency (MAF) and depth are used to produce the results. Therefore, two types of files are needed for each sample.

To use the wrapper or analyze samples inside the shiny application, path to **config file** and **metadata file** in correct formats need to be provided.

##### Config
Config parameter expects a path to an R script as its argument. This R script defines paths to the input directories for both count files and files with snv information and also to the output directory. The names of variables inside the script must be identical as those in the example below:
```
out_dir = "/Path/to/output_dir"
count_dir = "/Path/to/dir/with/count_files"
snv_dir = "/Path/to/dir/with/vcf_files"

```

##### Metadata
Metadata parameter expects a path to a file with comma/tab/space separated table with three columns. The first column contains sample names, the second countains count file names and the third contains vcf/custom table file names. **The table cannot have a header and the order of these columns must be kept as in the example below:**

[]()|  | 
--- | --- | ---
SJALL014946_D1 | SJALL014946_D1.HTSeq | SJALL014946_D1.vcf
SJALL014949_D1 | SJALL014949_D1.HTSeq | SJALL014949_D1.vcf
SJALL014950_D1 | SJALL014950_D1.HTSeq | SJALL014950_D1.vcf

##### Input file formats:

##### 1. Per-gene read count

An output of HTSeq or similar read counting software. The table must have two columns: first column with ensembl gene ids and second column with read count. The table should not have a header.

| []() |    |
|-----------------|-----|
| ENSG00000000005 | 0   |
| ENSG00000000419 | 94  |
| ENSG00000000457 | 128 |
| ENSG00000000460 | 133 |
| ENSG00000000938 | 171 |
| ENSG00000000971 | 5   |

##### 2. SNV MAF and depth information

Either vcf file or custom tabular data

- vcf file with correct format can be acquired by running [GATK pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-).

| #CHROM | POS   | ID | REF | ALT | QUAL   | FILTER | INFO                                                                                                                                                                              | FORMAT         | sample_name           |
|-------|-------|----|-----|-----|--------|--------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------|--------------------------|
| 1     | 14792 | .  | G   | A   | 156.77 | .      | AC=1;AF=0.500;AN=2;BaseQRankSum=-1.442;ClippingRankSum=0.000;DP=27;ExcessHet=3.0103;FS=1.690;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=5.81;ReadPosRankSum=-0.169;SOR=1.124 | GT:AD:DP:GQ:PL | 0/1:19,8:27:99:185,0,644 |
| 1     | 14907 | .  | A   | G   | 126.77 | .      | AC=1;AF=0.500;AN=2;BaseQRankSum=-0.015;ClippingRankSum=0.000;DP=10;ExcessHet=3.0103;FS=2.808;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=14.09;ReadPosRankSum=0.751;SOR=1.170 | GT:AD:DP:GQ:PL | 0/1:3,6:9:64:155,0,64    |
| 1     | 14930 | .  | A   | G   | 161.77 | .      | AC=1;AF=0.500;AN=2;BaseQRankSum=0.751;ClippingRankSum=0.000;DP=10;ExcessHet=3.0103;FS=2.808;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=17.97;ReadPosRankSum=2.515;SOR=1.170  | GT:AD:DP:GQ:PL | 0/1:3,6:9:64:190,0,64    |

- custom tabular data must have four required columns: chr - chromosome the SNV is located on, start - locus of the SNV, depth - read depth for this locus, maf - minor allele frequency, can be calculated as depth for alternative (minor) allele divided by overall read depth of the locus. The header names are required to respect the format below:

| chr | start | depth | maf    |
|-----|-------|-------|--------|
| 1   | 14599 | 40    | 0.5   |
| 1   | 14604 | 9    | 0.3333   |
| 1   | 14610 | 10   | 0.25 |

#### Output
The output (figures and tables) of both wrapper and the shiny app will be saved in the output directory as specified in the config file.

##### Main figure

![main figure 1](./README/main_fig_1.png)

The main figure consists of two panels.

The upper panel shows the visualization of per-chromosom expression level. Y axis shows log2 fold change of expression against reference samples, x axis is divided into 23 separete facets, on each for chromosomes from 1-X (chromosome Y purpousely excluded) and the position along the x axis represents the position of genes on a chromosome. For each chromosome a weighted boxplot is drawn based upon the distribution of normalized expression of genes on that chromosome. The median expression value of a chromosome is also represented by the color of a boxplot on scale from blue (low median of expression), white (median expression around 0) to red(high median of expression). For each chromosome a pair of random forrest models estimates the copy number. This estimation can be seen in the upper part of each facet. The red colour of this estimation label signifies lower confidence (quality) of CNV call. 

The bottom panel shows the density graphs of MAF for each chromosome. It is important to note that only MAF values in the interval from 0.05 to 0.9 were kept, since the SNVs with values out of this range are not helpful in determining CNVs. In the upper part of each density graph there is also a peak distance number. It is the distance on x axis between two highest peaks in the density graph. This can help in distinguishing between CNVs and also copy neutral loss of heterozygozity (LOH).

##### Arm-level figures

![arm level figure](./README/arm_level.png)

Users have the option to generate close up figures of each chromosome with either

```
RNAseqCNV_wrapper(config = "path/to/config", metadata = "path/to/metadata", snv_format = "vcf", arm_lvl = TRUE)
```
or ticking the appropriate box in the shiny app interface.

The large panel in the middle is a close up of the main figure, specific for one chromosome. In the upper part in addition to the random forest estimated alteration, there is also the percentage of trees in the model that agreed upon this alteration. On both sides of this panel, there are two MAF density graphs, one for p arm and one for q arm. For chromosomes without p arm there is only one side panel on the left.

##### Estimation table
The estimated gender, alterations and chromosome number are saved in a table in the output directory.

| sample       | gender | chrom_n | alterations                                                                                 |
|--------------|--------|---------|---------------------------------------------------------------------------------------------|
| SJHYPER141_D | male   | 59      | 1q+, 4+, 5+, 6+, 7+, ?8+, ?8+, 10+, 12+, 14+, 14+, ?16q, 17+, ?18+, ?18+, 21+, 21+, 22+, X+ |
 

