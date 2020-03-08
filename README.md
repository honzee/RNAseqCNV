# RNAseqCNV tutorial

To download the package from gitub, it is convenient to use the package devtools, namely function: install_github. The package installation can take a few minutes.
```
library(devtools)
install_github(repo = "honzee/RNAseqCNV")
```

## Input
The package uses expression level and SNV information to estimate CNV of chromosomes and chromosomal arms. Each samples must have two input files. 
1. An output of HTSeq or similar read counting software, with the output having two columns: one with ensembl gene names and the second column with count number.
2. .snv file. An SNV file can be acquired by running GATK pipeline to get .vcf file. .vcf file can be subsequntly converted into .snv with perl script included in the package: "/path_to/library/RNAseqCNVapp/inst/vcf_to_snv.pl.txt"

To run the shiny application or the wrapper for generating figures, **config file** and **metadata file** is needed.

### Metadata
Metadata parameter is a path to a file containing a table with three columns. The first column contains sample names, the second countains count file names and the third contains snv file names.


**Example**

sample | count | snv
--- | --- | ---
SJALL014946_D1 | SJALL014946_D1.HTSeq | SJALL014946_D1.snv
SJALL014949_D1 | SJALL014949_D1.HTSeq | SJALL014949_D1.snv
SJALL014950_D1 | SJALL014950_D1.HTSeq | SJALL014950_D1.snv
SJALL014951_D1 | SJALL014951_D1.HTSeq | SJALL014951_D1.snv
SJALL014946_D1 | SJALL014954_D1.HTSeq | SJALL014954_D1.snv

```
# Can be saved and used for testing purpouses
structure(list(sample = c("SJALL014946_D1", "SJALL014949_D1", 
"SJALL014950_D1", "SJALL014951_D1", "SJALL014946_D1", "SJALL015619_D1", 
"SJALL015927_D1", "SJALL015940_D1", "SJALL015942_D1", "SJALL015943_D1"
), count = c("SJALL014946_D1.HTSeq", "SJALL014949_D1.HTSeq", 
"SJALL014950_D1.HTSeq", "SJALL014951_D1.HTSeq", "SJALL014954_D1.HTSeq", 
"SJALL015619_D1.HTSeq", "SJALL015927_D1.HTSeq", "SJALL015940_D1.HTSeq", 
"SJALL015942_D1.HTSeq", "SJALL015943_D1.HTSeq"), snv = c("SJALL014946_D1.snv", 
"SJALL014949_D1.snv", "SJALL014950_D1.snv", "SJALL014951_D1.snv", 
"SJALL014954_D1.snv", "SJALL015619_D1.snv", "SJALL015927_D1.snv", 
"SJALL015940_D1.snv", "SJALL015942_D1.snv", "SJALL015943_D1.snv"
)), row.names = c(NA, 10L), class = "data.frame")

```

### Config

Config file parameter is a path to an R script, which is defining paths needed for the analysis. An example of a script below:
```
out_dir = "C:/Users/honza/Dropbox/St.Jude files/shiny/output_test"
count_dir = "C:/Users/honza/Dropbox/St.Jude files/HTSeq"
snv_dir = "C:/Users/honza/Dropbox/St.Jude files/snv"

```

## Running the shiny app

To run the shiny app, first load the RNAseqCNVapp package and then call function: launchApp()

```
library(RNAseqCNVapp)

launchApp()
```
In the app itself, load your metadata file and config file and then you should be able to run the analysis.
 
## Running the wrapper

The analysis has a wrapper, that allows to analyse the samples without running the shiny app. The two compuylsory parameters are config and metadata.

```
RNAseqCNV_wrapper(config = "/some/config", metadata = "some/metadata")
```

