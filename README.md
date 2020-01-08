---
title: "README"
author: "Honzik"
date: "January 7, 2020"
output: html_document
---

# RNAseqCNA(V)app tutorial for Zhaohui

To download the package from gitub, use devtools functions.
```
library(devtools)
install_github(repo = "honzee/RNAseqCNAapp")
```

## Input
To run the app or wrapper, config file and metadata file is needed.

Metadata parameter is a path to a file containing a table with three columns. The first column with sample names, the second with count file names and the third snv file names.

```
# Can be saved and used for testing purpouses
structure(list(sample = c("SJALL014946_D1", "SJALL014949_D1", 
"SJALL014950_D1", "SJALL014951_D1", "SJALL014954_D1", "SJALL015619_D1", 
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

Config file parameter is a path to an R script defining paths needed for the analysis. An example of a script below:
```
out_dir = "C:/Users/honza/Dropbox/St.Jude files/shiny/output_test"
count_dir = "C:/Users/honza/Dropbox/St.Jude files/HTSeq"
snv_dir = "C:/Users/honza/Dropbox/St.Jude files/snv"

```


## Running the shiny app

To run the shiny app:

```
library(RNAseqCNAapp)

launchApp()
```

## Running the wrapper

The analysis has also a wrapper

```
RNAseqCNA_wrapper(config = "/some/config", metadata = "some/metadata")
```

