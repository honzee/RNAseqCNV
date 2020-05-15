#' Model data preparation.
#'
#' This function prepares the expression and SNV data for modelling the CNA
#'
#' @param config R script assigning paths to needed directories into variables: 1. count_dir - path to a directory with count files, 2. snv_dir - path to a directory with files with snv information
#' (either vcf or custom tabular data), out_dir - path to an output directory. More detailed description can be found in the package README file.
#' @param metadata path to a metadata table with three columns. First colum: sample names, second column: file names of count files, third column: file names of snv files. There should be no header.
#'  More information is included in the package README file.
#' @param snv_format character string, either "vcf" or "custom". "vcf" argument should be used when vcf files with snv information are generated with GATK. Otherwise "custom" arguments can be used when input
#' with snv iformation has 4 required columns: chromosome, locus of the snv, overall sequencing depth of the locus and MAF. MAF is the proportion of alternative allele sequencing depth to overall sequencing depth of the locus.
#' @param CNV_data path to a table with known CNVs on arm-level. The information is used to generate CNV-gene expression level correlation coeffiecient, that are
#' subsequently used for gene weight calculation. The more samples are provided, the more robust the calculated correlation coeffiecient will be.
#' There are two optional input formats. Either in wide format, where each column represents one chromosome arm e.g. columns (1p,1q,2p,2q,3p,3q..) and
#' one additional column with sample name. The values in chromosome arm column are the CNV states on that chromosome arm, if the CNV state is unknown the value should be set to NA. Also the CNV_data file
#' can be in "long" format. Here the table is required to have three columns: sample: sample name, arm: chromosome arm (format: "1p", "19q", "Xp") and CNV: (CNV state, e.g. +1, -1, +2, 0)
#' @param CNV_data_format string, specifies in which format is the CNV_data provided, can be either "wide" or "long"
#' @param output_file path to a file, where the training data will be saved.
#' @param referData table, reference data for gene annotation with ensamble ids
#' @param keptSNP vector of realiable SNPs to keep for the MAF graphs
#' @param par_reg table with pseudoautosomal regions. These regions will be filtered out.
#' @param centr_refer table with chromosomal centromeric locations.
#' @param weight_tab table with per-gene weight for calculating weighted quantiles for the boxplots in the main figure.
#' @param generate_weights logical value, if TRUE, weights for calculating weighted quantiles will be contructed from variance and depth of the analyzed cohort of samples. If batch is TRUE, the weights will be analyzed
#' For more information, please see the README file.
#' @param chroms vector of chromosomes to be analyzed.
#' @param batch logical value, if TRUE, the samples will be normalized together as a batch, also gene expression median will be calculated from these samples
#' @param standard_samples character vector with sample names of samples which should be used as a standard for vst and log2 fold centering. The samples names must be included in the metadata table and batch analysis cannot be TRUE. If NULL (default), in-build standard samples will be used.
#' @param minDepth minimal depth of of SNV to be kept (default 20).
#' @param minReadCnt numeric value value used for filtering genes with low expression according to to formula: at least samp_prop*100 percent of samples have more reads than minReadCnt. (default 3)
#' @param samp_prop sample proportion which is required to have at least minReadCnt reads for a gene. The samples inlcude the diploid reference (from standard_samples parameter) and analyzed sample. (default 0.8)
#' @param weight_samp_prop proportion of samples with highest weight to be kept. default (1)

get_train_data <- function(config, metadata, snv_format, CNV_data, CNV_data_format, output_file, referData = RNAseqCNV:::refDataExp, keptSNP = RNAseqCNV:::keepSNP, par_reg = RNAseqCNV:::par_reg, centr_refer = centr_ref, weight_tab = RNAseqCNV:::weight_table, generate_weights = FALSE, chroms = RNAseqCNV:::chrs, batch = FALSE, standard_samples = NULL,
                          minDepth = 3, minReadCnt = 3, samp_prop = 0.8, weight_samp_prop = 1) {

  #check whether a snv_format was selected
  if (is.null(snv_format) | !snv_format %in% c("vcf", "custom")) {
    stop("snv_format parameter has to be either 'vcf' or 'custom'")
  }

  source(config, local = TRUE)
  if (is.null(out_dir) | is.null(count_dir) | is.null(snv_dir)) {
    stop("Incorrect config file format")
  } else if (!dir.exists(out_dir) | !dir.exists(snv_dir) | !dir.exists(count_dir)) {
    stop("Directory from config file does not exist")
  }

  #check metadata file
  metadata_tab = fread(metadata, header = FALSE)
  if (ncol(metadata_tab) != 3) {
    stop("The number of columns in metadata table should be 3")
  }

  if(!CNV_data_format %in% c("long", "wide")) {
    stop("CNV_data_format has to be set to either long or wide format")
  }


  CNV_table <- fread(CNV_data)

  if(CNV_data_format == "wide") {
    CNV_table <- pivot_longer(data = CNV_table, names_to = "arm", values_to = "CNV", cols = c(paste0(c(1:22, "X"), "q"), paste0(c(1:12, 16:20, "X"), "p")))
  }

  #Create sample table
  sample_table = metadata_tab %>% mutate(count_path = file.path(count_dir, pull(metadata_tab, 2)), snv_path = file.path(snv_dir, pull(metadata_tab, 3)))

  #check whether any of the files is missing
  count_check <- file.exists(sample_table$count_path)
  snv_check <- file.exists(sample_table$snv_path)

  if (any(!c(count_check, snv_check))) {
    count_miss <- sample_table$count_path[!count_check]
    snv_miss <- sample_table$snv_path[!snv_check]

    files_miss <- paste0(c(count_miss, snv_miss), collapse = ", ")

    stop(paste0("File/s: ", files_miss, " not found"))
  }

  if (!is.null(standard_samples)) {

    #Check whether both standard samples and batch analysis were not selected together
    if (batch == TRUE & !is.null(standard_samples)) {
      stop("Both batch analysis and single sample analysis with selected standard samples cannot be performed together. Either select batch as FALSE or do not input standard samples")
    }
    #Check whether the samples are present in the sample table
    if (all(standard_samples %in% sample_table[, 1]) == FALSE) {
      stop("The input standard samples are not in metadata table.")
    }
    #Create standard sample table
    standard_samples <- create_standard(standard_samples = standard_samples, sample_table = sample_table)

  } else {
    standard_samples <- RNAseqCNV:::diploid_standard[, c(1:20, 41)]
  }

  CNV_table <- fread(CNV_data)

  if(CNV_data_format == "wide") {
    CNV_table <- pivot_longer(data = CNV_table, names_to = "arm", values_to = "CNV", cols = c(paste0(c(1:22, "X"), "q"), paste0(c(1:12, 16:20, "X"), "p")))
  }

  #If batch analysis was selected normalize the input samples together
  if(batch == TRUE) {

    print("Normalizing gene expression and applying variance stabilizing transformation...")

    #calculate normalized count values with DESeq2 normalization method for batch of samples from the input
    count_norm <- get_norm_exp(sample_table = sample_table, sample_num = 1, standard_samples = standard_samples, minReadCnt = minReadCnt, samp_prop = samp_prop, weight_table = weight_tab, weight_samp_prop = weight_samp_prop, batch = TRUE, generate_weights)

    #calculate median gene expression across diploid reference and analyzed sample for batch of samples from the input
    pickGeneDFall <- get_med(count_norm = count_norm, refDataExp = referData, generate_weights = generate_weights, weight_table = weight_tab)

    if (generate_weights == TRUE) {
      #create weights based on variance of gene expression and expression depth of the batch of samples
      weight_tab <- create_weights(pickGeneDFall)
    }
  }



  #Run the analysis for every sample in the table
  for(i in 1:nrow(sample_table)) {

    sample_name <- as.character(sample_table[i, 1])

    #normalize the samples with in-build standard or standard from the input and calculate gene medians
    if(batch == FALSE) {

      #calculate normalized count values with DESeq2 normalization method
      count_norm <- get_norm_exp(sample_table = sample_table, sample_num = i, standard_samples = standard_samples, minReadCnt = minReadCnt, samp_prop = samp_prop, weight_table = weight_tab, weight_samp_prop = weight_samp_prop, batch, generate_weights)

      #calculate median gene expression across diploid reference and analyzed sample
      pickGeneDFall <- get_med(count_norm = count_norm, refDataExp = referData, generate_weights = generate_weights)

      if (generate_weights == TRUE) {
        #create weights based on variance of gene expression and expression depth of the batch of samples
        weight_tab <- create_weights(pickGeneDFall)
      }
    }

    # if count file format was incorrect print out a message and skip this sample
    if (is.character(count_norm)) {
      message(count_norm)
      next()
    }

    #load SNP data
    smpSNP <- prepare_snv(sample_table = sample_table, sample_num = i, centr_ref = centr_ref, chrs = chroms, snv_format = snv_format, minDepth = mindepth)

    # if SNV data format was incorrect print out a message and skip this sample
    if (is.character(smpSNP[[1]])) {
      message(smpSNP[[1]])
      next()
    }

    #filter SNP data base on dpSNP database
    smpSNPdata.tmp <- filter_snv(smpSNP[[1]], keepSNP = keptSNP, minDepth = minDepth)

    #analyze chromosome-level metrics (out-dated)
    smpSNPdata <- calc_chrom_lvl(smpSNPdata.tmp)

    #arm-level metrics
    smpSNPdata_a_2 <- calc_arm(smpSNPdata.tmp)

    #select sample
    count_norm_samp <- count_norm %>% select(!!quo(sample_name)) %>% mutate(ENSG = rownames(.))

    #join reference data and weight data
    count_ns <- count_transform(count_ns = count_norm_samp, pickGeneDFall, refDataExp = referData, weight_table = weight_tab)

    #remove PAR regions
    count_ns <- remove_par(count_ns = count_ns, par_reg = par_reg)

    #Calculate metrics for chromosome arms
    feat_tab <- get_arm_metr(count_ns = count_ns, smpSNPdata = smpSNPdata_a_2, sample_name = sample_names, centr_ref = centr_ref, chrs = chrs)

    # merge with standard
    train_s = cbind(sample = sample_name, feat_tab) %>% mutate(chr = as.character(chr), arm = paste0(chr, arm))

    CNV_table_s = CNV_table %>% filter(sample == sample_name)
    train_s = train_s %>% right_join(select(CNV_table_s, sample, arm, CNV), by = c("sample", "arm")) %>% filter(chr != "Y") %>% mutate(chr_status = ifelse(CNV == 0, "dipl", "no_dipl")) %>% metr_dipl()

    if (i == 1) {
      train_data <- train_s
    } else {
      train_data = rbind(train_data, train_s)
    }
  }

  write.table(train_data, file =  output_file, sep = "\t", row.names = FALSE)
}
