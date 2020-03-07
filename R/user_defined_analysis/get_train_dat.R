#' Model data preparation.
#'
#' This function prepares the expression and SNV data for modelling the CNA
#'
#' @param config
#' @param metadata
#' @param weight_table
#' @param standard
#' @param out_dir

get_train_dat <- function(config, metadata, weight_table = weight_table, standard) {

  #read input
  metadata <- fread(metadata, header = TRUE)
  source(config)

  #create sample_table
  sample_table <- metadata

  count_f = pull(sample_table, 2)
  snv_f = pull(sample_table, 3)
  #create paths to the files
  sample_table$count_path <- file.path(count_dir, count_f)
  sample_table$snv_path <-  file.path(snv_dir, snv_f)

  #training data frame
  train_data <- data.frame(sample = c(), chr = c(), arm = c(), medianvst_nor_med = c(), up_quart = c(), low_quart = c(), peak_max = c(), peak_m_dist = c(), peakdist = c(), sd = c(), sds_median = c(), sds_025 = c(), sds_075 = c(), n_02_04 = c(), n_04 = c(), alteration = c(), sds_median_dipl = c(), sds_025_dipl = c(), sds_075_dipl = c())

  for (i in 1:nrow(sample_table)) {

    # normalize gene expression
    count_norm <- get_norm_exp(sample_table = sample_table, sample_num = i, diploid_standard = dipl_standard, minReadCnt = minReadCnt, samp_prop = samp_prop, weight_table = weight_tab, weight_samp_prop = weight_samp_prop)

    pickGeneDFall <- get_med(count_norm = count_norm, refDataExp = refDataExp)

    smpSNP <- prepare_snv(sample_table = sample_table, sample_num = i, centr_ref = centr_ref, minDepth = 20, chrs = chrs)

    #filter SNP data base on dpSNP database
    smpSNPdata.tmp <- filter_snv(smpSNP[[1]], keepSNP = keepSNP)

    #analyze chromosome-level metrics
    smpSNPdata <- calc_arm(smpSNPdata.tmp)

    sample_name <- as.character(sample_table[i, 1])

    count_ns <- select(vst, !!quo(sample_name)) %>% mutate(ENSG = rownames(vst))

    #join reference data and weight data
    count_ns <- count_transform(count_ns = count_ns, pickGeneDFall, refData, weight_table)

    #remove PAR regions
    count_ns <- remove_par(count_ns = count_ns, par_reg = par_reg)

    summ_arm_train <- get_arm_metr(count_ns = count_ns, smpSNPdata = smpSNPdata, sample_name = sample_name, centr_ref = centr_ref, chrs = chrs)

    # merge with standard
    train_s = cbind(sample = sample_name, summ_arm_train) %>% mutate(chr = as.character(chr))

    standard_s = standard %>% filter(sample == sample_name)
    train_s = train_s %>% right_join(select(standard_s, sample, chr, arm, alteration), by = c("sample", "chr", "arm")) %>% filter(chr != "Y") %>% mutate(chr_status = ifelse(alteration == 0, "dipl", "no_dipl")) %>% metr_dipl()

    train_data = rbind(train_data, train_s)
  }

  write.table(train_data, file =  file.path(out_dir, "train_data"), sep = "\t", row.names = FALSE)
}
