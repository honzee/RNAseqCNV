#' Weight calculation
#'
#' Get weights for CNA analysis
#'
#' The function provides an option to create customized weight matrix based on specific data (cancer) type. The weights are based on the correlation of gene expression with CNV and expression variance across the cohort. To calculate gene expression-CNV correlation
#' reliable set of CNVs on arm level is required.  For the function to work properly, each chromosome arm should have
#' minimum two CNV states in the provided CNV data. If this is not the case the parameter fill_cor should be set as TRUE to artificially fill out the correlation coefficients.
#'
#' @param CNV_data path to a table with known CNVs on arm-level. The information is used to generate CNV-gene expression level correlation coeffiecient, that are
#' subsequently used for gene weight calculation. The more samples are provided, the more robust the calculated correlation coeffiecient will be.
#' There are two optional input formats. Either in wide format, where each column represents one chromosome arm e.g. columns (1p,1q,2p,2q,3p,3q..) and
#' one additional column with sample name. The values in chromosome arm column are the CNV states on that chromosome arm, if the CNV state is unknown the value should be set to NA. Also the CNV_data file
#' can be in "long" format. Here the table is required to have three columns: sample: sample name, arm: chromosome arm (format: "1p", "19q", "Xp") and CNV: (CNV state, e.g. +1, -1, +2, 0)
#' @param standard_samp character vector with standard sample names, which will be used for normalization, and adjustment of diploid level. These samples should not
#' carry any CNVs. Providing standard samples is important for datasets, where high numbers of CNVs (e.g. 10 gained chromosomes or 10 deleted chromosomes) of one type are present is some samples.
#'  The sample names must be included in the first column of metadata table. If no standard samples samples are input, the data will be normalized as a batch and no diploid level adjustment will
#'  be performed.
#' @param metadata path to a metadata table with three columns. First colum: sample names, second column: file names of count files, third column: file names of snv files. There should be no header.
#'  More information is included in the package README file.
#' @param CNV_data_format string, in which format is the CNV_data provided, can be either "wide" or "long"
#' @param fill_cor if TRUE (default), genes for which correlation coefficients cannot be estimated and and those for which adjusted p-value by Benjamini-Hochberg method is higher than 0.1, the correlation coefficient value will be set
#' to 0.1. This was determined as the mean gene expression-CNV correlation coefficient from qcute lymphoblastic leukemia data on all genes. This parameter may help in creating weights for datasets, for which the coverage by CNVs
#' of genome is low and thus correlation coeffiecients may be constructed only for small proportion of genes. If FALSE, genes without determined correlation coefficients and those for which adjusted p-value is higher than 0.1 will be excluded.
#' @param output_file path to a file, where gene weight matrix will be saved.
get_weights <- function(CNV_data, count_dir, standard_samp = NULL, metadata, CNV_data_format, fill_cor = TRUE) {

  if(!CNV_data_format %in% c("long", "wide")) {
    stop("CNV_data_format has to be set to either long or wide format")
  }


  CNV_table <- fread(CNV_data)

  if(CNV_data_format == "wide") {
    CNV_table <- pivot_longer(data = CNV_table, names_to = "arm", values_to = "CNV", cols = c(paste0(c(1:22, "X"), "q"), paste0(c(1:12, 16:20, "X"), "p")))
  }


  #function for geometric
  gm_mean = function(a){prod(a)^(1/length(a))}

  #read metadata table
  metadata <- fread(metadata, header = FALSE)

  #get file paths to count files
  metadata$count_path <- file.path(count_dir, pull(metadata, 2))

  if (any(duplicated(metadata[, 2]))) {
    stop("Duplicated samples")
  }

  if (!is.null(standard_samples)) {

    #Check whether the samples are present in the sample table
    if (all(standard_samp %in% pull(metadata, 1)) == FALSE) {
      stop("The input standard samples are not in metadata table.")
    }
    #Create standard sample table
    standard_table <- create_standard(standard_samples = standard_samp, sample_table = metadata)
  }

  if (!is.null(standard_samp)) {
    for (i in 1:nrow(metadata)) {
      if (i == 1) {
        count_norm <- get_norm_exp(sample_table = metadata, standard_samples = standard_table, minReadCnt = 3, sample_num = i, samp_prop = 0.8, batch = FALSE, weight_samp_prop = NULL,  weight_table = NULL, generate_weights = TRUE) %>% .[, !duplicated(colnames(.))] %>% select(1) %>%
          mutate(ENSG = rownames(.))
      } else {
        count_norm <- inner_join(count_norm, get_norm_exp(sample_table = metadata, standard_samples = standard_table, minReadCnt = 3, sample_num = i, samp_prop = 0.8, batch = FALSE, weight_samp_prop = NULL,  weight_table = NULL, generate_weights = TRUE) %>% .[, !duplicated(colnames(.))] %>%
                                   select(1) %>%  mutate(ENSG = rownames(.)))
      }
    }

  } else {
    count_norm <- get_norm_exp(sample_table = metadata, standard_samples = standard_table, minReadCnt = 3, sample_num = i, samp_prop = 0.8, batch = TRUE, weight_samp_prop = NULL,  weight_table = NULL, generate_weights = TRUE) %>% mutate(ENSG = row.names(.))
    }

  #annotate the normalized counts with chromosome and chromosome arm
  count_norm_an <-  count_norm %>% left_join(as.data.frame(refDataExp), by = "ENSG") %>% left_join(mutate(centr_ref, chr = as.character(chr)), by = "chr") %>% mutate(arm = ifelse(start < cstart, "p", ifelse(start > cend, "q", "centr"))) %>%
    mutate(arm = paste0(chr, arm)) %>% select(-start, -end, -cstart, -cend, -chr)

  if (!is.null(standard_samp)) {

    #calculate pseudo reference in diploid samples for second normalization
    count_dipl_norm <-  count_norm_an[, which(colnames(count_norm_an) %in% c(standard_samp, "arm", "ENSG"))]
    pseudo_ref_norm <- apply(X = select(count_dipl_norm, -arm, - ENSG), MARGIN = 1, function(x) gm_mean(x))
    ref_dipl = cbind(select(count_dipl_norm, arm, ENSG), pseudo_ref = as.numeric(pseudo_ref_norm)) %>% filter(pseudo_ref != 0)

    #second normalization, especially for samples with high numbers of CNA
    count_norm_an_cor <- count_norm_an
    non_dipl <- pull(metadata[!V1 %in% standard_samp], V1)
    for (i in non_dipl) {
      sample_vst <- count_norm_an[, c(non_dipl, "ENSG", "arm")]
      colnames(sample_vst)[1] <- "count"
      size_fac = sample_vst %>% right_join(ref_dipl, by = c("ENSG", "arm")) %>% inner_join(filter(CNV_table, sample == i), by = c("arm")) %>% filter(CNV == 0) %>% mutate(size_facs = count/pseudo_ref) %>% summarise(median(size_facs))
      count_norm_an_cor[, i] <- count_norm_an_cor[, i] / size_fac
    }
    count_norm_an <- count_norm_an_cor
  }

  #annotate the data with alteration information
  count_final <- gather(count_norm_an, key = "sample", value = "vst", which(!colnames(count_norm_an) %in% c("arm", "ENSG"))) %>% inner_join(select(CNV_table, sample, arm, CNV), by = c("sample", "arm"))

  #get per-gene CNA expression correlation
  genes <- unique(count_final$ENSG)

  for (i in 1:length(genes)) {
    gene_data <- count_final %>% filter(ENSG == genes[i])
    suppressWarnings(cor <- cor.test(x = gene_data$vst, y = as.numeric(gene_data$CNV)))

    if (i != 1) {
      gene_cor <- rbind(gene_cor, data.frame(ENSG = genes[i], pearson_r = cor$estimate, p = cor$p.value))
    } else {
      gene_cor <- data.frame(ENSG = genes[i], pearson_r = cor$estimate, p = cor$p.value, stringsAsFactors = FALSE)
    }
  }

  # calculate adjusted p values, filter reliable genes and optionally fill out the rest of genes with base correleation coefficient of 0.1
  gene_cor_padj <- gene_cor %>% filter(!is.na(pearson_r) & !is.na(p)) %>% mutate(p_adj = p.adjust(p, method = "hochberg"))
  gene_cor_filt <- gene_cor_padj %>% filter(p_adj < 0.1)

  if (nrow(gene_cor_filt) < 1000 & fill_cor == FALSE) {
    stop("There is not enough data to generate robust correlation coefficients for most of the genes. Either expand the data or set fill_cor to TRUE to artificially assign a baseline correlation coefficient.")
  }

  if (fill_cor == TRUE) {
    gene_fill <- gene_cor %>% filter(is.na(pearson_r) | is.na(p)) %>% mutate(cor = 0.1)
    gene_fill <- gene_fill %>% bind_rows(gene_cor_padj %>% filter(p_adj  > 0.1) %>% mutate(cor = 0.1) %>% select(-p_adj))
  }

  #get per-gene variance and depth
  var <- apply(select(count_dipl_norm, -arm, -ENSG), 1, var)
  depth <- apply(select(count_dipl_norm, -arm, -ENSG), 1, mean)
  var_depth_table <- data.frame(ENSG = count_dipl_norm$ENSG, var = var, depth = depth, arm = count_dipl_norm$arm, stringsAsFactors = FALSE)

  #calculate weights
  weight_table <- var_depth_table %>% right_join(gene_cor, by = "ENSG") %>% mutate(weight = scales::rescale(pearson_r^5, to = c(1,100))*scales::rescale(1/var, to = c(1,100))) %>% mutate(chromosome_name = sub("p|q", "", arm)) %>%
    arrange(desc(weight)) %>% dplyr::select(ENSG, weight, chromosome_name)

  write.table(weight_table, file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")
}
