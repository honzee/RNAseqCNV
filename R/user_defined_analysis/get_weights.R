#' Weight calculation
#'
#' Get weights for CNA analysis
#'
#' @param standard table with 4 compulsory columns: chr = chromosome, ENSG = ensembl gene id, sample = the sample name, alteration = CNA of that gene (+1 copy = 1, +2 = 2, -1 = -1, no CNV = 0, other = NA)
#' @param count_dir directory with count files
#' @param dipl_samp character vector with diploid samples
#' @param metadata same table as the input into RNAseqCNV_wrapper function
get_weights <- function(standard, count_dir, dipl_samp, metadata) {

  #function for geometric
  gm_mean = function(a){prod(a)^(1/length(a))}

  #read metadata table
  metadata <- fread(metadata, header = FALSE)

  if (any(duplicated(metadata[, 2]))) {
    stop("Duplicated samples")
  }

  #read in all count files
  count_table <- read.table(file = paste0(count_dir, metadata[1, 2]), header = FALSE, row.names = 1, stringsAsFactors = FALSE)
  for (i in 2:nrow(metadata)) {
    count_table <- cbind(count_table, read.table(file = paste0(count_dir, metadata[i, 2]), header = FALSE, row.names = 1, stringsAsFactors = FALSE))
  }
  colnames(count_table) <- pull(metadata, 1)
  count_table <- as.data.frame(count_table)

  ####filtering of samples###
  count_filt = as.data.frame(count_table) %>% mutate(keep_gene = apply(., MARGIN = 1, FUN = function(x) sum(x > 2) > (length(x) * 0.8)), ENSG = rownames(.), id = row_number()) %>% filter(keep_gene == TRUE) %>%
    select(-keep_gene, -id, -ENSG)

  ###normalize samples and calculate vst####
  count_norm <- count_filt
  count_dipl <- count_filt[colnames(count_filt) %in% dipl_samp]
  non_dipl <- pull(metadata, 1)[!pull(metadata, 1) %in% dipl_samp]
  for (i in 1:ncol(count_filt)) {
    to_trans <- count_filt %>% select(i) %>% bind_cols(count_dipl)
    count_matr <- as.matrix(to_trans)
    count_col <- as.data.frame(colnames(to_trans))

    dds <- DESeqDataSetFromMatrix(colData = count_col, countData = count_matr, design= ~ 1)
    dds_vst <- varianceStabilizingTransformation(dds, blind=T, fitType='local')
    count_norm[, i] <- assay(dds_vst)[, 1]
  }

  #annotate the normalized counts with chr
  count_norm_an <-  count_norm %>% mutate(ENSG = row.names(.)) %>% merge(refDataExp)

  #calculate pseudo reference in diploid samples for second normalization
  count_dipl_norm <-  count_norm_an[, which(colnames( count_norm_an) %in% c(dipl_samp, "chr", "ENSG"))]
  pseudo_ref_norm <- apply(X = select(count_dipl_norm, -chr, - ENSG), MARGIN = 1, function(x) gm_mean(x))
  ref_dipl = cbind(select(count_dipl_norm, chr, ENSG), pseudo_ref = as.numeric(pseudo_ref_norm)) %>% filter(pseudo_ref != 0)

  #second normalization, especially for samples with high numbers of CNA
  count_norm_an_cor <- count_norm_an
  for (i in non_dipl) {
    sample_vst <-  count_norm_an[, c(i, "chr", "ENSG")]
    colnames(sample_vst)[1] <- "count"
    size_fac = sample_vst %>% right_join(ref_dipl, by = c("ENSG", "chr")) %>% left_join(filter(standard, sample == i), by = c("ENSG", "chr")) %>% filter(alteration == 0) %>% mutate(size_facs = count/pseudo_ref) %>% summarise(median(size_facs))
    count_norm_an_cor[, i] <- count_norm_an_cor[, i] * size_fac
  }

  #annotate the data with alteration information
  count_final <- gather(count_norm_an_cor, key = "sample", value = "exp", which(!colnames(count_norm_an_cor) %in% c("chr", "ENSG"))) %>% left_join(select(standard, sample, ENSG, chr, alteration), by = c("ENSG", "sample", "chr"))

  #get per-gene CNA expression correlation
  genes <- unique(count_final$ENSG)
  gene_cor <- c()
  for (i in 1:10) {
    gene_data <- count_final %>% filter(ENSG == genes[i])
    cor <- cor.test(x = gene_data$exp, y = as.numeric(gene_data$alteration))

    if (!is.null(gene_cor)) {
      gene_cor <- rbind(gene_cor, data.frame(ENSG = genes[i], pearson_r = cor$estimate, p = cor$p.value))
    } else {
      gene_cor <- data.frame(ENSG = genes[i], pearson_r = cor$estimate, p = cor$p.value)
    }
  }

  #get per-gene variance and depth
  var <- apply(select(count_dipl_norm, -chr, -ENSG), 1, var)
  depth <- apply(select(count_dipl_norm, -chr, -ENSG), 1, mean)
  var_depth_table <- data.frame(ENSG = count_dipl_norm$ENSG, var = var, depth = depth, chr = count_dipl_norm$chr)

  #calculate weights
  wgt_tab <- var_depth_table %>% right_join(gene_cor, by = "ENSG") %>% mutate(weight = pearson_r^2*(1/var^2)*depth^2) %>% mutate(chromosome_name = chr) %>%
    arrange(desc(weight)) %>% dplyr::select(ENSG, weight, chromosome_name)

  return(weight_table)
}
