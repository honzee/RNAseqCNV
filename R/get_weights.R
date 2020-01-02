#' Weight calculation
#' 
#' Get weights for CNA analysis
#' 
get_weights <- function(standard, count_dir, dipl_samp, metadata) {
  
  metadata <- fread(metadata, header = TRUE)
  
  # create baseline for normalization
  base_col <- data.frame(sampleName = dipl_samp)
  base_matr <- fread(file.path(count_dir, metadata$count[metadata$sample == dipl_samp[1]]), header = FALSE)
  genes <- pull(base_matr, 1)
  base_matr = base_matr %>% select(-1)
  for (i in 2:length(dipl_samp)) {
    count_s <- fread(file.path(count_dir, metadata$count[metadata$sample == dipl_samp[i]]), header = FALSE)[, 2]
    base_matr <- cbind(base_matr, count_s)
  }
  base_matr <- base_matr %>% as.data.frame() %>% magrittr::set_rownames(genes)
  
  # normalize all samples from the standard
  count_all <- c()
  for (i in 1:3) {
    
    if (! unique(standard$sample)[i] %in% dipl_samp) {
  
      # read count file for a sample from standard
      count_f_s <- metadata[metadata$sample == unique(standard$sample)[i], "count"]
      count_s <- fread(file.path(count_dir, count_f_s), header = FALSE) %>% as.data.frame() %>% pull(2)
      
      # add it to a baseline matrix
      final_mat <- cbind(count_s, base_matr)
      final_col <- rbind(data.frame(sampleName = unique(standard$sample)[i]), base_col)
      
      #create DESeqDataSet object
      ddsHTSeq = DESeqDataSetFromMatrix(countData = final_mat, colData = final_col, design= ~ 1)
      print("Loading count files is done!")
      
      #normalize
      ddsHTSeq <- estimateSizeFactors(ddsHTSeq)
      ddsCount <- counts(ddsHTSeq, normalized = TRUE)
      
      if (!is.null(count_all)) {
        count_all <- cbind(count_all, ddsCount[, 1])
        colnames(count_all)[ncol(count_all)] <- as.character(final_col$sampleName[1])
      } else {
        count_all <- data.frame(ddsCount)
        colnames(count_all) <- as.character(final_col$sampleName)
      }
    }
  }
  
  #annotate the normalized counts with chr and arm
  count_all_an <- count_all %>% mutate(ENSG = row.names(.)) %>% merge(refDataExp) %>% merge(centr_ref) %>% mutate(arm = ifelse(start < cstart, "p", ifelse(start > cend, "q", "centr"))) %>%
    select(-start, -end, -cstart, -cend)
  
  #calculate median of expression 
  count_dipl <- count_all_an[, which(colnames(count_all_an) %in% c(dipl_samp, "chr", "ENSG", "arm"))]
  med = apply(select(count_dipl, -chr, -ENSG), 1, median)
  #establish diploid expression baseline per chromosome
  ref_dipl = cbind(select(count_dipl, chr, arm, ENSG), median = as.numeric(med)) %>% group_by(chr, arm) %>% summarise(sum_chr_dipl = sum(median))
  
  #adjust non-diploid samples
  count_wide <- count_long %>% spread(key = sample, value = norm_exp)
  #choose non-diploid samples
  count_wide_non_dipl <- count_wide[, which(!colnames(count_wide) %in% dipl_samp)]
  #choose which columns to adjust
  to_adj <- which(! colnames(count_wide_non_dipl) %in% c("chr", "arm", "ENSG"))
  count_adj <- count_wide_non_dipl
  for (i in to_adj) {
    exp_sample <- count_wide_non_dipl[, c("chr", "arm", colnames(count_wide_non_dipl)[i])] 
    colnames(exp_sample)[3] <- "exp"
    #get sums of counts per chromosome and calculate the sife factor compared to diploid samples
    dipl <- filter(standard, sample == colnames(count_wide_non_dipl)[i], alteration == "0") %>% select(chr, arm)
    chrom_sum_samp_dipl <- exp_sample %>% group_by(chr, arm) %>% summarise(sum_chr_nd = sum(exp, na.rm = TRUE)) %>% ungroup() %>% merge(ref_dipl, by = c("chr", "arm")) %>% mutate(sizefac = sum_chr_dipl / sum_chr_nd) %>% 
      right_join(dipl)
    #calculate final sample size factor
    final_sf <- median(chrom_sum_samp_dipl$sizefac)
    #adjust the sample
    count_adj[, i] <- count_adj[, i] * final_sf
  }
  
  #merge dipl and non-dipl
  count_all_adj <- merge(count_adj, count_dipl, by = c("chr", "arm", "ENSG"))
  
  #annotate the data with alteration information
  count_all_adj_an <- gather(count_all_adj, key = "sample", value = "exp", which(!colnames(count_all_adj) %in% c("chr", "arm", "ENSG"))) %>% right_join(select(standard, sample, chr, arm, alteration))
 
  #get per-gene CNA expression correlation
  genes <- unique(count_all_adj_an$ENSG)
  gene_cor <- c()
  for (i in 1:length(genes)) {
    gene_data <- count_all_adj_an %>% filter(ENSG == genes[i])
    cor <- cor.test(x = gene_data$exp, y = as.numeric(gene_data$alteration))
    
    if (!is.null(gene_cor)) {
      gene_cor <- rbind(gene_cor, data.frame(ENSG = genes[i], pearson_r = cor$estimate, p = cor$p.value))
    } else {
      gene_cor <- data.frame(ENSG = genes[i], pearson_r = cor$estimate, p = cor$p.value)
    }
  }
  
  #get per-gene variance
  count_var <- count_all_adj %>% mutate(var = apply(select(., -ENSG, -chr, -arm), 1, function(x) var(x, na.rm = TRUE))) %>% arrange(var)
  
  #calculate weights
  weight_table <- select(count_var, ENSG, var) %>% right_join(gene_cor, by = "ENSG") %>% mutate(pearson_r = ifelse(pearson_r < 0.05, 0, pearson_r), pearson_r_q = pearson_r^2,   var_div_q = 1/var^2, weight = pearson_r_q*var_div_q) %>%
    arrange(desc(weight)) %>% select(ENSG, weight)
  
  return(weight_table)
}
 