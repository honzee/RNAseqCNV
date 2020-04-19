#' Create model for CNA estimation
#'
#' Create random forest model for diploid level estimation and alteration estimation
#'
#' @param train_data path to training data created by get_train_data() function. However, the data should be checked for class imabalances and potentially oversamples/undersampled.
#' @param output_file path to a file, where the model will be saved as list in RDS format

#'
create_rf_models <- function(train_data){
  # Read training data
  train_d <- read.table(train_data, header = TRUE)

  # model for diploid/non-diploid estimation
  model_dipl <- randomForest(chr_status ~ n_02_04 + n_04 + sds_025 + chr + y_0.5 + sds_median + sd + peak_m_dist, data = train_d, type = "class")

  # model for alteration estimation
  model_alt <- randomForest(alteration ~ chr + sds_median_dipl + y_0.5 + sds_025_dipl + peak_m_dist + peak_max + mean_dipl + peakdist, data = train_d, type = "class")

  # model for alteration without SNV
  model_noSNV <- randomForest(alteration ~ sds_median_dipl + chr + sds_025_dipl + sd_dipl + mean_dipl, data = train_d, type = "class")

  saveRDS(model_alt, model_dipl, model_noSNV)
}
