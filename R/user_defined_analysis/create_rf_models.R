#' Create model for CNA estimation
#'
#' Create random forest model for diploid level estimation and alteration estimation
#'
#' @param train_data from get_train_dat function
#'
#'
create_rf_models <- function(train_data){
  # Read training data
  train_d <- read.table(file.path(out_dir, "train_data"), header = TRUE)

  # model for diploid/non-diploid estimation
  model_dipl <- randomForest(chr_status ~ chr + peak_max + peak_m_dist + peakdist + sd + sds_median + sds_025 + sds_075 + n_02_04 + n_04, data = train_d_over, importance = TRUE)
  train_d_over_alt_all <- train_d_over %>% mutate(alteration = as.factor(alteration)) %>% metr_dipl()

  # model for alteration estimation
  model_alt <- randomForest(alteration ~ chr + peak_max + peak_m_dist + peakdist + sds_median + sds_025 + sd_dipl + n_02_04 + n_04 + mean_dipl + sds_median_dipl + sds_025_dipl + sds_075_dipl, train_d_over_alt_all, type = "class")

  return(list(model_dipl, model_alt))
}
