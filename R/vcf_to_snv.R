#' Fucnction that converts vcf files to tabular input for CNV analysis
vcf_to_snv <- function(vcf_file) {

  #read the vcf files
  vcf_data <- vcfR::read.vcfR(vcf_file)

  #select importand variables from fic
  snv_table <- as.data.frame(vcf_data@fix[, c("CHROM", "POS", "REF", "ALT", "QUAL")])

  #mapping quality
  message("Extracting mapping quality..")
  snv_table$MQ = sub("MQ=", "", str_extract(vcf_data@fix[, "INFO"], "MQ=[0-9]+\\.[0-9]*"))

  #reference allele and alternative allele depths
  message("Extracting reference allele and alternative allele depths..")
  gt_dat <- str_split(vcf_data@gt[, 2], ":")
  depth_vec <- sapply(gt_dat, function(x) x[2])
  depth_list <- str_split(depth_vec, ",")
  snv_table$refDP <- as.numeric(sapply(depth_list, function(x) x[1]))
  snv_table$varDP <- as.numeric(sapply(depth_list, function(x) x[2]))
  message("Needed information from vcf extracted")

  #calculate overall depth and mutant allele frequency
  snv_table <- snv_table %>% mutate(depth = refDP+varDP, maf = varDP/depth, start = POS, end = POS)
  snv_table <- snv_table[, c("CHROM", "start", "REF", "ALT", "end", "QUAL", "depth", "refDP", "varDP", "MQ", "maf")]
  colnames(snv_table) <- c("chr", "start", "ref", 'var', 'end', 'qual', 'depth', 'refDp', 'varDp', 'mapQ', "maf")

  #filtering base on maf and depth
  snv_table <- snv_table %>% filter(maf > 0.01, depth > 5)
  message("Conversion completed")
  return(snv_table)
}
