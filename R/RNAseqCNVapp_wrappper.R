#' RNAseqCNA_wrapper
#'
#' Wrapper for generating figures for analysis and preview figure
#'
#' @param config R script assigning variables needed for the analysis
#' @param metadata path to a metadata table with three columns. First colum represents sample names, second file names of count files, third file names of snv files.
#' @export RNAseqCNV_wrapper
RNAseqCNV_wrapper <- function(config, metadata, adjust = TRUE, arm_lvl = TRUE, estimate = TRUE, referData = refDataExp, keptSNP = keepSNP, par_region = par_reg, centr_refer = centr_ref, weight_tab = weight_table, model_gend = model_gender, model_dip = model_dipl, model_alter = model_alt, chroms = chrs, base_matrix = base_matr, base_column = base_col, scale_cols = scaleCols, dpRatioChromEdge = dpRatioChrEdge, minDepth=20, minReadCnt = 30, q = 0.9) {

  #Check the format

  #Create sample table
  source(config)

  sample_table = fread(metadata, header = FALSE)
  HTSeq_f = pull(sample_table, 2)
  snv_f = pull(sample_table, 3)
  #create paths to the files
  sample_table$count_path <- file.path(count_dir, HTSeq_f)
  sample_table$snv_path <-  file.path(snv_dir, snv_f)

  #Create estimation table
  est_table <- data.frame(sample = character(),
                          gender = factor(levels = c("female", "male")),
                          chrom_n = integer(),
                          type = factor(levels = c("diploid", "high hyperdiploid", "low hyperdiploid", "low hypodiploid", "near haploid", "high hypodiploid", "near diploid")),
                          alterations = character(), stringsAsFactors = FALSE)

  #Create a log file to keep track of the analysis
  log <-  file.path(out_dir, "log.txt")
  cat("", file = log)

  #Run the analysis for every sample in the table
  for(i in 1:nrow(sample_table)) {

    sample_name <- as.character(sample_table[i, 1])

    #load SNP data
    smpSNP <- prepare_snv(sample_table = sample_table, sample_num = i, centr_ref = centr_ref, minDepth = minDepth, chrs = chroms)

    #calculate vst values with DESeq2
    vst <- get_vst2(sample_table = sample_table, minReadCnt = minReadCnt, q = q, sample_num = i, base_col = base_column, base_matr = base_matrix, weight_table = weight_tab, keep_perc = 0.8)

    #calculate medians for analyzed genes
    pickGeneDFall <- get_med(vst = vst, refDataExp = referData)

    #filter SNP data base on dpSNP database
    smpSNPdata.tmp <- filter_snv(smpSNP[[1]], keepSNP = keptSNP)

    #analyze chromosome-level metrics (out-dated)
    smpSNPdata <- calc_chrom_lvl(smpSNPdata.tmp)

    #arm-level metrics
    smpSNPdata_a_2 <- calc_arm(smpSNPdata.tmp)

    s_vst <- select(vst, !!quo(sample_name)) %>% mutate(ENSG = rownames(vst))

    #join reference data and weight data
    s_vst <- vst_norm(s_vst = s_vst, pickGeneDFall, refDataExp = referData, weight_table)

    #remove PAR regions
    s_vst <- remove_par(s_vst = s_vst, par_reg = par_region)

    #Calculate metrics for chromosome arms
    feat_tab <- get_arm_metr(s_vst = s_vst, smpSNPdata = smpSNPdata_a_2, sample_name = sample_names, centr_ref = centr_ref)

    #estimate gender
    s_vst_gend <- s_vst %>% filter(ENSG %in% c("ENSG00000114374", "ENSG00000012817", "ENSG00000260197", "ENSG00000183878")) %>%  select(ENSG, !!quo(sample_name)) %>% spread(key = ENSG, value = !!quo(sample_name))
    gender = ifelse(randomForest:::predict.randomForest(model_gend, newdata = s_vst_gend, type = "class") == 1, "male", "female")

    #preprocess data for karyotype estimation and diploid level adjustement
    if (adjust == TRUE | estimate == TRUE) {

      # model diploid level
      feat_tab$chr_status <- randomForest:::predict.randomForest(model_dip, feat_tab, type = "class")

      #exclude non-informative regions
      feat_tab_alt <- feat_tab %>%
        filter(arm != "p" | !chr %in% c(13, 14, 15, 21)) %>%
        metr_dipl()

      #model alteration on chromosome arms
      feat_tab_alt <- feat_tab_alt %>% mutate(alteration = as.character(randomForest:::predict.randomForest(model_alter, ., type = "class")))
      feat_tab_alt$alteration_prob <- apply(randomForest:::predict.randomForest(model_alter, feat_tab_alt, type = "prob"), 1, max)

      feat_tab_alt <- colour_code(feat_tab_alt) %>% group_by(chr) %>% mutate(alteration = as.character(alteration), chr_alt = as.character(ifelse(length(unique(alteration)) == 1, unique(alteration), "ab")))

    }

      #estimate karyotype
      if (estimate == TRUE) {

        kar_list <- gen_kar_list(feat_tab_alt = feat_tab_alt, sample_name = sample_name, gender = gender)

        est_table <- rbind(est_table, kar_list)
        write.table(x = est_table, file = file.path(out_dir, "estimation_table.tsv"), sep = "\t")
        write.table(x = cbind(est_table , status = "not checked", comments = "none"), file = file.path(out_dir, "manual_an_table.tsv"), sep = "\t")


      }

      #adjust for diploid level
      if (adjust == TRUE) {
        s_vst <-  adjust_dipl(feat_tab_alt, s_vst)
      }

      #calculate box plots for plotting
      box_wdt <- get_box_wdt(s_vst = s_vst, chrs = chroms, scaleCols = scale_cols)

      #adjust y axis limits
      ylim <- adjust_ylim(box_wdt = box_wdt, ylim = c(-0.4, 0.4))


      s_vst_final <- prep_expr(s_vst = s_vst, dpRatioChrEdge = dpRatioChromEdge, ylim = ylim, chrs = chroms)

      if(arm_lvl == TRUE) {

        centr_res <- rescale_centr(centr_ref, s_vst_final)

      }

      s_vst_final <- filter_expr(s_vst_final = s_vst_final, cutoff = 0.6)

      #plot arm-level figures
      if(arm_lvl == TRUE) {

          chr_dir = file.path(out_dir, sample_name)
          dir.create(path = chr_dir)
          chr_to_plot <- c(1:22, "X")

      }

      #plot every chromosome
      for (i in chr_to_plot) {

        gg_exp_zoom <- plot_exp_zoom(s_vst_final = s_vst_final, centr_res = centr_res, plot_chr = i,  estimate = estimate, feat_tab_alt = feat_tab_alt)

        yAxisMax_arm = get_yAxisMax(smpSNPdata = smpSNPdata, plot_chr = i)

        gg_snv_arm_p <- plot_snv_arm(smpSNPdata_a = smpSNPdata_a_2, plot_arm = "p", plot_chr = i, yAxisMax = yAxisMax_arm)

        gg_snv_arm_q <- plot_snv_arm(smpSNPdata_a = smpSNPdata_a_2, plot_arm = "q", plot_chr = i, yAxisMax = yAxisMax_arm)

        gg_arm <- chr_plot(p_snv = gg_snv_arm_p, q_snv = gg_snv_arm_q, arm_expr = gg_exp_zoom)

        ggsave(filename = file.path(chr_dir, paste0("chromosome_", i, ".png")), plot = gg_arm, device = "png", width = 20, height = 10)

      }

      gg_exp <- plot_exp(s_vst_final = s_vst_final, box_wdt = box_wdt, sample_name = sample_name, ylim = ylim, estimate = estimate, feat_tab_alt = feat_tab_alt)

      gg_snv <- plot_snv(smpSNPdata, chrs = chroms, sample_name = sample_name)

      fig <- arrange_plots(gg_exp = gg_exp, gg_snv = gg_snv)

      ggsave(plot = fig, filename = file.path(out_dir, paste0(sample_name, "_CNA_fig.png")), device = 'png', width = 16, height = 10)

  }
}


