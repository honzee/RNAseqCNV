# Wrapper for generating figures for analysis and preview figure
gen_fig_wrapper <- function(config, metadata, avail, sample_table, to_analyse, adjust, arm_lvl, estimate, refData, keepSNP, par_reg, centr_ref, weight_table, model_gender, model_dipl, model_alt, chrs, base_matr, base_col, scaleCols, dpRatioChrEdge, minDepth=20, minReadCnt = 3, samp_prop = 0.8) {

    #Is any neccessary input missing?
    if (all(metadata == "no_input")) {showNotification("A metadata file is needed to generate figures", duration = 5, id = "not_conf", type = "message"); return(NULL) }
    if (all(config == "no_input")) {showNotification("A config file is needed to generate figures", duration = 5, id = "not_count", type = "message"); return(NULL)}

    #Is the directory relevant?
    if (config["count_dir"] == FALSE) {showNotification("Could not find directory with count files", duration = 5, id = "not_valid_count", type = "message"); return(NULL)}
    if (config["snv_dir"] == FALSE) {showNotification("Could not find directory with snv files", duration = 5, id = "not_valid_snv", type = "message"); return(NULL)}
    if (config["out_dir"] == FALSE) {showNotification("Could not find the output directory", duration = 5, id = "not_out", type = "message"); return(NULL)}
    if (all(metadata == "incorrect_format")) {showNotification("Config file does not have the neccessary three columns", duration = 5, id = "not_valid_config", type = "message"); return(NULL)}

    if (avail != "all_present") {showNotification(avail , duration = NULL, id = "avail", type = "warning"); return(NULL)}


      #Create a table to write the estimation into
      est_table <- data.frame(sample = character(),
                              gender = factor(levels = c("female", "male")),
                              chrom_n = integer(),
                              type = factor(levels = c("diploid", "high hyperdiploid", "low hyperdiploid", "low hypodiploid", "near haploid", "high hypodiploid", "near diploid")),
                              alterations = character(), stringsAsFactors = FALSE)

      #Create a log file to keep track of the analysis
      log <-  paste0(config["out_dir"], "/", "log.txt")
      cat("", file = log)


    #Run the code with progress bar
    withProgress(message = "Analyzing..", value = 0, {

      # set the number of samples to be analysed
      samples <- 1:to_analyse

      for(i in samples) {

        sample_name <- as.character(sample_table[i, 1])

        incProgress(amount = 1/samples, message = paste0("Analyzing sample: ", sample_name))

        #Run with secondary progress indicator
        withProgress(message = "Currently:", value = 0, {

          incProgress(amount = 0.15, detail = "Prepring SNV and count data")

          #load SNP data
          smpSNP <- prepare_snv(sample_table = sample_table, sample_num = i, centr_ref = centr_ref, minDepth = minDepth, chrs = chrs)

          #check whether the snv file is in correct format
          if(is.null(smpSNP)) {

            writeLines(c(readLines(log), paste("File", sample_table$snv_path[i], "has an unsupported format"), paste("Sample", sample_name, "skipped")), con = log)
            incProgress(amount = 1/nrow(sample_table), message = "File skipped")
            next()

          }

          #calculate noemalized count values
          count_norm <- get_norm_exp(sample_table = sample_table, sample_num = i, diploid_standard = dipl_standard, minReadCnt = minReadCnt, samp_prop = samp_prop)

          #check whether the count file is in correct format
          if (is.null(count_norm)) {
              writeLines(c(readLines(log), paste("File", sample_table$count_path[i], "has an unsupported format"), paste("Sample", sample_name, "skipped")), con = log)
              incProgress(amount = 1/nrow(sample_table), message = "File skipped")
              next()
          }


          incProgress(amount = 0.5, detail = "Calculating expression medians for each gene")

          #calculate medians for analyzed genes
          pickGeneDFall <- get_med(count_norm = count_norm, refData = refData)

          #filter SNP data base on dpSNP database
          smpSNPdata.tmp <- filter_snv(smpSNP[[1]], keepSNP = keepSNP)

          #analyze chromosome-level metrics
          smpSNPdata <- calc_chrom_lvl(smpSNPdata.tmp)

          #analyze arm-level metrics
          smpSNPdata_a <- calc_arm_lvl(smpSNPdata.tmp)

          #arm-level metrics
          smpSNPdata_a_2 <- calc_arm(smpSNPdata.tmp)

          incProgress(amount = 0.1, detail = "Prepring expression data")

          count_ns <- select(count_norm, !!quo(sample_name)) %>% mutate(ENSG = rownames(count_norm))

          #join reference data and weight datacount_trans
          count_ns <- count_transform(count_ns = count_ns, pickGeneDFall, refData, weight_table)

          #remove PAR regions
          count_ns <- remove_par(count_ns = count_ns, par_reg = par_reg)

          feat_tab <- get_arm_metr(count_ns = count_ns, smpSNPdata = smpSNPdata_a_2, sample_name = sample_names, centr_ref = centr_ref, chrs = chrs)

          #estimate gender
          count_ns_gend <- count_ns %>% filter(ENSG %in% c("ENSG00000114374", "ENSG00000012817", "ENSG00000260197", "ENSG00000183878")) %>%  select(ENSG, !!quo(sample_name)) %>% spread(key = ENSG, value = !!quo(sample_name))
          gender = ifelse(randomForest:::predict.randomForest(model_gender, newdata = count_ns_gend, type = "class") == 1, "male", "female")

          #preprocess data for karyotype estimation and diploid level adjustement
          if (adjust == TRUE | estimate == TRUE) {

            incProgress(amount = 0.05, detail = "Estimating diploid baseline")

            feat_tab$chr_status <- randomForest:::predict.randomForest(model_dipl, feat_tab, type = "class")

            feat_tab_alt <- feat_tab %>%
              filter(arm != "p" | !chr %in% c(13, 14, 15, 21)) %>%
              metr_dipl()

            feat_tab_alt <- feat_tab_alt %>% mutate(alteration = as.character(randomForest:::predict.randomForest(model_alt, ., type = "class")))
            feat_tab_alt$alteration_prob <- apply(randomForest:::predict.randomForest(model_alt, feat_tab_alt, type = "prob"), 1, max)

            feat_tab_alt <- colour_code(feat_tab_alt) %>% group_by(chr) %>% mutate(alteration = as.character(alteration), chr_alt = as.character(ifelse(length(unique(alteration)) == 1, unique(alteration), "ab")))

          } else {
            incProgress(amount = 0.05)
          }


          #estimate karyotype
          if (estimate == TRUE) {

            incProgress(amount = 0.05, detail = "Estimating karyotype")

            kar_list <- gen_kar_list(feat_tab_alt = feat_tab_alt, sample_name = sample_name, gender = gender)

            est_table <- rbind(est_table, kar_list)
            write.table(x = est_table, file = paste0(config["out_dir"], "/", "estimation_table.tsv"), sep = "\t", quote = FALSE)
            write.table(x = cbind(est_table , status = "not checked", comments = "none"), file = paste0(config["out_dir"], "/", "manual_an_table.tsv"), sep = "\t", quote = FALSE)

          }

          #adjust for diploid level
          if (adjust == TRUE) {
            incProgress(amount = 0.05, detail = "Adjusting for diploid level")
            count_ns <-  adjust_dipl(feat_tab_alt, count_ns)
          } else {
            incProgress(amount = 0.05)
          }

          #calculate box plots for plotting
          box_wdt <- get_box_wdt(count_ns = count_ns, chrs = chrs, scaleCols = scaleCols)

          #adjust y axis limits
          ylim <- adjust_ylim(box_wdt = box_wdt, ylim = c(-0.4, 0.4))

          count_ns_final <- prep_expr(count_ns = count_ns, dpRatioChrEdge = dpRatioChrEdge, ylim = ylim, chrs = chrs)

          if(arm_lvl == TRUE) {

            centr_res <- rescale_centr(centr_ref, count_ns_final)

          }

          count_ns_final <- filter_expr(count_ns_final = count_ns_final, cutoff = 0.6)

          #plot arm-level figures
          if(arm_lvl == TRUE) {

            incProgress(amount = 0.08, detail = "Plotting chromosomes in detail")

            chr_dir <- paste0(config["out_dir"], "/", sample_name)
            dir.create(path = chr_dir)
            chr_to_plot <- c(1:22, "X")


            for (i in chr_to_plot) {

              gg_exp_zoom <- plot_exp_zoom(count_ns_final = count_ns_final, centr_res = centr_res, plot_chr = i,  estimate = estimate, feat_tab_alt = feat_tab_alt)

              yAxisMax_arm = get_yAxisMax(smpSNPdata = smpSNPdata_a, plot_chr = i)

              gg_snv_arm_p <- plot_snv_arm(smpSNPdata_a = smpSNPdata_a, plot_arm = "p", plot_chr = i, yAxisMax = yAxisMax_arm)

              gg_snv_arm_q <- plot_snv_arm(smpSNPdata_a = smpSNPdata_a, plot_arm = "q", plot_chr = i, yAxisMax = yAxisMax_arm)

              gg_arm <- chr_plot(p_snv = gg_snv_arm_p, q_snv = gg_snv_arm_q, arm_expr = gg_exp_zoom)

              ggsave(filename = paste0(config["out_dir"], "/", sample_name, "/", "chromosome_", i, ".png"), plot = gg_arm, device = "png", width = 20, height = 10, dpi = 100)

            }

          }


          incProgress(amount = 0.02, detail = "Plotting main figure")

          gg_exp <- plot_exp(count_ns_final = count_ns_final, box_wdt = box_wdt, sample_name = sample_name, ylim = ylim, estimate = estimate, feat_tab_alt = feat_tab_alt, gender = gender)

          gg_snv <- plot_snv(smpSNPdata, chrs = chrs, sample_name = sample_name)

          fig <- arrange_plots(gg_exp = gg_exp, gg_snv = gg_snv)

          ggsave(plot = fig, filename = paste0(chr_dir, "/", sample_name, "_CNV_main_fig.png"), device = 'png', width = 16, height = 10, dpi = 200)

          writeLines(c(readLines(log), paste("Sample", sample_name, "analyzed successfully")), con = log)

        })
      }
    })

      showNotification("Analysis complete", id = "an_compl", type = "message", closeButton = TRUE)

      if (to_analyse == 1) {

        chr_choices_prev <- reactive({

          chromosomes <- list.files(path = paste0(config["out_dir"], "/", sample_name, "/"), pattern = "^chromosome_.*.png$")

          if (length(chromosomes) == 0) {
              return(NULL)
            } else if (length(chromosomes) > 0) {
              choices <-  factor(gsub("_|.png", " ", chromosomes), levels = paste0("chromosome ", c(1:22, "X"), " "))
              choices <- choices[order(choices)]
              return(choices)
          }

        })

      }
}

