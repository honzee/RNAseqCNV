# Wrapper for generating figures for analysis and preview figure
gen_fig_wrapper <- function(config, metadata, snv_format, avail, sample_table, to_analyse, adjust, arm_lvl, estimate_lab, refData, keepSNP, par_reg, centr_ref, weight_table, generate_weights, model_gender, model_dipl, model_alt, model_noSNV, chrs, batch, standard_samples, scaleCols, dpRatioChrEdge, minDepth=20, minReadCnt = 3, samp_prop = 0.8, weight_samp_prop = 1) {


    #Is any neccessary input missing?
    if (all(metadata == "no_input")) {showNotification("A metadata file is needed to generate figures", duration = 5, id = "not_conf", type = "message"); return(NULL) }
    if (all(config == "no_input")) {showNotification("A config file is needed to generate figures", duration = 5, id = "not_count", type = "message"); return(NULL)}

    #Is the directory relevant?
    if (config["count_dir"] == FALSE) {showNotification("Could not find directory with count files", duration = 5, id = "not_valid_count", type = "message"); return(NULL)}
    if (config["snv_dir"] == FALSE) {showNotification("Could not find directory with snv files", duration = 5, id = "not_valid_snv", type = "message"); return(NULL)}
    if (config["out_dir"] == FALSE) {showNotification("Could not find the output directory", duration = 5, id = "not_out", type = "message"); return(NULL)}
    if (all(metadata == "incorrect_format")) {showNotification("Config file does not have the neccessary three columns", duration = 5, id = "not_valid_config", type = "message"); return(NULL)}

    if (avail != "all_present") {showNotification(avail , duration = NULL, id = "avail", type = "warning"); return(NULL)}

    #check whether a snv_format was selected
    if (is.null(snv_format) | !snv_format %in% c("vcf", "custom")) {
      stop("snv_format parameter has to be either 'vcf' or 'custom'")
    }

      #Create a table to write the estimation into
      est_table <- data.frame(sample = character(),
                              gender = factor(levels = c("female", "male")),
                              chrom_n = integer(),
                              alterations = character(), stringsAsFactors = FALSE)

    #Run the code with progress bar
    withProgress(message = "Analyzing..", value = 0, {

      #If batch analysis was selected normalize the input samples together
      if(batch == TRUE) {
        incProgress(detail = "Normalizing samples as a batch..")

        #calculate normalized count values with DESeq2 normalization method for batch of samples from the input
        count_norm <- get_norm_exp(sample_table = sample_table, sample_num = 1, standard_samples = standard_samples, minReadCnt = minReadCnt, samp_prop = samp_prop, weight_table = weight_table, weight_samp_prop = weight_samp_prop, batch = batch, generate_weights = generate_weights)

        #calculate median gene expression across diploid reference and analyzed sample for batch of samples from the input
        pickGeneDFall <- get_med(count_norm = count_norm, refDataExp = refDataExp, generate_weights = generate_weights)

        if (generate_weights == TRUE) {
          #create weights based on variance of gene expression and expression depth of the batch of samples
          weight_table <- create_weights(pickGeneDFall)
        }

        incProgress(detail = "")
      }

      # set the number of samples to be analysed
      samples <- 1:to_analyse

      for(i in samples) {

        sample_name <- as.character(sample_table[i, 1])

        incProgress(amount = 1/samples, message = paste0("Analyzing sample: ", sample_name))

        #Run with secondary progress indicator
        withProgress(message = "Currently:", value = 0, {

          #normalize the samples with in-build standard or standard from the input and calculate gene medians
          if(batch == FALSE) {

            incProgress(amount = 0.5, detail = "Normalizing gene expression and calculating medians for each gene")

            #calculate normalized count values with DESeq2 normalization method
            count_norm <- get_norm_exp(sample_table = sample_table, sample_num = i, standard_samples = standard_samples, minReadCnt = minReadCnt, samp_prop = samp_prop, weight_table = weight_table, weight_samp_prop = weight_samp_prop, batch = batch, generate_weights = generate_weights)

            #calculate median gene expression across diploid reference and analyzed sample
            pickGeneDFall <- get_med(count_norm = count_norm, refDataExp = refDataExp, generate_weights = generate_weights)

            if (generate_weights == TRUE) {
              #create weights based on variance of gene expression and expression depth of the batch of samples
              weight_table <- create_weights(pickGeneDFall)
            }

          } else {
            incProgress(amount = 0.5)
          }

          # if count file format was incorrect print out a message and skip this sample
          if (is.character(count_norm)) {
            showNotification(count_norm, duration = NULL, id = "not_valid_count_file", type = "message")
            next()
          }

          incProgress(amount = 0.15, detail = "Processing SNV and count data")

          #load SNP data
          smpSNP <- prepare_snv(sample_table = sample_table, sample_num = i, centr_ref = centr_ref, snv_format = snv_format, minDepth = minDepth)

          if (is.character(smpSNP[[1]])) {
            showNotification(paste0(sample_name, ": ", smpSNP[[1]]), duration = NULL, id = "not_valid_snv_info", type = "message")
            next()
          }

          #filter SNP data base on dpSNP database
          smpSNPdata.tmp <- filter_snv(smpSNP[[1]], keepSNP = keepSNP, minDepth = minDepth)

          #analyze chromosome-level metrics
          smpSNPdata <- calc_chrom_lvl(smpSNPdata.tmp)

          #analyze arm-level metrics
          smpSNPdata_a <- calc_arm_lvl(smpSNPdata.tmp)

          #arm-level metrics
          smpSNPdata_a_2 <- calc_arm(smpSNPdata.tmp)

          incProgress(amount = 0.1, detail = "Analysing expression data")

          #select sample
          count_norm_samp <- count_norm %>% select(!!quo(sample_name)) %>% mutate(ENSG = rownames(.))

          #join reference data and weight datacount_trans
          count_ns <- count_transform(count_ns = count_norm_samp, pickGeneDFall, refData, weight_table)

          #remove PAR regions
          count_ns <- remove_par(count_ns = count_ns, par_reg = par_reg)

          feat_tab <- get_arm_metr(count_ns = count_ns, smpSNPdata = smpSNPdata_a_2, sample_name = sample_names, centr_ref = centr_ref)

          #estimate gender
          count_ns_gend <- count_norm_samp %>% filter(ENSG %in% "ENSG00000012817") %>%  select(ENSG, !!quo(sample_name)) %>% spread(key = ENSG, value = !!quo(sample_name))
          gender <- ifelse(predict(model_gender, newdata = count_ns_gend, type = "response") > 0.5, "male", "female")

          #preprocess data for karyotype estimation and diploid level adjustement
            incProgress(amount = 0.05, detail = "Estimating diploid baseline")

            feat_tab$chr_status <- randomForest:::predict.randomForest(model_dipl, feat_tab, type = "class")

            #exclude non-informative regions  and
            #if the model was not able to call changes
            #(mainly due problematic density graphs on chromosome X) change the value to unknown
            feat_tab_dipl <- feat_tab %>%
              filter(arm != "p" | !chr %in% c(13, 14, 15, 21)) %>% mutate(chr_status = ifelse(is.na(chr_status), "unknown", as.character(chr_status))) %>%
              metr_dipl()

            #model alteration on chromosome arms an in case of problematic SNV graph, use model without this information included
            print(paste("Estimating chromosome arm CNV",":", sample_name))
            feat_tab_alt <- feat_tab_dipl %>% filter(chr_status != "unknown") %>% mutate(alteration = as.character(randomForest:::predict.randomForest(model_alt, ., type = "class")),
                                                                                         alteration_prob = apply(randomForest:::predict.randomForest(model_alt, ., type = "prob"), 1, max))
            if (any(feat_tab_dipl$chr_status == "unknown")) {
              feat_tab_alt <- feat_tab_dipl %>% filter(chr_status == "unknown") %>% mutate(alteration = as.character(randomForest:::predict.randomForest(model_noSNV, ., type = "class")),
                                                                                           alteration_prob = apply(randomForest:::predict.randomForest(model_noSNV, ., type = "prob"), 1, max)) %>%
                bind_rows(feat_tab_alt)
            }
            feat_tab_alt <- colour_code(feat_tab_alt, conf_tresh = 0.85) %>% group_by(chr) %>% mutate(alteration = as.character(alteration), chr_alt = as.character(ifelse(length(unique(alteration)) == 1, unique(alteration), "ab")))

          #estimate karyotype
          incProgress(amount = 0.05, detail = "Estimating karyotype")

            kar_list <- gen_kar_list(feat_tab_alt = feat_tab_alt, sample_name = sample_name, gender = gender)

            est_table <- rbind(est_table, kar_list)
            write.table(x = est_table, file = paste0(config["out_dir"], "/", "estimation_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
            write.table(x = cbind(est_table , status = "not checked", comments = "none"), file = paste0(config["out_dir"], "/", "manual_an_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

          #adjust for diploid level
          if (adjust == TRUE) {
            incProgress(amount = 0.05, detail = "Adjusting for diploid level")
            count_ns <-  adjust_dipl(feat_tab_alt, count_ns)
          } else {
            incProgress(amount = 0.05)
          }

          #calculate box plots for plotting
          box_wdt <- get_box_wdt(count_ns = count_ns, scaleCols = scaleCols)

          #adjust y axis limits
          ylim <- adjust_ylim(box_wdt = box_wdt, ylim = c(-0.4, 0.4))

          count_ns_final <- prep_expr(count_ns = count_ns, dpRatioChrEdge = dpRatioChrEdge, ylim = ylim)

          if(arm_lvl == TRUE) {

            centr_res <- rescale_centr(centr_ref, count_ns_final)

          }

          count_ns_final <- filter_expr(count_ns_final = count_ns_final, cutoff = 0.6)

          #Create per-sample folder for figures
          chr_dir = file.path(config["out_dir"], sample_name)
          dir.create(path = chr_dir)

          #Create and plot arm-level figure
          incProgress(amount = 0.02, detail = "Plotting main figure")

          gg_exp <- plot_exp(count_ns_final = count_ns_final, box_wdt = box_wdt, sample_name = sample_name, ylim = ylim, estimate = estimate_lab, feat_tab_alt = feat_tab_alt, gender = gender)

          gg_snv <- plot_snv(smpSNPdata, sample_name = sample_name, estimate = estimate_lab)

          fig <- arrange_plots(gg_exp = gg_exp, gg_snv = gg_snv)

          print(paste0("Plotting main figure: ", sample_name))

          ggsave(plot = fig, filename = paste0(chr_dir, "/", sample_name, "_CNV_main_fig.png"), device = 'png', width = 16, height = 10, dpi = 200)


          #plot arm-level figures
          if(arm_lvl == TRUE) {
            incProgress(amount = 0.05, detail = "Plotting arm-level figures")

            chr_to_plot <- c(1:22, "X")

            centr_res <- rescale_centr(centr_ref, count_ns_final)

            print(paste0("Plotting arm-level figures: ", sample_name))

            #plot every chromosome
            for (i in chr_to_plot) {

              print(paste0("Plotting chr ", i, " arm-level figure"))

              gg_exp_zoom <- plot_exp_zoom(count_ns_final = count_ns_final, centr_res = centr_res, plot_chr = i,  estimate = estimate_lab, feat_tab_alt = feat_tab_alt)

              yAxisMax_arm = get_yAxisMax_arm(smpSNPdata = smpSNPdata_a_2, plot_chr = i)

              gg_snv_arm_p <- plot_snv_arm(smpSNPdata_a = smpSNPdata_a_2, plot_arm = "p", plot_chr = i, yAxisMax = yAxisMax_arm)

              gg_snv_arm_q <- plot_snv_arm(smpSNPdata_a = smpSNPdata_a_2, plot_arm = "q", plot_chr = i, yAxisMax = yAxisMax_arm)

              gg_arm <- chr_plot(p_snv = gg_snv_arm_p, q_snv = gg_snv_arm_q, arm_expr = gg_exp_zoom)

              ggsave(filename = file.path(chr_dir, paste0("chromosome_", i, ".png")), plot = gg_arm, device = "png", width = 20, height = 10, dpi = 100)

            }

          } else {
            incProgress(amount = 0.05)
          }

          print(paste0("Analysis for sample: ", sample_name, " finished"))

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

