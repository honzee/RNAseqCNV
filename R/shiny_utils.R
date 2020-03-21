#' Importing packages
#' @import shiny
#' @import tidyr
#' @import dplyr
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import stringr
#' @importFrom data.table fread
#' @import ggpubr
#' @import shinyFiles

#### functions for deseq norm testing ####

# get normalized counts
get_norm_exp <- function(sample_table, sample_num, diploid_standard, minReadCnt, samp_prop, weight_table, weight_samp_prop) {

  #extract file names from sample table
  count_file <- pull(sample_table, "count_path")[sample_num]
  sample_n <- pull(sample_table, 1)[sample_num]

  #read the count table
  count_table <- fread(file = count_file)
  data.table::setnames(count_table, colnames(count_table), c("ENSG", "count"))

  #check the count file format
  if (ncol(count_table) != 2 | !is.numeric(count_table[, count])) {
    return(paste0("Incorrect count file format for the sample: ", sample_n))
  }

  #inner join reference and analyzed sample
  data.table::setkey(count_table, ENSG)
  data.table::setkey(diploid_standard, ENSG)
  final_mat <- as.data.frame(count_table[diploid_standard, nomatch = 0])

  #keep genes for determining gender for later
  gender_genes = final_mat %>% filter(ENSG %in% c("ENSG00000114374", "ENSG00000012817", "ENSG00000260197", "ENSG00000183878"))

  #filter genes based on reads count; top 1-q have read count > N, filter base on weight
  keepIdx = final_mat %>% mutate(keep_gene = apply(.[, -1], MARGIN = 1, FUN = function(x) sum(x > minReadCnt) > (length(x) * samp_prop)), id = row_number()) %>% filter(keep_gene == TRUE) %>%
    inner_join(weight_table, by = "ENSG") %>% group_by(chromosome_name) %>% mutate(weight_chr_quant = quantile(weight, 1 - weight_samp_prop)) %>% filter(weight > weight_chr_quant) %>% pull(id)

  # filter table for normalization, get rid of genes with 0 counts and keep genes for gender estimation
  count_filt <- final_mat %>% .[c(keepIdx), ] %>% .[pull(., 2) != 0, ] %>% bind_rows(gender_genes)
  ENSG <- count_filt$ENSG
  count_filt <- select(count_filt, -ENSG)

  #calculate per-gene standard geom_mean and remove zeros for size factor calculation
  pseudo_ref_des <- apply(X = count_filt, MARGIN = 1, function(x) gm_mean(x))
  low_exp <- which(pseudo_ref_des == 0)
  size_fac <- apply(count_filt[-low_exp, ], MARGIN = 2, function(x) median(x/pseudo_ref_des[-low_exp]))

  #Divide each column by size factor
  count_norm <- count_filt
  for (i in 1:ncol(count_filt)) {
    count_norm[, i] <- count_filt[, i]/size_fac[i]
  }

  print(paste0("Normalization for sample: ", sample_n, " completed"))

  #Modify table for downstream analysis
  rownames(count_norm) <- ENSG
  colnames(count_norm)[1] <- sample_n

  return(count_norm)
}

#calculate geometric mean
gm_mean = function(a){prod(a)^(1/length(a))}

####get median expression level####
get_med <- function(count_norm, refDataExp) {

  ENSG=rownames(count_norm)

  ####calculate median for all genes####
  pickGeneDFall=count_norm %>% mutate(ENSG=ENSG) %>% left_join(select(refDataExp, chr, ENSG), by = "ENSG") %>%
    mutate(med = apply(.[, -c(ncol(.) - 1, ncol(.))], 1, median)) %>%
    select(ENSG, med)

  return(pickGeneDFall)
}

####prepare snv files####
prepare_snv <- function(sample_table, centr_ref, sample_num, minDepth, chrs, snv_format = c("vcf", "custom")) {

  snv_file <- pull(sample_table, snv_path)[sample_num]
  sample_n <- pull(sample_table, 1)[sample_num]

  smpSNP=list()
  print(paste("Preparing file with snv information for:", sample_n))

  #prepare data from custom table
  if (snv_format == "custom") {
    #read the table
    snv_table_pre <- fread(snv_file)

    #check if all appropriate columns are present
    cols <- colnames(snv_table_pre)
    chr <- str_which(cols, "^#Chromosome$|#CHROM$|^CHR$|^chr$|^Chr$|^CHROM$|^chrom$|^Chrom$|^CHROMOSOME$|^chromosome$|^Chromosome$")
    start <- str_which(cols, "^START$|^start$|^Start$|^POS$|^pos$|^Pos$")
    depth <- str_which(cols, "^DEPTH$|^depth$|^Depth$|^DP$|^dp$|^Dp$")
    maf <- str_which(cols, "^MAF$|^maf$|^Maf$|^HET$|^het$|^Het$")
    to_keep <- as.numeric(c(chr, start, depth, maf))
    if (length(to_keep) != 4) {
      smpSNP[[sample_n]] <- "Incorrect column name in a custom snv file."
      return(smpSNP)
    }
    snv_table <- snv_table_pre[, to_keep, with = FALSE]
    data.table::setnames(snv_table, colnames(snv_table), c("chr", "start", "depth", "maf"))

    #Check some column parameters
    if (is.numeric(snv_table[, start]) == FALSE | is.numeric(snv_table[, depth]) == FALSE | is.numeric(snv_table[, maf]) == FALSE) {
      smpSNP[[sample_n]] <- "Incorrect type of a column in a custom file with snv information."
      return(smpSNP)
    }

  } else if (snv_format == "vcf") {
    snv_table <- vcf_to_snv(snv_file)
    if(is.character(snv_table)) {
      smpSNP[[sample_n]] <- snv_table
      return(smpSNP)
    }
  }

  smpSNP[[sample_n]] <- snv_table %>% filter(chr %in% chrs) %>%
    mutate(chr = factor(chr, levels=chrs), ID=paste0(chr,"-", start), sampleID=sample_n) %>% left_join(centr_ref, by = "chr") %>% mutate(arm = ifelse(start < cstart, "p", ifelse(start > cend, "q", "centr")))

  return(smpSNP)
}


####filter SNVs of interest for samples###
filter_snv <- function(one_smpSNP, keepSNP, minDepth) {
  smpSNPdata.tmp= one_smpSNP %>% dplyr::select(sampleID, ID, maf, chr, start, depth, arm) %>%  filter(
    data.table::inrange(maf, 0.05, 0.9),
    depth > minDepth,
    ID %in% keepSNP) %>% filter(chr != "Y")
  return(smpSNPdata.tmp)
}


####Calculate peak statistics for whole chromosome####
calc_chrom_lvl <- function(smpSNPdata.tmp) {
  smpSNPdata <- smpSNPdata.tmp %>% group_by(chr) %>% arrange(chr, desc(depth) ) %>%
    mutate(snvOrd=1:n()) %>% filter(snvOrd<=1000) %>%
    mutate(snvNum=n(), peak_max=densityMaxY(maf),
           peak=findPeak(maf), peakCol=ifelse(between(peak, 0.42, 0.58), 'black', 'red'), peakdist = find_peak_dist(maf)) %>%
   ungroup() %>% mutate(chr = factor(chr, levels = c(1:22, "X")))
  return(smpSNPdata)
}

####Calculate peak statistics for arms separately#s###
calc_arm_lvl <- function(smpSNPdata.tmp) {
  smpSNPdata_a=smpSNPdata.tmp %>% group_by(chr, arm) %>% arrange(chr, desc(depth) ) %>%
    mutate(snvNum=n(), peak_max=densityMaxY(maf),
           peak=findPeak(maf), peakCol=ifelse(between(peak, 0.42, 0.58), 'black', 'red'), peakdist = find_peak_dist(maf)) %>%
    ungroup() %>% mutate(chr = factor(chr, levels = c(1:22, "X")))
  return(smpSNPdata_a)
}

####normalize normalized counts (against median of expression for each gene) and join weight values####
#beta-needs cleaning
count_transform <- function(count_ns, pickGeneDFall, refDataExp, weight_table) {
  count_ns_tmp = count_ns %>% left_join(pickGeneDFall, by = "ENSG") %>%
    mutate(count_nor_med=log2(.[, 1] / med) ) %>% filter(med != 0)
  sENSGinfor=refDataExp[match(count_ns_tmp$ENSG, refDataExp$ENSG), ] %>% select(chr, end, start)

  #keeping only the genes which have weights calculated for geom_poit and boxplot
  count_ns = cbind(sENSGinfor, count_ns_tmp) %>% inner_join(weight_table, by = "ENSG")
}

####filter out par regions####
remove_par <- function(count_ns, par_reg) {
  #### get rid of PAR regions

  parX = filter(par_reg, chr == "X")
  count_ns = count_ns %>% filter(chr != "X" | ! data.table::inrange(start, parX[1, 1], parX[1, 2]) & ! data.table::inrange(start, parX[2, 1], parX[2, 2]) &
                                 ! data.table::inrange(end, parX[1, 1], parX[1, 2]) & !  data.table::inrange(end, parX[2, 1], parX[2, 2]))
  parY = filter(par_reg, chr == "Y")
  count_ns = count_ns %>% filter(chr != "Y" | ! data.table::inrange(start, parY[1, 1], parY[1, 2]) & ! data.table::inrange(start, parY[2, 1], parY[2, 2]) &
                             ! data.table::inrange(end, parY[1, 1], parY[1, 2]) & !  data.table::inrange(end, parY[2, 1], parY[2, 2]))
}

####Calculate weighted boxplot values####
get_box_wdt <- function(count_ns, chrs, scaleCols) {
  box_wdt <- count_ns %>% filter(chr %in% c(1:22, "X"))  %>% group_by(chr) %>% mutate(med_weig = spatstat::weighted.median(x = count_nor_med,w = weight, na.rm = TRUE), low = spatstat::weighted.quantile(x = count_nor_med, w = weight, probs = 0.25),
                                                                                   high = spatstat::weighted.quantile(x = count_nor_med, w = weight, probs = 0.75), IQR = abs(high - low), max = high + IQR*1.5, min = low - IQR*1.5)

  colours <- c()
  for(i in 1:nrow(box_wdt)) {
    colours[i] <- scaleCols$colour[which.min(abs(box_wdt$med_weig[i] - scaleCols$med_expr))]
  }

  box_wdt <- box_wdt %>% ungroup() %>% mutate(medianCol = colours) %>%
    distinct(chr, .keep_all = TRUE) %>% select(chr, med_weig, low, high, min, max, medianCol) %>%
    mutate(chr = factor(x = chr, levels = c(1:22, "X")), pos = 0.5)

  return(box_wdt)
}

#### change ylim accordingly values in graph
adjust_ylim <- function(box_wdt, ylim) {
  box_wdt_noy <- filter(box_wdt, !chr %in% "Y")
  if (any(box_wdt_noy$max > ylim[2])) {
    ylim[2] <- max(box_wdt_noy$max)*1.1
  }

  if (any(box_wdt_noy$min < ylim[1])) {
    ylim[1] <- min(box_wdt_noy$min)*1.1
  }
  return(ylim)
}

####Apply limit to datapoints and normalize point position####
prep_expr <- function(count_ns, dpRatioChrEdge, ylim, chrs) {
  count_ns_final= count_ns %>% select(chr, end, count_nor_med, weight) %>%
    bind_rows(dpRatioChrEdge) %>% filter(chr %in% c(1:22, "X"), between(count_nor_med, ylim[1], ylim[2]) ) %>%
    mutate(chr=factor(chr, levels = c(1:22, "X"))) %>%
    arrange(chr, end) %>% group_by(chr) %>%
    mutate(normPos=scales::rescale(end, from = range(end, na.rm = TRUE)))
  return(count_ns_final)
}

filter_expr <- function(count_ns_final, cutoff = 0.6) {
  count_ns_final <- count_ns_final %>% filter(weight > quantile(weight, cutoff, na.rm = TRUE))
}

####plot expression boxplot and point plot####
plot_exp <- function(count_ns_final, box_wdt, sample_name, ylim, estimate, feat_tab_alt, gender) {
  gp_expr <- ggplot() + ylim(ylim) + ylab("Normalized expression") +
    scale_fill_identity()+
    geom_point(data = count_ns_final, aes(x = normPos, y = count_nor_med, size = weight), alpha = 0.32, show.legend = FALSE)+
    scale_size(range = c(2, 6)) +
    #scale_alpha(range = c(0.22, 0.4)) +
    geom_boxplot(data = box_wdt, aes(ymin = min, lower = low, middle = med_weig, upper = high, ymax = max, fill=medianCol, x = pos), alpha=0.75, outlier.colour = NA, stat = "identity", show.legend = FALSE)+
    geom_hline(yintercept = 0, colour = "red")+
    labs(title = paste0(sample_name),
    subtitle = paste0("estimated gender: ", gender))
  if (estimate == TRUE) {
    gp_expr <- gp_expr +
      geom_point(data = data.frame(x = c(0.5, 0.5), y = c(ylim[2], ylim[2]), point_col = c("low", "high"), chr = factor(c(1, 1), levels = c(1:22, "X"))), mapping = aes(x = x, y = y, color = point_col), shape = 0, size = 4, stroke = 2) +
      geom_label(data = distinct(feat_tab_alt, chr, colour_chr, chr_alt), aes(x = 0.5, y = ylim[2], color = colour_chr, label = chr_alt), label.size = 2, show.legend = FALSE) +
      scale_color_manual(limits = c("low", "high"), values=c("orangered", "black")) +
      guides(color = guide_legend(
        title = "Quality"
      ))
    }

  gp_expr <- gp_expr +
    facet_grid(.~chr) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.title.y=element_text(size = 12),
          axis.ticks = element_blank(),
          legend.justification = "top",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17)) +
    guides(size = FALSE)
}

####plot snv density plots####
plot_snv <- function(smpSNPdata, chrs, sample_name, estimate) {

  missedChr=c(1:22, "X")[table(smpSNPdata$chr) < 15]
  if(length(missedChr) > 0){
    tmpSNPdata=data.frame(sampleID = sample_name, ID=paste0(missedChr, "-1"), maf=0.5, chr=factor(missedChr, levels = c(1:22, "X")), start=1,
                          depth=100, snvOrd=1, snvNum=1, peak_max=0, peak=0, peakCol="red", stringsAsFactors = F)
    smpSNPdata = bind_rows(smpSNPdata, tmpSNPdata)
  }
  if(nrow(smpSNPdata)<500){
    gp.maf=ggplot()+annotate("text", x = 1, y = 1, label = "Low coverage")+theme_void()
  }else{
    snvNumDensityMaxY=smpSNPdata %>% select(chr, snvNum, peak_max, peakdist) %>% unique()
    yAxisMax=snvNumDensityMaxY %>% filter(snvNum > 100) %>% .$peak_max %>% max()
    snvNumDF = snvNumDensityMaxY %>% mutate(x=0.5, y=yAxisMax*1.05)
    peakdist_dat = snvNumDensityMaxY %>% mutate(x = 0.5, y = yAxisMax*1.15, label = round(peakdist, 3))
    gp.maf=ggplot(data=smpSNPdata) + xlab("Mutant allele frequency") + ylab("Density") +
      geom_density(aes(maf, color=peakCol), show.legend = FALSE) +
      geom_text(data = peakdist_dat, aes(x, y, label = label), vjust=0)+
      geom_vline(xintercept = c(1/3, 0.5, 2/3), alpha = 0.4, size = 0.5)+
      scale_color_identity(guide = guide_legend(override.aes = list(color = "white")))+
      scale_x_continuous(breaks = round(c(1/3, 2/3), 3), labels = c("1/3", "2/3"), minor_breaks = NULL, limits = c(0,1)) +
      scale_y_continuous(breaks = c(seq(from = 1, to = floor(yAxisMax)), yAxisMax*1.15), labels = c(seq(from = 1, to = floor(yAxisMax)), "peak dist."), limits = c(0, yAxisMax*1.2)) +
      facet_grid(.~chr, scales="free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 40, vjust = 0.5),
            axis.ticks = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.title.y = element_text(size = 12),
            axis.text.y = element_text(size = 7),
            plot.margin = unit(c(0,1,1,1), "lines")
      )
    if (estimate == TRUE) {
      #need to add identical legend as for the expression in order for the graphs to align correctly with ggarrange
      gp.maf <- gp.maf +
        geom_point(data = data.frame(x = c(0.5, 0.5), y = c(0, 0), point_col = c("white", "white"), chr = factor(c(1, 1), levels = c(1:22, "X"))), mapping = aes(x = x, y = y, color = point_col), shape = 20, size = 5, alpha = 0) +
        guides(color = guide_legend(
          title = "Quality"
        )) +
        theme(
          legend.text = element_text(color = "white", size = 15),
          legend.title = element_text(color = "white", size = 17),
          legend.key = element_blank(),
          legend.justification = "top"
        )
    }
  }
  return(gp.maf)
}

####arrange expression and snv graph####
arrange_plots <- function(gg_exp, gg_snv) {
  fig <- ggarrange(plotlist =list(expr=gg_exp, maf=gg_snv)
                , ncol = 1, nrow = 2, heights = c(3, 1), align='v')
  return(fig)
}


####arm level utilities####
#### rescale centromeric region ####
rescale_centr <- function(centr_ref, count_ns_final) {
  count_ns_range <- count_ns_final %>% group_by(chr) %>% summarise(chr_end = max(end)) %>% mutate(chr_start = 1)
  centr_res <- centr_ref  %>% filter(chr != "Y") %>% mutate(chr = factor(chr, levels = c(1:22, "X"))) %>% left_join(count_ns_range, by = "chr") %>% mutate(p_mid = cstart/2, q_mid = cend + (chr_end - cend)/2)
  rescaled <- data.frame(cstartr = c(), cendr = c(), p_midr = c(), q_midr = c())

  for (i in 1:nrow(centr_res)) {
    cstartr_chr <- scales::rescale(centr_res$cstart[i], from = c(centr_res$chr_start[i], centr_res$chr_end[i]))
    cendr_chr <- scales::rescale(centr_res$cend[i], from = c(centr_res$chr_start[i], centr_res$chr_end[i]))
    p_midr <- scales::rescale(centr_res$p_mid[i], from = c(centr_res$chr_start[i], centr_res$chr_end[i]))
    q_midr <- scales::rescale(centr_res$q_mid[i], from = c(centr_res$chr_start[i], centr_res$chr_end[i]))
    rescaled <- rbind(rescaled, c(cstartr_chr, cendr_chr, p_midr, q_midr))
  }
  colnames(rescaled) <- c("cstartr", "cendr", "p_midr", "q_midr")
  centr_res <- cbind(centr_res, rescaled) %>% select(chr, cstartr, cendr, p_midr, q_midr)

  return(centr_res)
}

####draw zoomed in expression graph####
plot_exp_zoom <- function(count_ns_final, centr_res, plot_chr, estimate, feat_tab_alt) {

  #filter only for chromosome of interest

  count_ns_chr <- filter(count_ns_final, chr == plot_chr, count_nor_med < 2.1 & count_nor_med > -2.1)

  gg_expr_zoom = ggplot(data=count_ns_chr) + ylim(c(-2.1, 2.1)) + ylab("Normalized expression") +
    geom_point(aes(x = normPos, y = count_nor_med, size = weight), alpha=0.6) +
    scale_size(range = c(1,5)) +
    geom_smooth(aes(x = normPos, y = count_nor_med, weight = weight), alpha = 0.5, size = 0.5, method = "loess", formula = 'y ~ x') +
    annotate("segment", x = 0, xend = 1, y = 0, yend = 0,
             colour = "red", alpha = 0.85) +
    scale_x_continuous(expand = c(0,0)) +
    geom_vline(data = filter(centr_res, chr == plot_chr), mapping = aes(xintercept = cstartr)) +
    geom_vline(data = filter(centr_res, chr == plot_chr), mapping = aes(xintercept = cendr)) +
    ggtitle(paste0("chromosome ", plot_chr)) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 17),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank())

  if (estimate == TRUE) {

    q_alt <- feat_tab_alt %>% filter(arm == "q", chr == plot_chr)
    p_alt <- feat_tab_alt %>% filter(arm == "p", chr == plot_chr)

    gg_expr_zoom <- gg_expr_zoom +
      geom_label(data = q_alt, aes(x = centr_res$q_midr[centr_res$chr == plot_chr], y = 0.26, label = paste0(alteration, ", " , alteration_prob*100, "%"), color = colour_arm, size = 15000), nudge_y = 1) +
      scale_color_manual(limits = c("low", "high"), values=c("orangered", "black")) +
      if(nrow(p_alt) > 0) {
        geom_label(data = p_alt, aes(x = centr_res$p_midr[centr_res$chr == plot_chr], y = 0.26, label = paste0(alteration, ", " , alteration_prob*100, "%"), color = colour_arm, size = 15000), nudge_y = 1)
      }
  }



  return(gg_expr_zoom)
}


####calculate per chromosome max of y axis for arm level MAF graphs
get_yAxisMax_arm <- function(smpSNPdata_arm, plot_chr) {

  smpSNP_arm_chr <- smpSNPdata_arm %>% filter(chr == plot_chr)
  yAxisMax_arm <- smpSNP_arm_chr %>% distinct(arm, peak_max) %>% summarise(y_max = max(peak_max)) %>% pull(y_max)
  return(yAxisMax_arm)

}

####draw snv graph for chromosome arm####
plot_snv_arm <- function(smpSNPdata_a, plot_arm, plot_chr, yAxisMax) {

  smpSNP_arm <- filter(smpSNPdata_a, chr == plot_chr, arm == plot_arm)

  SNP_ann <- smpSNP_arm %>% select(chr, snvNum, peak_max, peakdist) %>% unique()

  peakdist_dat = SNP_ann %>% mutate(x = 0.5, y = yAxisMax*1.15, label = round(peakdist, 3))

  if (nrow(smpSNP_arm) < 10) {
    gg_snv_arm = ggplot() + ylim(0, yAxisMax*1.2)+
      scale_x_continuous(limits = c(0,1)) +
      annotate("text", x = 0.5, y = 0.8 * yAxisMax, label = "Low coverage")+
      ggtitle(paste0(plot_arm, " arm")) +
      theme_void() +
      theme(plot.margin = unit(c(0,1,1,1), "lines"), plot.title = element_text(hjust = 0.5, size = 20))
  } else {
    gg_snv_arm <- ggplot(data=smpSNP_arm) + xlab("Mutant allele frequency") + ylab("Density") + ylim(0, yAxisMax*1.2) +
      geom_density(aes(maf, color=peakCol)) +
      geom_text(data = peakdist_dat, aes(x, y, label = paste0("peak dist. = ", label)), vjust=0, size = 6)+
      geom_vline(xintercept = c(1/3, 0.5, 2/3), alpha = 0.4, size = 0.5)+
      scale_color_identity()+
      scale_x_continuous(breaks = round(c(1/3, 0.5, 2/3), 2), minor_breaks = NULL, limits = c(0,1), labels = c("1/3", "1/2", "2/3")) +
      ggtitle(paste0(plot_arm, " arm")) +
      theme_bw() +
      theme(axis.text = element_text(size = 16), axis.ticks = element_blank(),
            strip.background = element_blank(), strip.text.x = element_blank(),
            plot.margin = unit(c(0,1,1,1), "lines"), plot.title = element_text(hjust = 0.5, size = 20),
            axis.title = element_text(size = 17)
      )
  }

  # switcht axis side and remove q arm y axis title
  if (plot_arm == "q") {
    suppressMessages(gg_snv_arm <- gg_snv_arm +
                       theme(axis.title.y.left = element_blank(),
                             axis.text.y.left = element_blank(),
                             axis.ticks.y.left = element_blank()) +
                       scale_y_continuous(position = "right") +
                       ylab("Density"))
  }

  return(gg_snv_arm)
}

####arrange arm graphs####
chr_plot <- function(p_snv, q_snv, arm_expr) {
  chr_arm_level <- ggarrange(plotlist = list(p_snv, arm_expr, q_snv), ncol = 3, widths = c(1.2, 4, 1.2), align = "h")
  return(chr_arm_level)
}

####modify decimals####
specify_decimal <- function(x, k) {
  trimws(format(round(x, k), nsmall=k))
}

####find peak in a vector####
findPeak <- function(vec){
  vec=vec[between(vec, 0.1, 0.9)]
  len=length(vec)
  if(len < 10){
    return(0)
  }else{
    d=density(vec)
    maxVec=max(d$y)
    maxPos=d$x[d$y == maxVec]
    # peak=d$x[which(d$y==max(d$y[which(diff(sign(diff(d$y) ))==-2)]))]
    return(maxPos[1])
  }
}

####find peak distance####
find_peak_dist <- function(vec) {
  len=length(vec)
  if(len < 10){
    return(0)
  }else{
    d=density(vec)
    #
    peaks_max = d$x[which(d$y %in% sort(d$y[which(diff(sign(diff(d$y) ))==-2)], decreasing = TRUE)[1:2])]
    #making sure the distance is measured between symetric peaks
    #symmetry treshold set at 1.2
    if (length(na.omit(peaks_max)) > 1 & data.table::inrange(max(peaks_max - 0.5), 0.42 - min(peaks_max), 0.58 - min(peaks_max))) {
      dist_peak <- abs(peaks_max[1] - peaks_max[2])
      return(dist_peak)
    } else {
      return(0)
    }
  }
}

####find max on y axis from density vector of allele frequency####
densityMaxY <- function(vec){
  len=length(vec)
  if(len < 10){
    return(0)
  }else{
    d=density(vec)
    return(max(d$y))
  }
}

####adjust for diploid level based on diploid chromosomes####
adjust_dipl <- function(feat_tab_alt, count_ns) {

  if (sum(feat_tab_alt$chr_status == "dipl") < 1) {
    print("Unable to adjust expression level")
    return(count_ns)
  }

  baseline_shift = median(feat_tab_alt$arm_med[feat_tab_alt$chr_status == "dipl"])
  count_ns$count_nor_med <- count_ns$count_nor_med - baseline_shift
  return(count_ns)
}

####get arm metrics####
get_arm_metr <- function(count_ns, smpSNPdata, sample_name, centr_ref, chrs) {

  #calculate weighted median for every chromosome and use only 1:22
  summ_arm <- count_ns %>% filter(!is.infinite(count_nor_med)) %>% mutate(chr = factor(chr, levels = chrs)) %>% filter(chr %in% c(1:22, "X")) %>% left_join(centr_ref, by = "chr") %>%
    mutate(arm = ifelse(end < cstart, "p", ifelse(end > cend, "q", "centr"))) %>% group_by(chr, arm) %>%

    # get rid of 21 p arm
    filter(!(chr == 21 & arm == "p")) %>%

    mutate(arm_med = spatstat::weighted.median(x = count_nor_med,w = weight, na.rm = TRUE),
           up_quart = spatstat::weighted.quantile(x = count_nor_med, w = weight, probs = 0.75, na.rm = TRUE),
           low_quart = spatstat::weighted.quantile(x = count_nor_med, w = weight, probs = 0.25, na.rm = TRUE)) %>%

    ungroup() %>% distinct(chr, arm, arm_med, up_quart, low_quart) %>%
    left_join(distinct(.data = smpSNPdata, chr, arm, peakdist, peak_m_dist, peak_max), by = c("chr", "arm")) %>%

    mutate(chr = factor(chr, levels = c(1:22, "X")), sd = sd(arm_med), mean_all = mean(arm_med)) %>% ungroup() %>%

    mutate(sds_median = (arm_med - mean_all)/sd, sds_025 = (low_quart - mean_all)/sd, sds_075 = (up_quart-mean_all)/sd,
           n_02_04 = sum(data.table::inrange(peakdist, 0.2, 0.4) & peak_m_dist > 0.08), n_04 = sum(data.table::inrange(peakdist, 0.4, 0.9) & peak_m_dist > 0.08)) %>%

    arrange(chr) %>%

    select(chr, arm, arm_med, up_quart, low_quart, peak_max, peak_m_dist, peakdist, sd, sds_median, sds_025, sds_075, n_02_04, n_04)
  return(summ_arm)

}

#### calculate chromosomal statistics ####
calc_arm <- function(smpSNPdata.tmp) {
  smpSNPdata <- smpSNPdata.tmp %>% group_by(chr, arm) %>% arrange(chr, desc(depth) ) %>%
    mutate(snvOrd=1:n()) %>%
    mutate(snvNum=n(), peak_max=densityMaxY(maf),
           peak=findPeak(maf), peak_m_dist = abs(peak - 0.5), peakdist = find_peak_dist(maf), peakCol=ifelse(between(peak, 0.42, 0.58), 'black', 'red')) %>%
    ungroup() %>% mutate(chr = factor(chr, levels = c(1:22, "X", "Y")))
  return(smpSNPdata)
}

#### calculate statistics with diploid knowldedge ####
metr_dipl <- function(data) {

  if (sum(data$chr_status == "dipl") > 3) {
    sd_chr = "dipl"
  } else {
    sd_chr = "no_dipl"
  }

  sd_dipl <- data %>% filter(chr_status == sd_chr) %>% mutate(sd_dipl = sd(arm_med)) %>% pull(sd_dipl) %>% unique()
  mean_dipl <- data %>% filter(chr_status == "dipl") %>% mutate(mean_dipl = mean(arm_med)) %>% pull(mean_dipl) %>% unique()

  data_mod <- data %>% mutate(sd_dipl = sd_dipl, mean_dipl = mean_dipl, sds_median_dipl = (arm_med - mean_dipl)/sd_dipl, sds_025_dipl = (low_quart - mean_dipl)/sd_dipl, sds_075_dipl = (up_quart-mean_dipl)/sd_dipl)

  return(data_mod)
}

### Create colour coding for estimation values ####
colour_code <- function(data, conf_tresh) {
  data_col <- data %>% mutate(colour_arm = factor(ifelse(alteration_prob < conf_tresh, "low", "high"), levels = c("low", "high"))) %>% group_by(chr) %>% mutate(min_prob = min(alteration_prob)) %>% ungroup() %>% mutate(colour_chr = factor(ifelse(min_prob < conf_tresh, "low", "high"), levels = c("low", "high"), ordered = TRUE)) %>%
    select(-min_prob)
  return(data_col)
}

### Generate report for sample ###
gen_kar_list <- function(feat_tab_alt, sample_name, gender) {

  ##extract only the chromosomes which have only one call and have high quality
  tab_mod <- feat_tab_alt %>% select(chr, arm, alteration, colour_arm) %>% group_by(chr) %>% mutate(alt_types = length(unique(alteration)), qual_types = length(unique(colour_arm))) %>% ungroup()

  # whole chromosome changes
  alt_chr_whole <- tab_mod %>% filter(alt_types == 1 & qual_types == 1 & (colour_arm != "high" | alteration != 0)) %>%
    mutate(alt_sign = ifelse(alteration %in% c(1, 2), "+", ifelse(alteration == -1, "-", "")), qual_sign = ifelse(colour_arm == "high", "", "?"), alt_str = paste0(qual_sign, chr, alt_sign)) %>%
    distinct(chr, alt_str, alteration)
  double_chr_whole <- alt_chr_whole %>% filter(alteration == 2)
  alt_chr_whole_fin <- alt_chr_whole %>% bind_rows(double_chr_whole) %>% arrange(chr)

  # arm only changes
  alt_arm <- tab_mod %>% filter((alt_types == 2 | qual_types == 2) & (colour_arm != "high" | alteration != 0)) %>%
    mutate(alt_sign = ifelse(alteration %in% c(1, 2), "+", ifelse(alteration == -1, "-", "")), qual_sign = ifelse(colour_arm == "high", "", "?"), alt_str = paste0(qual_sign, chr, arm, alt_sign)) %>%
    distinct(chr, alt_str, alteration, arm)
  double_arm <- alt_arm %>% filter(alteration == 2)
  alt_arm_fin <- alt_arm %>% bind_rows(double_arm) %>% arrange(chr, arm) %>% select(-arm)

  #fuse the tables
  alterations <- bind_rows(alt_chr_whole_fin, alt_arm_fin) %>% arrange(chr) %>% pull(alt_str) %>% paste0(collapse = ", ")

  #chromosome number
  chrom_diff <- tab_mod %>% filter(alt_types == 1, qual_types == 1, colour_arm == "high") %>% distinct(chr, alteration) %>% pull(alteration) %>% as.numeric(.) %>% sum(.)
  chrom_n = 46 + chrom_diff

  #fill in empty alteration vector
  if (alterations == "") {
    alterations <- "none"
  }


  kar_table <- data.frame(sample = sample_name,
                          gender = factor(gender, levels = c("female", "male")),
                          chrom_n = chrom_n,
                          alterations = alterations, stringsAsFactors = FALSE)

  return(kar_table)
}

#' Fucnction that converts vcf files to tabular input for CNV analysis
vcf_to_snv <- function(vcf_file, maf_tresh = 0.01, depth_tresh = 5) {

  #read the vcf files
  message("Reading in vcf file..")
  vcf_data <- fread(vcf_file)
  #check whether the input is in correct format
  if (dim(vcf_data)[2] < 10) {
    vcf_final <- "Incorrect vcf file format. Incorrect number of columns"
    return(vcf_final)
  }
  if (!identical(colnames(vcf_data)[2:8], c("POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))) {
    vcf_final <- "Incorrect vcf file format."
    return(vcf_final)
  }
  if (str_detect(vcf_data[1, INFO], "DP=[0-9]+") != TRUE) {
    vcf_final <- "Incorrect vcf file format. No depth information in INFO column"
    return(vcf_final)
  }
  if (str_detect(vcf_data[1, 9], "AD") != TRUE) {
    vcf_final <- "Incorrect vcf file format. No allele depth (AD) in FORMAT column"
    return(vcf_final)
  }


  vcf_data <- vcf_data[, 1:10]
  data.table::setnames(x = vcf_data, old = colnames(vcf_data), new = c("chr", "start", "ID", "ref", "var", "qual", "FILTER", "INFO", "FORMAT", "INFO2"))

  # Getting depth out of the INFO column
  message("Extracting depth..")
  vcf_data[, depth := as.numeric(sub("DP=", "", str_extract(INFO, "DP=[0-9]+")))]

  #reference allele and alternative allele depths
  message("Extracting reference allele and alternative allele depths..")
  vcf_data[, REF_ALT_num := sapply(str_split(INFO2, ":"), function(x) x[2])]
  #extract count for alternative allele
  vcf_data[, varDp := as.numeric(sapply(str_split(REF_ALT_num, ","), function(x) x[2]))]
  #mutant allele frequency
  vcf_data[, maf := varDp/depth]

  message("Needed information from vcf extracted")

  #return needed columns
  vcf_final <- vcf_data[, c("chr", "start", 'depth', "maf")]

  message("Finished reading vcf")
  return(vcf_final)
}
