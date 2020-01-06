#' Importing packages
#' @import shiny
#' @import tidyr
#' @import dplyr
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import stringr
#' @importFrom data.table fread
#' @import ggpubr
#' @import DESeq2
#' @importFrom spatstat weighted.median
#' @importFrom spatstat weighted.quantile

#### get vst values ####
get_vst <- function(sample_table, minReadCnt, q, sample_num, base_col, base_matr, weight_table, keep_perc) {

  count_file <- pull(sample_table, count_path)[sample_num]
  sample_n <- pull(sample_table, 1)[sample_num]

  #read the tables in
  count_table <- read.table(file = count_file, header = FALSE, row.names = 1, stringsAsFactors = FALSE)

  #check the count file format
  if(ncol(count_table) != 1 | typeof(count_table[, 1]) != "integer") return(NULL)

  colnames(count_table) <- sample_n

  #create coldata
  colData = data.frame(sampleName = sample_n, fileName = count_file)

  #bind with baseline
  final_mat <- cbind(count_table, base_matr)
  final_col <- rbind(colData, base_col)
  colnames(final_mat) <- as.character(final_col$sampleName)

  #calculate vst
  ddsHTSeq=DESeqDataSetFromMatrix(countData = final_mat, colData = final_col, design= ~ 1)
  print("Loading HTSeq files is done!")

  ddsCount <- counts(ddsHTSeq)

  #filter genes based on reads count; top 1-q have read count > N and filter based on weight
  keepIdx = as.data.frame(ddsCount) %>% mutate(id = row_number(), keep_count = apply(., MARGIN = 1, FUN = function(x) quantile(x, q)), ENSG = rownames(.)) %>% filter(keep_count >= minReadCnt) %>%
    left_join(weight_table) %>% filter(weight > quantile(.$weight, 1-keep_perc, na.rm = TRUE)) %>% pull(id)

  ddsHTSeq <- ddsHTSeq[keepIdx, ]
  rldHTSeq <- varianceStabilizingTransformation(ddsHTSeq, blind=T, fitType='local')

  print("VST transformation is done!")
  vst = round(data.frame(SummarizedExperiment::assay(rldHTSeq)), digits = 2)

  return(vst)
}



####get median expression level####
get_med <- function(vst, refDataExp) {

  ENSG=rownames(vst)

  ####calculate median for all genes####
  pickGeneDFall=vst %>% mutate(ENSG=ENSG) %>% left_join(select(refDataExp, chr, ENSG), by = "ENSG") %>%
    mutate(med= ifelse(chr != "Y", apply(.[, -c(ncol(.) - 1, ncol(.))], 1, median), apply(.[, -c(ncol(.) - 1, ncol(.))], 1, function(x) quantile(x = x, 0.60)))) %>%
    select(ENSG, med)

  return(pickGeneDFall)
}

####prepare snv files####
prepare_snv <- function(sample_table, centr_ref, sample_num, minDepth, chrs) {

  header <- c("chr", "start", "ref", 'var', 'end', 'qual', 'depth', 'refDp', 'varDp', 'mapQ', "maf")

  snv_file <- pull(sample_table, snv_path)[sample_num]
  sample_n <- pull(sample_table, 1)[sample_num]

  refSNP=c()
  smpSNP=list()
  print(paste("preparing MAF file:", sample_n) )
  rnaData <- fread(snv_file)
  if (ncol(rnaData) != 11 | all(sapply(rnaData, class) == c("character", "integer", "character", "character", "integer", "numeric", "integer", "integer", "integer", "numeric", "numeric")) != TRUE) return(NULL)
  rnaData <- rnaData %>% select(1:11) %>% magrittr::set_colnames(header)
  smpSNP[[sample_n]] <- rnaData %>% filter(chr %in% chrs, depth > minDepth ) %>%
    mutate(chr = factor(chr, levels=chrs), ID=paste0(chr,"-", start), sampleID=sample_n) %>% left_join(centr_ref, by = "chr") %>% mutate(arm = ifelse(start < cstart, "p", ifelse(start > cend, "q", "centr")))
  return(smpSNP)
}


####filter SNVs of interest for samples###
filter_snv <- function(one_smpSNP, keepSNP) {
  smpSNPdata.tmp= one_smpSNP %>% dplyr::select(sampleID, ID, maf, chr, start, depth, arm) %>%  filter(
    data.table::inrange(maf, 0.05, 0.9),
    ID %in% keepSNP) %>% filter(chr != "Y")
  return(smpSNPdata.tmp)
}


####Calculate peak statistics for whole chromosome####
calc_chrom_lvl <- function(smpSNPdata.tmp) {
  smpSNPdata <- smpSNPdata.tmp %>% group_by(chr) %>% arrange(chr, desc(depth) ) %>%
    mutate(snvOrd=1:n()) %>% filter(snvOrd<=1000) %>%
    mutate(snvNum=n(), densityMaxY=densityMaxY(maf),
           peak=findPeak(maf), peakCol=ifelse(between(peak, 0.42, 0.58), 'black', 'red'), peakdist = find_peak_dist(maf)) %>%
   ungroup() %>% mutate(chr = factor(chr, levels = c(1:22, "X")))
  return(smpSNPdata)
}

####Calculate peak statistics for arms separately#s###
calc_arm_lvl <- function(smpSNPdata.tmp) {
  smpSNPdata_a=smpSNPdata.tmp %>% group_by(chr, arm) %>% arrange(chr, desc(depth) ) %>%
    mutate(snvNum=n(), densityMaxY=densityMaxY(maf),
           peak=findPeak(maf), peakCol=ifelse(between(peak, 0.42, 0.58), 'black', 'red'), peakdist = find_peak_dist(maf)) %>%
    ungroup() %>% mutate(chr = factor(chr, levels = c(1:22, "X")))
  return(smpSNPdata_a)
}

####calculate normalized vst values and join weight values####
#beta-needs cleaning
vst_norm <- function(s_vst, pickGeneDFall, refDataExp, weight_tab_q) {
  s_vst_tmp= s_vst %>% left_join(pickGeneDFall, by = "ENSG") %>%
    mutate(vst_nor_med=log2(.[, 1] / med) )
  sENSGinfor=refDataExp[match(s_vst_tmp$ENSG, refDataExp$ENSG), ] %>% select(chr, end, start)

  #keeping only the genes which have weights calculated for geom_poit and boxplot
  s_vst = cbind(sENSGinfor, s_vst_tmp) %>% left_join(weight_tab_q, by = "ENSG")
  return(s_vst)
}

####filter out par regions####
remove_par <- function(s_vst, par_reg) {
  #### get rid of PAR regions

  parX = filter(par_reg, chr == "X")
  s_vst = s_vst %>% filter(chr != "X" | ! data.table::inrange(start, parX[1, 1], parX[1, 2]) & ! data.table::inrange(start, parX[2, 1], parX[2, 2]) &
                                 ! data.table::inrange(end, parX[1, 1], parX[1, 2]) & !  data.table::inrange(end, parX[2, 1], parX[2, 2]))
  parY = filter(par_reg, chr == "Y")
  s_vst = s_vst %>% filter(chr != "Y" | ! data.table::inrange(start, parY[1, 1], parY[1, 2]) & ! data.table::inrange(start, parY[2, 1], parY[2, 2]) &
                             ! data.table::inrange(end, parY[1, 1], parY[1, 2]) & !  data.table::inrange(end, parY[2, 1], parY[2, 2]))
}

####Calculate weighted boxplot values####
get_box_wdt <- function(s_vst, chrs, scaleCols) {
box_wdt <- s_vst %>% filter(chr %in% c(1:22, "X"))  %>% group_by(chr) %>% mutate(med_weig = weighted.median(x = vst_nor_med,w = weight, na.rm = TRUE), low = weighted.quantile(x = vst_nor_med, w = weight, probs = 0.25),
                                                                           high = weighted.quantile(x = vst_nor_med, w = weight, probs = 0.75), IQR = abs(high - low), max = high + IQR*1.5, min = low - IQR*1.5,
                                                                           wdt_median_vst= specify_decimal(med_weig, 2) , medianCol=scaleCols[wdt_median_vst]) %>%

    distinct(chr, .keep_all = TRUE) %>% select(chr, med_weig, low, high, min, max, medianCol) %>% ungroup() %>%
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

####Add edge datapoint and normalize point position and filter for low weight genes####
prep_expr <- function(s_vst, dpRatioChrEdge, ylim, chrs) {
  s_vst_final= s_vst %>% select(chr, end, vst_nor_med, weight) %>%
    bind_rows(dpRatioChrEdge) %>% filter(chr %in% c(1:22, "X"), between(vst_nor_med, ylim[1], ylim[2]) ) %>%
    mutate(chr=factor(chr, levels=chrs)) %>%
    arrange(chr, end) %>% group_by(chr) %>%
    mutate(normPos=scales::rescale(end, from = range(end, na.rm = TRUE)))
  return(s_vst_final)
}

filter_expr <- function(s_vst_final, cutoff = 0.6) {
  s_vst_final <- s_vst_final %>% filter(weight > quantile(weight, cutoff, na.rm = TRUE))
}

####plot expression boxplot and point plot####
plot_exp <- function(s_vst_final, box_wdt, sample_name, ylim, estimate, feat_tab_alt) {
  gp_expr <- ggplot() + ylim(ylim) + ylab("Normalized vst") + ggtitle(paste0(sample_name)) +
    scale_fill_identity()+
    geom_point(data = s_vst_final, aes(x = normPos, y = vst_nor_med, size = weight), alpha = 0.32)+
    scale_size(range = c(2, 6)) +
    #scale_alpha(range = c(0.22, 0.4)) +
    geom_boxplot(data = box_wdt, aes(ymin = min, lower = low, middle = med_weig, upper = high, ymax = max, fill=medianCol, x = pos), alpha=0.75, outlier.colour = NA, stat = "identity")+
    geom_hline(yintercept = 0, colour = "red")

  if (estimate == TRUE) {
    gp_expr <- gp_expr +
    geom_label(data = distinct(feat_tab_alt, chr, colour_chr, chr_alt), aes(x = 0.5, y = ylim[2], color = colour_chr, label = chr_alt)) +
    scale_color_manual(limits = c("bad", "good", "very_good"), values=c("red", "blueviolet", "blue"))
  }

  gp_expr <- gp_expr +
    facet_grid(.~chr) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank())
}

####plot snv density plots####
plot_snv <- function(smpSNPdata, chrs, sample_name) {

  missedChr=c(1:22, "X")[!c(1:22, "X") %in% smpSNPdata$chr]
  if(length(missedChr) > 0){
    tmpSNPdata=data.frame(sampleID = sample_name, ID=paste0(missedChr, "-1"), maf=0.5, chr=factor(missedChr, levels = chrs), start=1,
                          depth=100, snvOrd=1, snvNum=1, densityMaxY=0, peak=0, peakCol="red", stringsAsFactors = F)
    smpSNPdata = bind_rows(smpSNPdata, tmpSNPdata)
  }
  if(nrow(smpSNPdata)<500){
    gp.maf=ggplot()+annotate("text", x = 1, y = 1, label = "No MAF plot. Less than 1k SNVs!")+theme_void()
  }else{
    snvNumDensityMaxY=smpSNPdata %>% select(chr, snvNum, densityMaxY, peakdist) %>% unique()
    yAxisMax=snvNumDensityMaxY %>% filter(snvNum > 100) %>% .$densityMaxY %>% max()
    snvNumDF = snvNumDensityMaxY %>% mutate(x=0.5, y=yAxisMax*1.05, label=paste0("n=", snvNum))
    peakdist_dat = snvNumDensityMaxY %>% mutate(x = 0.5, y = yAxisMax*1.15, label = round(peakdist, 3))
    gp.maf=ggplot(data=smpSNPdata) + xlab("Mutant allele frequency") + ylab("Density") + ylim(0, yAxisMax*1.2) +
      geom_density(aes(maf, color=peakCol)) +
      geom_text(data=snvNumDF, aes(x, y, label=label), vjust=0)+
      geom_text(data = peakdist_dat, aes(x, y, label = label), vjust=0)+
      geom_vline(xintercept = c(1/3, 0.5, 2/3), alpha = 0.4, size = 0.5)+
      scale_color_identity()+
      scale_x_continuous(breaks = round(c(1/3, 2/3), 3), labels = c("1/3", "2/3"), minor_breaks = NULL, limits = c(0,1)) +
      facet_grid(.~chr, scales="free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 40, vjust = 0.5), axis.ticks = element_blank(),
            strip.background = element_blank(), strip.text.x = element_blank(),
            plot.margin = unit(c(0,1,1,1), "lines"))

  }
}

####arrange expression and snv graph####
arrange_plots <- function(gg_exp, gg_snv) {
  fig <- ggarrange(plotlist =list(expr=gg_exp, maf=gg_snv)
                , ncol = 1, nrow = 2, heights = c(3, 1), align='v')
  return(fig)
}


####arm level utilities####
#### rescale centromeric region ####
rescale_centr <- function(centr_ref, s_vst_final) {
  s_vst_range <- s_vst_final %>% group_by(chr) %>% summarise(chr_end = max(end)) %>% mutate(chr_start = 1)
  centr_res <- centr_ref %>% filter(chr != "Y") %>% left_join(s_vst_range, by = "chr") %>% mutate(p_mid = cstart/2, q_mid = cend + (chr_end - cend)/2)
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
plot_exp_zoom <- function(s_vst_final, centr_res, plot_chr, estimate, feat_tab_alt) {
  #filter only for chromosome of interest
  s_vst_chr <- filter(s_vst_final, chr == plot_chr)
  gg_expr_zoom = ggplot(data=s_vst_chr) + ylim(c(-0.3, 0.3)) + ylab("Normalized vst") +
    geom_point(aes(x = normPos, y = vst_nor_med, size = weight), alpha=0.6) +
    scale_size(range = c(1,5)) +
    geom_smooth(aes(x = normPos, y = vst_nor_med, weight = weight), alpha = 0.5, size = 0.5) +
    annotate("segment", x = 0, xend = 1, y = 0, yend = 0,
             colour = "red", alpha = 0.85) +
    scale_x_continuous(expand = c(0,0)) +
    geom_vline(data = filter(centr_res, chr == plot_chr), mapping = aes(xintercept = cstartr)) +
    geom_vline(data = filter(centr_res, chr == plot_chr), mapping = aes(xintercept = cendr)) +
    ggtitle(paste0("chromosome ", plot_chr)) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank())

  if (estimate == TRUE) {

    q_alt <- feat_tab_alt %>% filter(arm == "q", chr == plot_chr)
    p_alt <- feat_tab_alt %>% filter(arm == "p", chr == plot_chr)

    gg_expr_zoom <- gg_expr_zoom +
      geom_label(data = q_alt, aes(x = centr_res$q_midr[centr_res$chr == plot_chr], y = 0.26, label = paste0(alteration, ", " , alteration_prob*100, "%"), color = colour_arm)) +
      scale_color_manual(limits = c("bad", "good", "very_good"), values=c("red", "blueviolet", "blue")) +
      if(nrow(p_alt) > 0) {
        geom_label(data = p_alt, aes(x = centr_res$p_midr[centr_res$chr == plot_chr], y = 0.26, label = paste0(alteration, ", " , alteration_prob*100, "%"), color = colour_arm))
      }
  }



  return(gg_expr_zoom)
}


####calculate per chromosome max of y axis for arm level MAF graphs
get_yAxisMax <- function(smpSNPdata, plot_chr) {

  smpSNP_chr <- smpSNPdata %>% filter(chr == plot_chr)
  yAxisMax <- smpSNP_chr %>% select(snvNum, densityMaxY, peakdist) %>% unique() %>% .$densityMaxY %>% max()
  return(yAxisMax)

}

####draw snv graph for chromosome arm####
plot_snv_arm <- function(smpSNPdata_a, plot_arm, plot_chr, yAxisMax) {

  smpSNP_arm <- filter(smpSNPdata_a, chr == plot_chr, arm == plot_arm)

  SNP_ann <- smpSNP_arm %>% select(chr, snvNum, densityMaxY, peakdist) %>% unique()

  snvNumDF = SNP_ann %>% mutate(x=0.5, y=yAxisMax*1.05, label=paste0("n=", snvNum))
  peakdist_dat = SNP_ann %>% mutate(x = 0.5, y = yAxisMax*1.15, label = round(peakdist, 3))


  if (nrow(smpSNP_arm) < 10) {
    gg_snv_arm = ggplot() + ylim(0, yAxisMax*1.2)+
      scale_x_continuous(limits = c(0,1)) +
      annotate("text", x = 0.5, y = 0.8 * yAxisMax, label = "Low n. of SNVs")+
      ggtitle(paste0(plot_arm, " arm")) +
      theme_void() +
      theme(plot.margin = unit(c(0,1,1,1), "lines"), plot.title = element_text(hjust = 0.5, size = 20))
  } else {
    gg_snv_arm <- ggplot(data=smpSNP_arm) + xlab("Mutant allele frequency") + ylab("Density") + ylim(0, yAxisMax*1.2) +
      geom_density(aes(maf, color=peakCol)) +
      geom_text(data=snvNumDF, aes(x, y, label=label), vjust=0)+
      geom_text(data = peakdist_dat, aes(x, y, label = paste0("peak_dist = ", label)), vjust=0)+
      geom_vline(xintercept = c(1/3, 0.5, 2/3), alpha = 0.4, size = 0.5)+
      scale_color_identity()+
      scale_x_continuous(breaks = round(c(1/3, 0.5, 2/3), 3), minor_breaks = NULL, limits = c(0,1)) +
      ggtitle(paste0(plot_arm, " arm")) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 0.8), axis.ticks = element_blank(),
            strip.background = element_blank(), strip.text.x = element_blank(),
            plot.margin = unit(c(0,1,1,1), "lines"), plot.title = element_text(hjust = 0.5, size = 20))
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
adjust_dipl <- function(feat_tab_alt, s_vst) {

  if (sum(feat_tab_alt$chr_status == "dipl") < 1) {
    print("Unable to adjust expression level")
    return(s_vst)
  }

  baseline_shift = median(feat_tab_alt$medianvst_nor_med[feat_tab_alt$chr_status == "dipl"])
  s_vst$vst_nor_med <- s_vst$vst_nor_med - baseline_shift
  return(s_vst)
}

####get arm metrics####
get_arm_metr <- function(s_vst, smpSNPdata, sample_name, centr_ref) {

  #calculate weighted median for every chromosome and use only 1:22
  summ_arm <- s_vst %>% filter(chr %in% c(1:22, "X")) %>% left_join(centr_ref, by = "chr") %>%
    mutate(arm = ifelse(end < cstart, "p", ifelse(end > cend, "q", "centr"))) %>% group_by(chr, arm) %>%

    mutate(medianvst_nor_med = weighted.median(x = vst_nor_med,w = weight, na.rm = TRUE),
           up_quart = weighted.quantile(x = vst_nor_med, w = weight, probs = 0.75, na.rm = TRUE),
           low_quart = weighted.quantile(x = vst_nor_med, w = weight, probs = 0.25, na.rm = TRUE)) %>%

    ungroup() %>% distinct(chr, arm, medianvst_nor_med, up_quart, low_quart) %>%
    left_join(distinct(.data = smpSNPdata, chr, arm, peakdist, peak_m_dist, peak_max), by = c("chr", "arm")) %>%

    mutate(chr = factor(chr, levels = c(1:22, "X", "Y")), sd = sd(medianvst_nor_med), mean_all = mean(medianvst_nor_med)) %>% ungroup() %>%

    mutate(sds_median = (medianvst_nor_med - mean_all)/sd, sds_025 = (low_quart - mean_all)/sd, sds_075 = (up_quart-mean_all)/sd,
           n_02_04 = sum(data.table::inrange(peakdist, 0.2, 0.4) & peak_m_dist > 0.08), n_04 = sum(data.table::inrange(peakdist, 0.4, 0.9) & peak_m_dist > 0.08)) %>%

    arrange(chr) %>%
    select(chr, arm, medianvst_nor_med, up_quart, low_quart, peak_max, peak_m_dist, peakdist, sd, sds_median, sds_025, sds_075, n_02_04, n_04)
  return(summ_arm)

}

#### calculate chromosomal statistics ####
calc_arm <- function(smpSNPdata.tmp) {
  smpSNPdata <- smpSNPdata.tmp %>% group_by(chr, arm) %>% arrange(chr, desc(depth) ) %>%
    mutate(snvOrd=1:n()) %>%
    mutate(snvNum=n(), peak_max=densityMaxY(maf),
           peak=findPeak(maf), peak_m_dist = abs(peak - 0.5), peakdist = find_peak_dist(maf)) %>%
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

  sd_dipl <- data %>% filter(chr_status == sd_chr) %>% mutate(sd_dipl = sd(medianvst_nor_med)) %>% pull(sd_dipl) %>% unique()
  mean_dipl <- data %>% filter(chr_status == "dipl") %>% mutate(mean_dipl = mean(medianvst_nor_med)) %>% pull(mean_dipl) %>% unique()

  data_mod <- data %>% mutate(sd_dipl = sd_dipl, mean_dipl = mean_dipl, sds_median_dipl = (medianvst_nor_med - mean_dipl)/sd_dipl, sds_025_dipl = (low_quart - mean_dipl)/sd_dipl, sds_075_dipl = (up_quart-mean_dipl)/sd_dipl)

  return(data_mod)
}

### Create colour coding for estimation values ####
colour_code <- function(data) {
  data_col <- data %>% mutate(colour_arm = factor(ifelse(alteration_prob < 0.75, "bad", ifelse(alteration_prob > 0.9, "very_good", "good")), levels = c("bad", "good", "very_good"))) %>% group_by(chr) %>% mutate(min_prob = min(alteration_prob)) %>% ungroup() %>% mutate(colour_chr = factor(ifelse(min_prob < 0.75, "bad", ifelse(min_prob > 0.9, "very_good", "good")), levels = c("bad", "good", "very_good"), ordered = TRUE)) %>%
    select(-min_prob)
  return(data_col)
}

### Generate report for sample ###
gen_kar_list <- function(feat_tab_alt, sample_name, gender) {

  #extract only the chromosomes which have only one call
  alt_table <- feat_tab_alt %>% distinct(chr, alteration) %>% group_by(chr) %>% mutate(alt_types = length(unique(alteration))) %>% filter(alt_types == 1)

  #gain vector
  single_gain = alt_table$chr[alt_table$alteration == 1]
  gains = alt_table$chr[alt_table$alteration == 2] %>% rep(2) %>% c(single_gain) %>% sort()

  #del vector
  dels = alt_table$chr[alt_table$alteration == -1]

  #chromosome number
  chrom_n = 46 + length(gains) - length(dels)

  #call type
  type <- if (chrom_n == 46) {
    "diploid"
  } else if (chrom_n >= 51) {
    "high hyperdiploid"
  } else if (chrom_n >= 47 & chrom_n <= 50) {
    "low hyperdiploid"
  } else if (chrom_n <= 39 & chrom_n >= 32) {
    "low hypodiploid"
  } else if (chrom_n <= 31 & chrom_n >= 24) {
    "near haploid"
  } else if (chrom_n <= 43 & chrom_n >= 40) {
    "high hypodiploid"
  } else if (chrom_n <= 45 & chrom_n >= 44) {
    "near diploid"
  }

  if(length(dels) > 0) {
    dels_str <- dels %>% paste0(., "-")
  } else {
    dels_str <- c()
  }

  if(length(gains) > 0) {
    gains_str <- gains %>% paste0(., "+")
  } else {
    gains_str <- c()
  }

  alterations = paste0(c(gains_str, dels_str), collapse = ",")

  kar_table <- data.frame(sample = sample_name,
                          gender = factor(gender, levels = c("female", "male")),
                          chrom_n = chrom_n,
                          type = factor(type, levels = c("diploid", "high hyperdiploid", "low hyperdiploid", "low hypodiploid", "near haploid", "high hypodiploid", "near diploid")),
                          alterations = alterations, stringsAsFactors = FALSE)

  return(kar_table)
}
