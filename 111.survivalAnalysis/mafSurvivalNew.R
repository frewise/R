library(maftools)

mafReheader <- function(maf){
	library(stringr)
	library(dplyr)
	# --------------------------- column replacement
	colnames(maf)[1]="Chromosome"
	colnames(maf)[2]="Start_Position"
	colnames(maf)[3]="End_Position"
	colnames(maf)[4]="Reference_Allele"
	colnames(maf)[5]="Tumor_Seq_Allele2"
	colnames(maf)[6]="Hugo_Symbol"
	colnames(maf)[8]="Variant_Classification"
	colnames(maf)[11]="Variant_Type"

	# --------------------------- variant classification replacement

	# maf[maf$Variant_Classification=="nonsynonymous SNV",]$Variant_Classification="Missense_Mutation"
	# maf[maf$Variant_Classification=="frameshift deletion",]$Variant_Classification="Frame_Shift_Del"
	# maf[maf$Variant_Classification=="frameshift insertion",]$Variant_Classification="Frame_Shift_Ins"
	# maf[maf$Variant_Classification=="nonframeshift deletion",]$Variant_Classification="In_Frame_Del"
	# maf[maf$Variant_Classification=="nonframeshift insertion",]$Variant_Classification="In_Frame_Ins"
	# maf[maf$Variant_Classification=="stopgain",]$Variant_Classification="Nonsense_Mutation"
	# maf[maf$Variant_Classification=="stoploss",]$Variant_Classification="Nonstop_Mutation"
	# maf[maf$Variant_Classification==".",]$Variant_Classification="Splice_Site"

	# maf.txt$Variant_Classification <- str_replace(maf.txt$Variant_Classification, "stoploss","Nonstop_Mutation")

	maf <- filter(maf, Variant_Classification != "synonymous SNV")



	maf$Variant_Classification <- maf$Variant_Classification %>%
	  str_replace("\\bnonsynonymous SNV\\b","Missense_Mutation") %>% 
	  str_replace("\\bnonframeshift deletion\\b","In_Frame_Del") %>% 
	  str_replace("\\bnonframeshift insertion\\b","In_Frame_Ins") %>% 
	  str_replace("\\bframeshift deletion\\b","Frame_Shift_Del") %>%
	  str_replace("\\bframeshift insertion\\b","Frame_Shift_Ins") %>%
	  str_replace("\\bstopgain\\b","Nonsense_Mutation") %>% 
	  str_replace("\\bstoploss\\b","Nonstop_Mutation") %>% 
	  str_replace('\\.',"Splice_Site")


	return(maf)
}

mafSurvival.new2 <- function (maf, genes = NULL, samples = NULL, clinicalData = NULL, 
          time = "Time", Status = "Status", groupNames = c("Mutant", "WT"), 
          textSize = 12, fn = NULL, isPrint = TRUE,
          width = 12, height = 10,
          xlab = "Survival", ylab = "Time", pvalue = 0.05,
          legendTitle = "Genotype",isRiskTable = TRUE,
          isUnion = TRUE,
          verbose = TRUE) 
{
  if(1){
    
    if (is.null(genes) & is.null(samples)) {
      stop("Either provide Gene names or Sample names to group by.")
    }
    if (!is.null(genes) & !is.null(samples)) {
      stop("Either provide Gene names or Sample names to group by. Not both!")
    }
    if (is.null(clinicalData)) {
      if(verbose) print(clinicalData)
      message("Looking for clinical data in annoatation slot of MAF..")
      clinicalData = getClinicalData(x = maf)
      clinicalData = data.table::setDT(clinicalData)
      #print(clinicalData)
    }
    if (!"Tumor_Sample_Barcode" %in% colnames(clinicalData)) {
      if(verbose) print(colnames(clinicalData))
      stop("Column Tumo_Sample_Barcode not found in clinical data. Check column names and rename it to Tumo_Sample_Barcode if necessary.")
    }
    if (length(colnames(clinicalData)[colnames(clinicalData) %in% 
                                      time]) == 0) {
      if(verbose) print(colnames(clinicalData))
      stop(paste0(time, " not found in clinicalData. Use argument time to povide column name containing time to event."))
    }
    else {
      colnames(clinicalData)[colnames(clinicalData) %in% time] = "Time"
      #print("time, has been replaced by Time")
      #print(clinicalData)
    }
    if (length(colnames(clinicalData)[colnames(clinicalData) %in% 
                                      Status]) == 0) {
      if(verbose) print(colnames(clinicalData))
      stop(paste0(Status, " not found in clinicalData. Use argument Status to povide column name containing events (Dead or Alive)."))
    }
    else {
      colnames(clinicalData)[colnames(clinicalData) %in% Status] = "Status"
      #print("status, has been replaced by Status")
      #print(clinicalData)
    }
    if (!is.null(genes)) {
      genesTSB = genesToBarcodes(maf = maf, genes = genes, 
                                 justNames = TRUE)
      genesTSB = genesTSB[sapply(genesTSB, FUN = function(x) length(x) != 
                                   0)]
      barcodeWmut = c()
      message("Number of mutated samples for given genes: ")
      if(verbose) print(sapply(genesTSB, FUN = length))
      barcodeWmut = unique(as.character(unlist(genesTSB)))
      genesMissing = genes[!genes %in% names(genesTSB)]
      if (length(genesMissing) > 0) {
        genes = genes[!genes %in% genesMissing]
        genesMissing = paste(genesMissing, collapse = ", ")
        message(paste0("genes ", genesMissing, " does not seeem to be mutated. Removing them."))
      }
      if (length(genes) == 0) {
        stop("None of the given genes are mutated!")
      }
      else {
        genes = paste(genes, collapse = ", ")
      }
      if(isUnion==TRUE){
        genesTSB = unique(as.character(unlist(genesTSB)))
        
      }else{
        genesTSB = Reduce(intersect, genesTSB)
      }
      
    }
    else {
      genesTSB = samples
    }
    data.table::setDT(clinicalData)
    clinicalData$Time = suppressWarnings(as.numeric(as.character(clinicalData$Time)))
    clinicalData$Status = suppressWarnings(as.integer(as.character(clinicalData$Status)))
    clinicalData=mutate(clinicalData, Group="SingleMut")
    clinicalData[clinicalData$Tumor_Sample_Barcode %in% genesTSB,]$Group = groupNames[1]
    clinicalData[!(clinicalData$Tumor_Sample_Barcode %in% barcodeWmut),]$Group = groupNames[2]
    data.table::setDT(clinicalData)
    #clinicalData$Group = ifelse(test = clinicalData$Tumor_Sample_Barcode %in% 
    #                              genesTSB, yes = groupNames[1], no = groupNames[2])
    clin.mut.dat = clinicalData[, .(medianTime = median(Time, 
                                                        na.rm = TRUE), N = .N), Group][order(Group)]
    #print(clinicalData)
    message("Median survival..")
    if(verbose) print(clin.mut.dat)
    clinicalData$Time = ifelse(test = is.infinite(clinicalData$Time), 
                               yes = 0, no = clinicalData$Time)
    surv.km = survival::survfit(formula = survival::Surv(time = Time, 
                                                         event = Status) ~ Group, data = clinicalData, conf.type = "log-log")
    #print(surv.km)
    strata=names(surv.km$strata)
    strata=gsub(pattern = "Group=",
                replacement = "",
                x = strata)
    names(surv.km$strata)=strata
    
    
    
    
    
    
    res = summary(surv.km)
    surv.diff = survival::survdiff(formula = survival::Surv(time = Time, 
                                                            event = Status) ~ Group, data = clinicalData)
    
    surv.diff.pval = round(1 - pchisq(surv.diff$chisq, length(surv.diff$n) - 
                                        1), digits = 5)
    surv.dat = data.frame(Group = res$strata, Time = res$time, 
                          survProb = res$surv, survUp = res$upper, survLower = res$lower)
    surv.dat$Group = gsub(pattern = "Group=", replacement = "", 
                          x = surv.dat$Group)
    if(surv.diff.pval <=pvalue){
      surv.gg = survminer::ggsurvplot(fit = surv.km, 
                                      data = clinicalData,#奇怪，这里的surv.km已经包含data信息，但data设置为null就报错
                                      pval = FALSE, 
                                      xlab = xlab,
                                      ylab = ylab,
                                      legend.title = legendTitle,
                                      break.x.by = 5, 
                                      palette = ezfun::msk_palette("contrast"), 
                                      censor = TRUE,
                                      risk.table = isRiskTable,
                                      risk.table.y.text = TRUE)
      
      surv.gg$plot = surv.gg$plot + 
        ggtitle(label = paste0(paste(genes, collapse="_"), " ", groupNames[1], " v/s ", groupNames[2]), 
                subtitle = paste0("P-value: ", surv.diff.pval))
      
      
      if (!is.null(fn)) {
        ggsave(paste0(fn, ".pdf"), plot = print(surv.gg), height = height, width = width)
      }
      if (isPrint){
        print(surv.gg)
      }
    }
    #print(clin.mut.dat[[1]])
    #print(clin.mut.dat[[2]])
    #print(clin.mut.dat[[3]])
    data.frame(gene=genes,pvalue=surv.diff.pval,
               Mut_median=clin.mut.dat$medianTime[clin.mut.dat$Group=="Mutant"],
               WT_median=clin.mut.dat$medianTime[clin.mut.dat$Group=="WT"],
               Mut_case=clin.mut.dat$N[clin.mut.dat$Group=="Mutant"],
               WT_case=clin.mut.dat$N[clin.mut.dat$Group=="WT"])
    #write.table(clinicalData,file="")
    #  if (!require("openxlsx")) install.packages("openxlsx")
    #if (!dir.exists("data_output/")) {
    #	dir.create("data_output/")
    #}
    #write.xlsx(clinicalData, paste0(fn,".xlsx"))
  }
  
}

fisherCorrection = function(fc){
  #Adding a tiny value to zero/Inf odds ratios (zeros become Inf's when log converted)
  fc = as.data.frame(fc)
  fc$or = ifelse(test = is.infinite(fc$or), yes = fc$ci.up, no = fc$or )
  fc$or = ifelse(test = fc$or == 0, yes = fc$ci.low, no = fc$or)
  
  fc$ci.up = ifelse(test = fc$ci.up == 0, yes = fc$or, no = fc$ci.up)
  fc$ci.low = ifelse(test = is.infinite(fc$ci.low), yes = fc$ci.up, no = fc$ci.low)
  
  return(data.table::data.table(fc))
  
}


forestPlot.2<-function (mafCompareRes, pVal = 0.05, fdr = NULL, color = NULL, 
          geneFontSize = 1.2, titleSize = 1.2, lineWidth = 2.2, file = NULL, 
          width = 5, height = 6) 
{
  res = mafCompareRes$results
  if (is.null(fdr)) {
    m.sigs = res[pval < pVal]
  }
  else {
    m.sigs = res[adjPval < fdr]
  }
  m1Name = mafCompareRes$SampleSummary[1, Cohort]
  m2Name = mafCompareRes$SampleSummary[2, Cohort]
  m1.sampleSize = mafCompareRes$SampleSummary[1, SampleSize]
  m2.sampleSize = mafCompareRes$SampleSummary[2, SampleSize]
  if (nrow(m.sigs) < 1) {
    stop("No differetially mutated genes found !")
  }
  m.sigs = fisherCorrection(fc = m.sigs)
  m.sigs$Hugo_Symbol = factor(x = m.sigs$Hugo_Symbol, levels = rev(m.sigs$Hugo_Symbol))
  m.sigs[, `:=`(log10OR, log10(or))]
  m.sigs[, `:=`(log10OR_high, log10(ci.up))]
  m.sigs[, `:=`(log10OR_low, log10(ci.low))]
  m.sigs$flow = ifelse(test = m.sigs$log10OR < 0, yes = m2Name, 
                       no = m1Name)
  m.sigs$statRight = paste(m2Name, ":", m.sigs[, 3, with = FALSE][[1]], 
                           sep = "")
  m.sigs$statLeft = paste(m1Name, ":", m.sigs[, 2, with = FALSE][[1]], 
                          sep = "")
  m.sigs = m.sigs[order(pval, decreasing = TRUE)]
  lim = max(abs(c(log10(m.sigs$ci.up), log10(m.sigs$ci.low)))) + 
    1
  lim = round(lim, digits = 2)
  if (!is.null(color)) {
    color = color
    names(color) = c(m1Name, m2Name)
  }
  else {
    color = c("royalblue", "maroon")
    names(color) = c(m1Name, m2Name)
  }
  layout(mat = matrix(c(1, 2, 3, 4, 5, 5, 5, 5), byrow = TRUE, 
                      ncol = 4, nrow = 2), widths = c(4, 1, 1), heights = c(6, 
                                                                            1.2))
  par(mar = c(3, 1, 3, 5))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(-lim, 
                                                      lim), axes = FALSE, pch = NA, xlab = "", ylab = "", 
       ylim = c(0.5, nrow(m.sigs)))
  segments(x0 = m.sigs$log10OR_low, y0 = 1:nrow(m.sigs), x1 = m.sigs$log10OR_high, 
           y1 = 1:nrow(m.sigs), lwd = lineWidth, color[m.sigs$flow])
  points(m.sigs$log10OR, 1:nrow(m.sigs), pch = 16, cex = 0.5 * 
           (lineWidth))
  axis(side = 1, at = c(-lim, 0, lim), lwd = 2.2, font = 2, 
       pos = 0.5, cex.axis = 1.3)
  segments(x0 = 0, y0 = 0.5, x1 = 0, y1 = nrow(m.sigs) + 0.2, 
           col = "gray70", lwd = 2, lty = 2)
  mtext(text = m.sigs$Hugo_Symbol, side = 4, line = 0.2, at = 1:nrow(m.sigs), 
        font = 4, las = 2, cex = geneFontSize, adj = 0)
  mtitle = paste(m2Name, " (n = ", m2.sampleSize, ")", " v/s ", 
                 m1Name, " (n = ", m1.sampleSize, ")", sep = "")
  title(main = mtitle, font = 2, adj = 0, cex.main = titleSize)
  mtext(text = "Log odds ratio", side = 1, line = 2, font = 2, 
        cex = 0.7 * (titleSize))
  par(mar = c(3, 0, 3, 0))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(0, 1), 
       axes = FALSE, pch = NA, xlab = "", ylab = "", ylim = c(0.5, 
                                                              nrow(m.sigs)))
  text(x = 0.5, y = 1:nrow(m.sigs), labels = unlist(m.sigs[, 
                                                                      2]), adj = 0, font = 2, cex = 1.4 * (geneFontSize))
  title(main = m1Name, cex.main = titleSize)
  par(mar = c(3, 0, 3, 0))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(0, 1), 
       axes = FALSE, pch = NA, xlab = "", ylab = "", ylim = c(0.5, 
                                                              nrow(m.sigs)))
  text(x = 0.5, y = 1:nrow(m.sigs), labels = unlist(m.sigs[, 
                                                                      3]), adj = 0, font = 2, cex = 1.4 * (geneFontSize))
  title(main = m2Name, cex.main = titleSize)
  m.sigs$significance = ifelse(test = as.numeric(m.sigs$pval) < 
                                 0.001, yes = "***", no = ifelse(test = as.numeric(m.sigs$pval) < 
                                                                   0.01, yes = "**", no = ifelse(test = as.numeric(m.sigs$pval) < 
                                                                                                   0.05, yes = "*", no = "NS")))
  par(mar = c(3, 0, 3, 0))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(0, 1), 
       axes = FALSE, pch = NA, xlab = "", ylab = "", ylim = c(0.5, 
                                                              nrow(m.sigs)))
  text(x = 0.5, y = 1:nrow(m.sigs), labels = m.sigs$significance, 
       adj = 0, font = 2, cex = 1.4 * (geneFontSize))
  title(main = "p-value", cex.main = titleSize)
  plot.new()
  par(mar = c(0, 0, 0, 0))
  legend(x = "top", legend = names(color), lwd = lineWidth, 
         col = color[c(m1Name, m2Name)], border = NA, bty = "n", 
         cex = 1.1 * (titleSize), horiz = TRUE, text.font = 2, 
         xpd = TRUE, pch = 16)
  if (!is.null(file)) {
    cowplot::save_plot(filename = paste(file, "pdf", sep = "."), 
                       plot = gg.fp, base_height = height, base_width = width, 
                       paper = "special", bg = "white")
  }
}

validateMaf.2<-function (maf, rdup = TRUE, isTCGA = isTCGA, chatty = TRUE) 
{
  required.fields = c("Hugo_Symbol", "Chromosome", "Start_Position", 
                      "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                      "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode")
  for (i in 1:length(required.fields)) {
    colId = suppressWarnings(grep(pattern = paste("^", required.fields[i], 
                                                  "$", sep = ""), x = colnames(maf), ignore.case = TRUE))
    if (length(colId) > 0) {
      colnames(maf)[colId] = required.fields[i]
    }
  }
  missing.fileds = required.fields[!required.fields %in% colnames(maf)]
  if (length(missing.fileds) > 0) {
    missing.fileds = paste(missing.fileds[1], sep = ",", 
                           collapse = ", ")
    stop(paste("missing required fields from MAF:", missing.fileds))
  }
  if (rdup) {
    maf = maf[, `:=`(variantId, paste(Chromosome, Start_Position, 
                                      Tumor_Sample_Barcode, sep = ":"))]
    if (nrow(maf[duplicated(variantId)]) > 0) {
      if (chatty) {
        message("NOTE: Removed ", nrow(maf[duplicated(variantId)]), 
                " duplicated variants")
      }
      maf = maf[!duplicated(variantId)]
    }
    maf[, `:=`(variantId, NULL)]
  }
  if (nrow(maf[Hugo_Symbol %in% ""]) > 0) {
    if (chatty) {
      message("NOTE: Found ", nrow(maf[Hugo_Symbol %in% 
                                         ""]), " variants with no Gene Symbols.")
      print(maf[Hugo_Symbol %in% "", required.fields, 
                with = FALSE])
      message("Annotating them as 'UnknownGene' for convenience")
    }
    maf$Hugo_Symbol = ifelse(test = maf$Hugo_Symbol == "", 
                             yes = "UnknownGene", no = maf$Hugo_Symbol)
  }
  if (nrow(maf[is.na(Hugo_Symbol)]) > 0) {
    if (chatty) {
      message("NOTE: Found ", nrow(maf[is.na(Hugo_Symbol) > 
                                         0]), " variants with no Gene Symbols.")
      print(maf[is.na(Hugo_Symbol), required.fields, with = FALSE])
      message("Annotating them as 'UnknownGene' for convenience")
    }
    maf$Hugo_Symbol = ifelse(test = is.na(maf$Hugo_Symbol), 
                             yes = "UnknownGene", no = maf$Hugo_Symbol)
  }
  if (isTCGA) {
    maf$Tumor_Sample_Barcode = substr(x = maf$Tumor_Sample_Barcode, 
                                      start = 1, stop = 12)
  }
  silent = c("3'UTR", "5'UTR", "3'Flank", "Targeted_Region", 
             "Silent", "Intron", "RNA", "IGR", "Splice_Region", "5'Flank", 
             "lincRNA", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", 
             "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del")
  vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
                   "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", 
                   "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")
  vt = c("SNP", "DNP", "TNP", "ONP", "INS", "DEL")
  maf.vcs = unique(as.character(maf[, Variant_Classification]))
  maf.vts = unique(as.character(maf[, Variant_Type]))
  if (length(maf.vcs[!maf.vcs %in% c(silent, vc.nonSilent)] > 
             0)) {
    message("NOTE: Non MAF specific values in Variant_Classification column:")
    #print(maf.vcs[!maf.vcs %in% c(silent, vc.nonSilent)])
  }
  if (length(maf.vts[!maf.vts %in% vt] > 0)) {
    message("NOTE: Non MAF specific values in Variant_Type column:")
    #print(maf.vts[!maf.vts %in% vt])
  }
  maf$Tumor_Sample_Barcode = as.factor(maf$Tumor_Sample_Barcode)
  return(maf)
}



read.maf.2<-function (maf, clinicalData = NULL, removeDuplicatedVariants = TRUE, 
          useAll = TRUE, gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL, 
          gisticDelGenesFile = NULL, gisticScoresFile = NULL, cnLevel = "all", 
          cnTable = NULL, isTCGA = FALSE, vc_nonSyn = NULL, verbose = TRUE) 
{
  if (is.data.frame(x = maf)) {
    maf = data.table::setDT(maf)
  }
  else {
    if (verbose) {
      message("reading maf..")
    }
    if (as.logical(length(grep(pattern = "gz$", x = maf, 
                               fixed = FALSE)))) {
      if (Sys.info()[["sysname"]] == "Windows") {
        maf.gz = gzfile(description = maf, open = "r")
        suppressWarnings(maf <- data.table::as.data.table(read.csv(file = maf.gz, 
                                                                   header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                                                                   comment.char = "#")))
        close(maf.gz)
      }
      else {
        maf = suppressWarnings(data.table::fread(input = paste("zcat <", 
                                                               maf), sep = "\t", stringsAsFactors = FALSE, 
                                                 verbose = FALSE, data.table = TRUE, showProgress = TRUE, 
                                                 header = TRUE, fill = TRUE, skip = "Hugo_Symbol"))
      }
    }
    else {
      suppressWarnings(maf <- data.table::fread(input = maf, 
                                                sep = "\t", stringsAsFactors = FALSE, verbose = FALSE, 
                                                data.table = TRUE, showProgress = TRUE, header = TRUE, 
                                                fill = TRUE, skip = "Hugo_Symbol"))
    }
  }
  maf = validateMaf.2(maf = maf, isTCGA = isTCGA, rdup = removeDuplicatedVariants, 
                    chatty = verbose)
  if (!useAll) {
    message("--Using only `Somatic` variants from Mutation_Status. Set useAll = TRUE to include everything.")
    if (length(colnames(maf)[colnames(x = maf) %in% "Mutation_Status"]) > 
        0) {
      maf = maf[Mutation_Status %in% "Somatic"]
      if (nrow(maf) == 0) {
        stop("No more Somatic mutations left after filtering for Mutation_Status! Maybe set useAll to TRUE ?")
      }
    }
    else {
      message("---Oops! Mutation_Status not found. Assuming all variants are Somatic and validated.")
    }
  }
  if (is.null(vc_nonSyn)) {
    vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", 
                     "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", 
                     "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", 
                     "Missense_Mutation")
  }
  else {
    vc.nonSilent = vc_nonSyn
  }
  maf.silent = maf[!Variant_Classification %in% vc.nonSilent]
  if (nrow(maf.silent) > 0) {
    maf.silent.vc = maf.silent[, .N, .(Tumor_Sample_Barcode, 
                                       Variant_Classification)]
    maf.silent.vc.cast = data.table::dcast(data = maf.silent.vc, 
                                           formula = Tumor_Sample_Barcode ~ Variant_Classification, 
                                           fill = 0, value.var = "N")
    summary.silent = data.table::data.table(ID = c("Samples", 
                                                   colnames(maf.silent.vc.cast)[2:ncol(maf.silent.vc.cast)]), 
                                            N = c(nrow(maf.silent.vc.cast), colSums(maf.silent.vc.cast[, 
                                                                                                       2:ncol(maf.silent.vc.cast), with = FALSE])))
    maf = maf[Variant_Classification %in% vc.nonSilent]
    if (verbose) {
      message(paste0("silent variants: ", nrow(maf.silent)))
      print(summary.silent)
    }
  }
  if (nrow(maf) == 0) {
    stop("No non-synonymous mutations found\nCheck `vc_nonSyn`` argumet in `read.maf` for details")
  }
  if (!is.null(gisticAllLesionsFile)) {
    gisticIp = readGistic(gisticAllLesionsFile = gisticAllLesionsFile, 
                          gisticAmpGenesFile = gisticAmpGenesFile, gisticDelGenesFile = gisticDelGenesFile, 
                          isTCGA = isTCGA, gisticScoresFile = gisticScoresFile, 
                          cnLevel = cnLevel)
    gisticIp = gisticIp@data
    suppressWarnings(gisticIp[, `:=`(id, paste(Hugo_Symbol, 
                                               Tumor_Sample_Barcode, sep = ":"))])
    gisticIp = gisticIp[!duplicated(id)]
    gisticIp[, `:=`(id, NULL)]
    maf = rbind(maf, gisticIp, fill = TRUE)
    maf$Tumor_Sample_barcode = factor(x = maf$Tumor_Sample_barcode, 
                                      levels = unique(c(levels(maf$Tumor_Sample_barcode), 
                                                        unique(as.character(gisticIp$Tumor_Sample_barcode)))))
  }
  else if (!is.null(cnTable)) {
    if (verbose) {
      message("Processing copy number data..")
    }
    if (is.data.frame(cnTable)) {
      cnDat = data.table::setDT(cnTable)
    }
    else {
      cnDat = data.table::fread(input = cnTable, sep = "\t", 
                                stringsAsFactors = FALSE, header = TRUE, colClasses = "character")
    }
    colnames(cnDat) = c("Hugo_Symbol", "Tumor_Sample_Barcode", 
                        "Variant_Classification")
    if (isTCGA) {
      cnDat[, `:=`(Tumor_Sample_Barcode, substr(x = cnDat$Tumor_Sample_Barcode, 
                                                start = 1, stop = 12))]
    }
    cnDat$Variant_Type = "CNV"
    suppressWarnings(cnDat[, `:=`(id, paste(Hugo_Symbol, 
                                            Tumor_Sample_Barcode, sep = ":"))])
    cnDat = cnDat[!duplicated(id)]
    cnDat[, `:=`(id, NULL)]
    maf = rbind(maf, cnDat, fill = TRUE)
    maf$Tumor_Sample_barcode = factor(x = maf$Tumor_Sample_barcode, 
                                      levels = unique(c(levels(maf$Tumor_Sample_barcode), 
                                                        unique(as.character(cnDat$Tumor_Sample_barcode)))))
  }
  maf$Variant_Classification = as.factor(as.character(maf$Variant_Classification))
  maf$Variant_Type = as.factor(as.character(maf$Variant_Type))
  if (verbose) {
    message("Summarizing..")
  }
  mafSummary = maftools:::summarizeMaf(maf = maf, anno = clinicalData, 
                            chatty = verbose)
  m = maftools:::MAF(data = maf, variants.per.sample = mafSummary$variants.per.sample, 
          variant.type.summary = mafSummary$variant.type.summary, 
          variant.classification.summary = mafSummary$variant.classification.summary, 
          gene.summary = mafSummary$gene.summary, summary = mafSummary$summary, 
          maf.silent = maf.silent, clinical.data = mafSummary$sample.anno)
  if (verbose) {
    message("Done !")
  }
  return(m)
}
