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
          isUnion = TRUE) 
{
  if(1){
    
    if (is.null(genes) & is.null(samples)) {
      stop("Either provide Gene names or Sample names to group by.")
    }
    if (!is.null(genes) & !is.null(samples)) {
      stop("Either provide Gene names or Sample names to group by. Not both!")
    }
    if (is.null(clinicalData)) {
      print(clinicalData)
      message("Looking for clinical data in annoatation slot of MAF..")
      clinicalData = getClinicalData(x = maf)
      clinicalData = data.table::setDT(clinicalData)
      #print(clinicalData)
    }
    if (!"Tumor_Sample_Barcode" %in% colnames(clinicalData)) {
      print(colnames(clinicalData))
      stop("Column Tumo_Sample_Barcode not found in clinical data. Check column names and rename it to Tumo_Sample_Barcode if necessary.")
    }
    if (length(colnames(clinicalData)[colnames(clinicalData) %in% 
                                      time]) == 0) {
      print(colnames(clinicalData))
      stop(paste0(time, " not found in clinicalData. Use argument time to povide column name containing time to event."))
    }
    else {
      colnames(clinicalData)[colnames(clinicalData) %in% time] = "Time"
      #print("time, has been replaced by Time")
      #print(clinicalData)
    }
    if (length(colnames(clinicalData)[colnames(clinicalData) %in% 
                                      Status]) == 0) {
      print(colnames(clinicalData))
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
      print(sapply(genesTSB, FUN = length))
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
    print(clin.mut.dat)
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
	write.xlsx(clinicalData, paste0("data_output/",fn,".xlsx"))
  }
  
}
