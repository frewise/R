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