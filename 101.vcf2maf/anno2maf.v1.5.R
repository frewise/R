ver=as.character("1.5")
update=c("1. Add tumor variant allele frequency",
         "2. remove unknown variants",
         "3. Add normal frequency filter",
         "4. Optionally select whether output filtered individuals",
         "5. Add NA test for conditions when there is not control sample (N.Freq is empty)",
         "6. Generate gene frequency file: gene.freq.txt")
library('getopt')

command=matrix(c("infile","i",1,"character",
				"pon","p",1,"character",
                 "tfrequpper","t",1,"character",
                 "tfreqlower","l",1,"character",
                 "tdp","d",1,"character",
                 "nfrequency","n",1,"character",
                 "unknown","u",0,"logical",
                 "tmp","m",0,"logical",
                 "outfile","o",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if (!is.null(args$help) || is.null(args$infile) || is.null(args$pon) || is.null(args$tfreqlower) || is.null(args$tfreqlower) || is.null(args$tdp) || is.null(args$outfile)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}else{
	logfile<-paste(Sys.Date(),".log",sep="")
	cat(paste(getopt(command), "\n"))
	write.table(paste(getopt(command), "\n"),file=logfile,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE)
  pon=args$pon
  tfrequpper=args$tfrequpper
  tfreqlower=args$tfreqlower
  tdp=as.integer(args$tdp)#这里很奇怪，如果没有这一步转换，结果异常
  nfrequency=ifelse(!is.null(args$nfrequency),args$nfrequency,1)#如果指定则使用指定值，如果不指定则取1
  infile=args$infile
  outfile=args$outfile
  cosmic<-c()
  gene_freq_df=read.table(text="Gene Freq",header = TRUE)
  #输出文件格式
  df <- data.frame(Chr=character(),
                 Start=character(), 
                 End=character(), 
                 Ref=character(), 
                 Alt=character(), 
                 Gene.refGene=character(), 
                 GeneDetail.refGene=character(), 
                 ExonicFunc.refGene=character(), 
                 AAChange.refGene=character(), 
                 Tumor_Sample_Barcode=character(), 
                 Func.refGene=character(), 
                 t_vaf=character(),
                 stringsAsFactors=FALSE) 
  write.table(df,file=outfile,quote=FALSE,row.names = FALSE,sep="\t")

  #TMB
  tmb_out <- data.frame(patient=character(),
  						mtCNT=character(),
  						tmb=character(),
  						stringsAsFactors=FALSE)
  write.table(tmb_out,file="tmb.txt",quote=FALSE,row.names = FALSE,col.names = TRUE,sep="\t")


  for(i in list.files(".",pattern=paste("*.",infile,sep=""))){
  	cat("Processing",i,"...\n")
  	write.table(paste("Processing",i,"...\n"),file=logfile,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE)
	  x<-read.table(i,sep="\t", quote="\"",header=TRUE, fill=TRUE, stringsAsFactors = FALSE)
	  cat("  There are ",dim(x)[1],"mutations before filter","\n")
	  write.table(paste("  There are ",dim(x)[1],"mutations before filter"),file=logfile,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE)
	  x<-x[x$PoN<=pon,]
	  cat("  There are ",dim(x)[1],"mutations after filter PoN","\n")
      write.table(paste("  There are ",dim(x)[1],"mutations after filter PoN"),file=logfile,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE)
  	  x<-x[x$T.Freq>=tfreqlower & x$T.Freq<=tfrequpper,]
  	  cat("  There are ",dim(x)[1],"mutations after filter T.Freq","\n")
  	  write.table(paste("  There are ",dim(x)[1],"mutations after filter T.Freq"),file=logfile,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE)
	  x<-x[x$T.DP>=tdp,]
	  cat("  There are ",dim(x)[1],"mutations after filter T.DP","\n")
	  write.table(paste("  There are ",dim(x)[1],"mutations after filter T.DP"),file=logfile,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE)
    x<-x[x$N.Freq<=nfrequency|is.na(x$N.Freq),]
    cat("  There are ",dim(x)[1],"mutations after filter N.Freq","\n")
    write.table(paste("  There are ",dim(x)[1],"mutations after filter N.Freq"),file=logfile,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE)
    # Remove unknown variants
    if(!is.null(args$unknown)){
      x<-subset(x,ExonicFunc.refGene!="unknown")
      cat("  There are ",dim(x)[1],"mutations after filter unknown variants","\n")
      write.table(paste("  There are ",dim(x)[1],"mutations after filter unknown variants"),file=logfile,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE)
    }




	  #cat(x$T.DP,"\n")
	  cosmic<-c(cosmic,x[grepl("COSM",x$cosmic70),]$Gene.refGene)
	  #freq<-as.array(x$T.Freq)
	  #cnt<-length(as.array(x$T.Freq))
	  patient <- strsplit(i,"\\.")[[1]][1]
	  x$Tumor_Sample_Barcode=rep(patient,nrow(x))
	  y<-x[,c("Chr","Start","End","Ref","Alt","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene","Tumor_Sample_Barcode","Func.refGene","T.Freq")]
	  #write.table(y,file=paste(patient,"cosmic.anno.txt",sep="."),quote=FALSE,row.names = FALSE,sep="\t")
	  #write.table(colnames(y),file=outfile,quote=FALSE,row.names = FALSE,sep="\t")
	  #cat(paste(patient,dim(y)[1]/1.66,sep="\t"),"\n")
    if(!is.null(args$tmp)){
      write.table(x,file=paste(patient,".filtered.xls.txt",sep=""),quote=FALSE,row.names=FALSE,col.names = FALSE,sep="\t")
    }
	  write.table(y,file=outfile,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
	  #rm(list=ls())
    #TMB
    mtCNT=dim(x)[1]
    tmb=mtCNT/1.66
    tmb_out=cbind(patient,mtCNT,tmb)
    write.table(tmb_out,file="tmb.txt",append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE,sep="\t")


    if(dim(x)[1]>0){
    	k=aggregate(x[,"T.Freq"],list(Gene.freq=x$Gene.refGene),mean)
    	colnames(k)=c("Gene",patient)
    	gene_freq_df=merge(gene_freq_df,k,by = "Gene",all=TRUE)
    }
    # Get gene mutation frequency
    
	}
  cosmic<-unique(cosmic)
  write.table(cosmic,file="cosmic_list.txt",quote=FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
  # gene frequency
  gene_freq_df=gene_freq_df[-2]
  write.table(gene_freq_df,file="gene.freq.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")


 
}