ver=as.character("1.4")
updates=paste('1. Auto resize the figure according to sampleSize','2. Add tiff output', '3. Add "Nonstop" category in variant type legend',sep="\n")

library('getopt')

command=matrix(c("maf","m",1,"character",
                 "clinical","c",1,"character",
                 "outfile","o",1,"character",
                 "group","g",1,"character",
                 "math","y",0,"logical",
                 "help","h",0,"logical",
                 "update","u",0,"logical"),byrow=T,ncol=4)
args=getopt(command)
if (!is.null(args$update)) {
  cat("Version:",ver,"\n","Updates:","\n",updates,"\n",sep="")
  q()
}
if (!is.null(args$help) || is.null(args$maf) || is.null(args$outfile)){
  cat(paste(getopt(command, usage = T), "\n", "Version:",ver,sep=""))
  q()
}
library(maftools)
group<-""
if(!is.null(args$maf) && !is.null(args$outfile)){
  maf=args$maf
  outfile=args$outfile
  if(!is.null(args$clinical)){
    #group<-ifelse(!is.null(args$group),unlist(strsplit(args$group,',')),group)
    if(!is.null(args$group)){
      group<-unlist(strsplit(args$group,','))
        }else{
          group<-group
        }
    cat("Clinical features used:",group,"\n")
    anno=read.table(args$clinical,header = TRUE)
    laml<-annovarToMaf(maf,sampleAnno=anno,refBuild="hg19",sep="\t",MAFobj = TRUE)
  }else{
    laml<-annovarToMaf(maf,refBuild="hg19",sep="\t",MAFobj = TRUE)
  }
}




#color
col = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Mutation','Fusion')




summary<-getGeneSummary(laml)
summary<-summary[summary$AlteredSamples>=2,]
sampleSize=dim(getSampleSummary(laml))[1]
width=round(sampleSize/5)

genes<-as.vector(summary$Hugo_Symbol)
if(file.exists("cosmic_list.txt")){
	genes.cosmic<-read.table("cosmic_list.txt",header=FALSE)#cosmic_list.txt由filter.R生成，记录所有含有cosmic ID的基因
	genes.cosmic<-as.array(genes.cosmic$V1)
	genes<-append(as.vector(genes.cosmic),genes)
  print("There cosmic_list.txt file")

  #
  #pdf output
  #

  #top 40 genes, including genes with cosmic ID
  pdf(file=paste(outfile,"pdf",sep="."),width=width,height = 10)
  #oncoplot(maf = laml, clinicalFeatures = c('FAB_classification','Overall_Survival_Status'),top = 40, colors=col,showTumorSampleBarcodes = TRUE,fontSize = 8)
  oncoplot(maf = laml, clinicalFeatures = if(nchar(group)>1) group else NULL,top = 40, SampleNamefontSize = 6,colors=col,showTumorSampleBarcodes = TRUE,fontSize = 8)
  dev.off()
  #show only genes with cosmic ID
  pdf(file=paste(outfile,"cosmic.pdf",sep="."),width=width,height = 13)
  oncoplot(maf = laml, clinicalFeatures = if(nchar(group)>1) group else NULL, genes=genes, SampleNamefontSize = 6,colors=col,showTumorSampleBarcodes = TRUE,legendFontSize = 8,
           annotationFontSize = 8, annotationTitleFontSize = 8,fontSize = 8)
  dev.off()
  #mutation statistics
  pdf(file=paste(outfile,"vis.pdf",sep="."))
  plotmafSummary(maf = laml, color=col, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()

  #
  #tiff output
  #

  #top 40 genes, including genes with cosmic ID
  tiff(file=paste(outfile,"tiff",sep="."),units="in", width=width,height = 10, res=200)
  #oncoplot(maf = laml, clinicalFeatures = c('FAB_classification','Overall_Survival_Status'),top = 40, colors=col,showTumorSampleBarcodes = TRUE,fontSize = 8)
  oncoplot(maf = laml, clinicalFeatures = if(nchar(group)>1) group else NULL,top = 40, SampleNamefontSize = 6,colors=col,showTumorSampleBarcodes = TRUE,fontSize = 8)
  dev.off()
  #show only genes with cosmic ID
  tiff(file=paste(outfile,"cosmic.tiff",sep="."),units="in", width=width,height = 13, res=200)
  oncoplot(maf = laml, clinicalFeatures = if(nchar(group)>1) group else NULL, genes=genes, SampleNamefontSize = 6,colors=col,showTumorSampleBarcodes = TRUE,legendFontSize = 8,
           annotationFontSize = 8, annotationTitleFontSize = 8,fontSize = 8)
  dev.off()
  #mutation statistics
  tiff(file=paste(outfile,"vis.tiff",sep="."))
  plotmafSummary(maf = laml, color=col, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()



}else{
  #There is not cosmic files
  genes<-as.vector(summary$Hugo_Symbol)
  #landscape

  #
  #pdf output
  #

  pdf(file=paste(outfile,"pdf",sep="."),width=width,height = 10)
  oncoplot(maf = laml, clinicalFeatures = if(nchar(group)>1) group else NULL, top = 40, SampleNamefontSize = 6,colors=col,showTumorSampleBarcodes = TRUE,fontSize = 8)
  dev.off()
  #mutation statistics
  pdf(file=paste(outfile,"vis.pdf",sep="."))
  plotmafSummary(maf = laml, color=col, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()


  #
  #tiff output
  #
  tiff(file=paste(outfile,"tiff",sep="."),units="in",width=width,height = 10, res=200)
  oncoplot(maf = laml, clinicalFeatures = if(nchar(group)>1) group else NULL, top = 40, SampleNamefontSize = 6,colors=col,showTumorSampleBarcodes = TRUE,fontSize = 8)
  dev.off()
  #mutation statistics
  tiff(file=paste(outfile,"vis.pdf",sep="."))
  plotmafSummary(maf = laml, color=col, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
}

if(!is.null(args$math)){
  #Hetergeneity calculation: MATH
  dir.create("math")
  mathoutputdir<-file.path(getwd(),"math")
  setwd(mathoutputdir)
  mathfile<-"math.txt"

  write.table(paste("Sample","MATH",sep="\t"),file = mathfile,quote = FALSE,sep="\t",row.names = FALSE,col.names = FALSE)
  for(i in getSampleSummary(laml)$Tumor_Sample_Barcode){
    i.het<-inferHeterogeneity(maf = laml, tsb = i, vafCol = 't_vaf')
    write.table(unique(i.het$clusterData[,c(5,8)]),file=mathfile,
                append = TRUE,quote = FALSE,sep="\t",row.names = FALSE,col.names = FALSE)
    plotClusters(clusters = i.het,savePlot = TRUE)
  }
}
#lollipopPlot(maf = maf, gene = 'TP53', AACol = 'AAChange.refGene', showMutationRate = TRUE)
#lollipopPlot(maf = maf, gene = 'APC', AACol = 'AAChange.refGene', showMutationRate = TRUE)
#lollipopPlot(maf = maf, gene = 'KRAS', AACol = 'AAChange.refGene', showMutationRate = TRUE)
