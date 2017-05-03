#!/usr/bin/env Rscript
library("beeswarm")
library("ComplexHeatmap")
library(circlize)

grp1tumors=c('SRR1001842','SRR1002713','SRR999423','SRR1001466','SRR1002670','SRR1001823','SRR999489','SRR1002343','SRR1002722','SRR1002656','SRR1002929','SRR999438','SRR1001915','SRR999594','SRR1001868','SRR1001635')
grp2tumors=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")

grp1dir="/fh/scratch/delete30/dai_j/mutect4"
grp2dir="/fh/scratch/delete30/dai_j/henan/mutect5"

load("compare_globals_2grps_mutect4.RData")
load("mutationrate.RData")
computep=function(x,y)
{
  
  print(paste0("meanx:",mean(x,na.rm = T)))
  print(paste0("meany:",mean(y,na.rm = T)))
  print(paste0("medianx:",median(x,na.rm = T)))
  print(paste0("mediany:",median(y,na.rm = T)))
  print(paste0("rangex:",range(x)))
  print(paste0("rangey:",range(y)))
  p=2*(1-pnorm(abs((mean(x)-mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))))
  res=list(min1=min(x,na.rm=T),max1=max(x,na.rm=T),min2=min(y,na.rm=T),max2=max(y,na.rm=T),
           mean1=mean(x),mean2=mean(y),median1=median(x),median2=median(y),pvalue=p)
  print(p)
  return(res)
}
computep_wilcox=function(x,y)
{
  # x=unique(x)
  # y=unique(y)
  # z=intersect(x,y)
  # x=x[! x %in% z]
  # y=y[! y %in% z]
  print(paste0("meanx:",mean(x,na.rm = T)))
  print(paste0("meany:",mean(y,na.rm = T)))
  print(paste0("medianx:",median(x,na.rm = T)))
  print(paste0("mediany:",median(y,na.rm = T)))
  print(paste0("rangex:",range(x)))
  print(paste0("rangey:",range(y)))
  p=wilcox.test(x,y)$p.value
  res=list(min1=min(x,na.rm=T),max1=max(x,na.rm=T),min2=min(y,na.rm=T),max2=max(y,na.rm=T),
           mean1=mean(x),mean2=mean(y),median1=median(x),median2=median(y),pvalue=p)
}


restable=data.frame(matrix(NA,nrow=length(grp1tumors)+length(grp2tumors),ncol=2))
colnames(restable)=c("dataset","sample")
restable[1:length(grp1tumors),1]="US-EAC"
restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),1]="CH-EAC"
restable[,2]=c(grp1tumors,grp2tumors)
#reorder group
restable$dataset=factor(restable$dataset,c("US-EAC","CH-EAC"))

#get data
restable=cbind(restable,totalrate=c(grp1mutationrate$mutationrate$total,grp2mutationrate$mutationrate$total))
counttransrate=function(transratemat,tumors,opt="a2c")
{
  res=data.frame(matrix(NA,ncol=2,nrow=length(tumors)))
  colnames(res)=c("tumor","mutationrate")
  idx=which(rownames(transratemat)==opt)
  for (numpair in 1:length(tumors))
  {
    res[numpair,1]=tumors[numpair]
    res[numpair,2]=transratemat[idx,numpair]
  }
  return(res)
}

#mutation a>c
grp1a2c=counttransrate(dulak_res$transratemat,grp1tumors)
grp2a2c=counttransrate(henan_res$transratemat,grp2tumors)
restable=cbind(restable,c(grp1a2c[,2],grp2a2c[,2]))
colnames(restable)[ncol(restable)]="a2c"
grp1a2g=counttransrate(dulak_res$transratemat,grp1tumors,opt="a2g")
grp2a2g=counttransrate(henan_res$transratemat,grp2tumors,opt="a2g")
restable=cbind(restable,c(grp1a2g[,2],grp2a2g[,2]))
colnames(restable)[ncol(restable)]="a2g"
#a2t
grp1a2t=counttransrate(dulak_res$transratemat,grp1tumors,opt="a2t")
grp2a2t=counttransrate(henan_res$transratemat,grp2tumors,opt="a2t")
restable=cbind(restable,c(grp1a2t[,2],grp2a2t[,2]))
colnames(restable)[ncol(restable)]="a2t"
grp1c2a=counttransrate(dulak_res$transratemat,grp1tumors,opt="c2a")
grp2c2a=counttransrate(henan_res$transratemat,grp2tumors,opt="c2a")
restable=cbind(restable,c(grp1c2a[,2],grp2c2a[,2]))
colnames(restable)[ncol(restable)]="c2a"
#c2g
grp1c2g=counttransrate(dulak_res$transratemat,grp1tumors,opt="c2g")
grp2c2g=counttransrate(henan_res$transratemat,grp2tumors,opt="c2g")
restable=cbind(restable,c(grp1c2g[,2],grp2c2g[,2]))
colnames(restable)[ncol(restable)]="c2g"
#c2t
grp1c2t=counttransrate(dulak_res$transratemat,grp1tumors,opt="c2t")
grp2c2t=counttransrate(henan_res$transratemat,grp2tumors,opt="c2t")
restable=cbind(restable,c(grp1c2t[,2],grp2c2t[,2]))
colnames(restable)[ncol(restable)]="c2t"

#copynumber-------------------
#length computed by compute_cnvlength.R
countcnvlength=function(cnvlenfile,tumors)
{
  res=data.frame(matrix(NA,ncol=4,nrow=length(tumors)))
  colnames(res)=c("tumor","lengain","lenloss","lenloh")
  res[,1]=tumors
  tmp=read.table(cnvlenfile,header=T)
  res[,2]=rowSums(tmp[,2:25])/10^6
  res[,3]=rowSums(tmp[,26:49])/10^6
  res[,4]=rowSums(tmp[,50:73])/10^6
  return(res)
}
cnvlenfile1="/fh/scratch/delete30/dai_j/freec/dulakcnvlengthcount.txt"
grp1cnvlen=countcnvlength(cnvlenfile1,grp1tumors)
cnvlenfile2="/fh/scratch/delete30/dai_j/henan/freec/henancnvlengthcount.txt"
grp2cnvlen=countcnvlength(cnvlenfile2,grp2tumors)
restable=cbind(restable,c(grp1cnvlen[,2],grp2cnvlen[,2]))
colnames(restable)[ncol(restable)]="lengain"
restable=cbind(restable,c(grp1cnvlen[,3],grp2cnvlen[,3]))
colnames(restable)[ncol(restable)]="lenloss"
restable=cbind(restable,c(grp1cnvlen[,4],grp2cnvlen[,4]))
colnames(restable)[ncol(restable)]="lenloh"

#beeswarm plot
addpvlue=function(data1,data2,ymax,x1,x2)
{
  pvalue=computep_wilcox(data1,data2)$pvalue
  pvalue=format(pvalue,digits=2)
  ymax1=max(data1,na.rm=T)
  ymax2=max(data2,na.rm=T)
  ymax12=max(ymax1,ymax2)
  segments(x1,ymax12*1.05+ymax*0.05,x2,ymax12*1.05+ymax*0.05,lwd=1.2)
  segments(x1,ymax12*1.05+ymax*0.05,x1,ymax12*1.05,lwd=1.2)
  segments(x2,ymax12*1.05+ymax*0.05,x2,ymax12*1.05,lwd=1.2)
  text(mean(c(x1,x2)),ymax12*1.05+ymax*0.1,paste0("p=",pvalue),cex=1)
}

plotbeeswarm=function(mat,formula,colors=c("skyblue", "forestgreen","black", "red", "gray", "tan3"),xlabel,ylabel,labels)
{
  ymax=max(mat$value,na.rm=T)
  numgroups=length(unique(mat[,2]))
  colors1=rep(colors,each=2)
  ats_c=rep(0,numgroups)
  for (i in 1:numgroups) ats_c[i]=4*i
  
  ats=rep(0,2*numgroups)
  for (i in 1:numgroups)
  {
    ats[(i-1)*2+1]=ats_c[i]-1
    ats[(i-1)*2+2]=ats_c[i]+1
  }
  if (numgroups==1)
  {
    colors1=colors
    beeswarm(as.formula(formula), data = mat, labels=c("US-EAC","CH-EAC"),
             col = colors1, pch = 19,add = FALSE, cex=1.2,ylim=c(0,ymax*1.25),main="",
             xlab=xlabel,ylab=ylabel,cex.axis=1.2,cex.main=1.2,cex.lab=1.2,at=ats)
  }else
  {
    beeswarm(as.formula(formula), data = mat, xaxt='n',
             col = colors1, pch = 19,add = FALSE, cex=1.2,ylim=c(0,ymax*1.25),main="",
             xlab=xlabel,ylab=ylabel,cex.axis=1.2,cex.main=1.2,cex.lab=1.2,at=ats,corral = "wrap")
    axis(side=1, at=ats_c, labels=labels,
         cex.axis=1.2,cex.lab=1.2,line=F,lwd.ticks=1.2)
  }
  #add median
  for (i in 1:numgroups)
  {
    groupdata1=mat[mat[,1]==unique(mat[,1])[1] & mat[,2]==unique(mat[,2])[i] ,3]
    groupdata2=mat[mat[,1]==unique(mat[,1])[2] & mat[,2]==unique(mat[,2])[i] ,3]
    median1=median(groupdata1,na.rm=T)
    median2=median(groupdata2,na.rm=T)
    if (numgroups==1)
    {
      segments(ats_c[i]-1-1.1,median1,ats_c[i]-1+1.1,median1,col=colors1[(i-1)*2+1],lwd=1.3,lty=1)
      segments(ats_c[i]+1-1.1,median2,ats_c[i]+1+1.1,median2,col=colors1[(i-1)*2+2],lwd=1.3,lty=1)
    }else
    {
      segments(ats_c[i]-1-1.6,median1,ats_c[i]-1+1.6,median1,col=colors1[(i-1)*2+1],lwd=1.3,lty=1)
      segments(ats_c[i]+1-1.6,median2,ats_c[i]+1+1.6,median2,col=colors1[(i-1)*2+2],lwd=1.3,lty=1)
    }
    addpvlue(groupdata1,groupdata2,ymax,ats_c[i]-1,ats_c[i]+1)
  }
}

formmatrix=function(selidx,types)
{
  numgroups=length(selidx)
  mat=data.frame(matrix(NA,nrow=2*numgroups,ncol=length(grp1tumors)))
  for (i in 1:numgroups)
  {
    groupdata1=restable[restable$dataset==unique(restable$dataset)[1],selidx[i]]
    groupdata2=restable[restable$dataset==unique(restable$dataset)[2],selidx[i]]
    mat[(i-1)*2+1,]=groupdata1
    mat[(i-1)*2+2,1:length(groupdata2)]=groupdata2
  }
  alllocations=rep(c("US","CH"),numgroups)
  alltypes=rep(types,each=2)
  mat$locations=alllocations
  mat$types=alltypes
  mat$locations=factor(mat$locations,c("US","CH"))
  mat$types=factor(mat$types,types)
  library(reshape2)
  mat=melt(mat, id = c('locations', 'types'))
  mat=mat[,-3]
  mat$value=as.numeric(mat$value)
  return(mat)
}

#plot fig1a
plotfig1a=function()
{
  postscript(file="fig1a.ps",width = 6,height = 6,paper= "special",horizontal=F)
  par(mar=c(4.1,4.5,4.1,2.1))
  idx=which(colnames(restable)=="totalrate")
  mat=formmatrix(idx,types="totalrate")
  plotbeeswarm(mat,formula="value ~ locations + types",xlabel="",ylabel="Mutations/Mb")
  dev.off()
}
plotfig1a()

#plot fig1b
plotfig1b=function()
{
  postscript(file="fig1b.ps",width = 8,height = 6,paper= "special",horizontal=F)
  par(mar=c(4.1,4.5,4.1,2.1))

  idx=match(c("c2a","c2g","c2t","a2g","a2c","a2t"),colnames(restable))
  mat=formmatrix(idx,types=c("c2a","c2g","c2t","a2g","a2c","a2t"))
  plotbeeswarm(mat,formula="value ~ locations + types",xlabel="Type of mutation",ylabel="Fraction of mutation",labels=c("C>A","C>G","C>T","A>G","A>C","A>T"))
  dev.off()
}

#plot fig1c
plotfig1c=function()
{
  # library(png)
  # img1c = readPNG("figure1c.png")
  # plot.new()
  # usr<-par("usr")
  # rasterImage(img1c, usr[1], usr[3], usr[2], usr[4],interpolate = F)
  #fig1c.png was created in ppt
  library("imager")
  img1c=load.image("figure1c.png")
  par(mar=c(1.1,1.5,1.1,1.1))
  plot(img1c,axes=F)
}

#figure 2
plotsignature3=function (signature, name, y_limit=NULL) 
{
  #sort signature subtypes
  mutationsubtypes=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/mutationsignature/mutationsubtypes.txt",stringsAsFactors = F)
  trinucleotide=rep("",nrow(mutationsubtypes))
  for (i in 1:nrow(mutationsubtypes))
  {
    tmp=unlist(strsplit(mutationsubtypes[i,1],"[",fixed = T))
    tmp1=unlist(strsplit(tmp[2],"]"))
    trinucleotide[i]=paste0(tmp[1],substr(tmp[2],1,1),tmp1[2])
  }
  #for lego data
  sortsinature=function(signature,mutationsubtypes)
  {
    signature_colnames=colnames(signature)
    signature_names=rep(NA,length(signature_colnames))
    for (i in 1:length(signature_colnames))
    {
      tmp=unlist(strsplit(signature_colnames[i],"_"))
      tmp1=unlist(strsplit(tmp[2],"x"))
      tmp2=NA
      if (tmp[1]=="C.T.G.A") tmp2="C>T"
      if (tmp[1]=="C.A.G.T") tmp2="C>A"
      if (tmp[1]=="C.G.G.C") tmp2="C>G"
      if (tmp[1]=="T.C.A.G") tmp2="T>C"
      if (tmp[1]=="T.G.A.C") tmp2="T>G"
      if (tmp[1]=="T.A.A.T") tmp2="T>A"
      signature_names[i]=paste0(tmp1[1],"[",tmp2,"]",tmp1[2])
    }
    idx=match(mutationsubtypes[,1],signature_names)
    res=signature[,idx]
  }
  if (class(signature)=="data.frame")
    if (grepl(".",colnames(signature)[1],fixed=T) & grepl("x",colnames(signature)[1],fixed=T))
    {
      signature=sortsinature(signature,mutationsubtypes)
    }
  if (max(signature)>1)
  {
    signature=signature/sum(signature)
  }
  
  if (is.null(y_limit))
    y_limit <- 1.4 * max(signature)
  
  barplot(signature, names.arg = "",
          cex.names = 0.7, las = 2, col = NA, ylim = c(-0.02, y_limit),
          border = NA, xaxt = "n", ann = FALSE, yaxt = "n", space = 0.3)
  # #box()
  mutationtypes=c('C>A','C>G','C>T','A>T','A>G','A>C')
  mutationtypecols=c("skyblue", "black", "red", "gray", "forestgreen", "tan3")
  
  cols=rep(mutationtypecols,each=16)
  
  x=barplot(signature, names.arg ="", 
            cex.names = 0.7, cex.lab=1.3,las = 2, col = cols, 
            ylim = c(0, y_limit), border = NA, space = 0.3, main = name, 
            ylab = "Mutation type fraction", add = T,cex.axis=1.5, axes=F, font.axis = 2)
  #write labels
  for (i in 1:length(x))
  {
    text(x[i],-0.002,substr(trinucleotide[i],1,1),col="black",srt=90,cex=0.5,font=2)
    text(x[i],-0.003,substr(trinucleotide[i],2,2),col=cols[i],srt=90,cex=0.5,font=2)
    text(x[i],-0.004,substr(trinucleotide[i],3,3),col="black",srt=90,cex=0.5,font=2)
  }
  axis(side=2, at=seq(0,0.05,0.01),cex.axis=1.2,cex.lab=1.2,line=F,lwd.ticks=1.1)
  for (i in seq(from = 0, to = y_limit-0.01, by = 0.01))
  {
    segments(par("usr")[1],i,125,i,col="gray66")
  }
  # abline(h = seq(from = 0, to = y_limit-0.01, by = 0.01), col="gray66")
  y1=max(signature)+0.005+y_limit*0.03
  # y2=y1+(y_limit-max(signature))/4
  y2=y1+0.002
  print((y1-max(signature)-0.005)/y_limit)
  for (i in 1:6)
  {
    rect(x[(i-1)*16+1],max(signature)+0.005,x[(i-1)*16+15],y1,col=mutationtypecols[i],border=F)
    text(x[(i-1)*16+8],y2,labels=mutationtypes[i],col=mutationtypecols[i],cex=1.5)
  }
  
  return(signature)
}

#stability=.8921    0.8090    0.7262    0.6603 (processStabAvg)
signaturemat=NULL
for (i in 1:length(grp2tumors))
{
  tmp=plotsignature2(signature=henan_res$legomat[i,], name=grp2tumors[i]) 
  signaturemat=rbind.data.frame(signaturemat,tmp)
}
save(signaturemat,file="henan_signaturedata.RData")

signaturemat1=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/mutationsignature/henan_mutect5_2_signatures.txt")

#plot figure2
plotfig2a=function()
{
  postscript(file="fig2a.ps",width = 8,height = 6,paper= "special",horizontal=F)
  par(mar=c(0.1,5.1,1.1,2.1))
  tmp=plotsignature3(signature = signaturemat1[,1],name="")
  dev.off()
}
plotfig2b=function()
{
  postscript(file="fig2b.ps",width = 8,height = 6,paper= "special",horizontal=F)
  par(mar=c(0.1,5.1,1.1,2.1))
  tmp=plotsignature3(signature = signaturemat1[,2],name="")
  dev.off()
}
plotfig2c=function()
{
  postscript(file="fig2c.ps",width = 8,height = 6,paper= "special",horizontal=F)
  par(mar=c(0.1,5.1,1.1,2.1))
  tmp=plotsignature3(signature = signaturemat1[,3],name="")
  dev.off()
}
outfig="fig2a.ps"
postscript(outfig, width = 8, height = 8)

par(mfrow=c(2,1))
par(mar=c(1.1,5.1,1.1,2.1))
mtexts=c("(a)","(b)")
for (i in 1:2)
{
  tmp=plotsignature3(signature = signaturemat1[,i],name="")
  text(60,-0.015,mtexts[i],cex=1.5)
}
dev.off()

#figure3
readmutation1=function(data,wgstumors)
{
  #uniq_genes=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/MutSigCV_1.4/uniqgenes_full192.txt")
  #uniq_genes=uniq_genes[,1]
  uniq_genes=unique(data$Hugo_Symbol)
  samples=wgstumors
  res=data.frame(matrix(" ",nrow=length(uniq_genes),ncol=length(samples)))
  colnames(res)=samples
  rownames(res)=uniq_genes
  for (i in 1:ncol(res))
  {
    res[,i]=as.character(res[,i])
  }
  for (i in 1:length(uniq_genes))
  {
    gene=uniq_genes[i]
    tmptable=data[data$Hugo_Symbol==gene,]
    uniq_samples=unique(tmptable$Tumor_Sample_Barcode)
    for (j in 1:length(uniq_samples))
    {
      idx=which(tmptable$Tumor_Sample_Barcode==uniq_samples[j])
      idx1=which(samples==uniq_samples[j])
      classes=c()
      for (k in 1:length(idx))
      {
        classes[k]=tmptable$Variant_Classification[idx[k]]
      }
      classes=unique(classes)
      classes=paste0(classes,collapse = ";")
      classes=paste0(classes,";")
      res[i,idx1]=classes
    }
    res$count[i]=sum(res[i,1:length(samples)] != " ")
  }
  #process indels
  for (i in 1:ncol(res))
  {
    res[,i]=gsub("In_Frame_Ins","Inframe_Indel",res[,i])
    res[,i]=gsub("Frame_Shift_Ins","Inframe_Indel",res[,i])
    res[,i]=gsub("Frame_Shift_Del","Inframe_Indel",res[,i])
  }
  res1=res[order(as.integer(res$count),decreasing=T),]
  res1=res1[,1:length(samples)]
}
formdata=function(segfiles)
{
  dataall=NULL
  for (i in 1:length(segfiles))
  {
    tmp=read.table(file=segfiles[i],header=T,sep="\t",stringsAsFactors=F,quote="")
    dataall=rbind.data.frame(dataall,tmp)
  }
  # includedtypes1=c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
  #                 "Silent","Splice_Site","Translation_Start_Site","Nonstop_Mutation","RNA","Targeted_Region","De_novo_Start_InFrame","De_novo_Start_OutOfFrame")
  includedtypes2=c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                   "Splice_Site","Translation_Start_Site","Nonstop_Mutation","Targeted_Region","De_novo_Start_InFrame","De_novo_Start_OutOfFrame")
  data=dataall[dataall$Variant_Classification %in% includedtypes2,]
  idx=which(data$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins","De_novo_Start_OutOfFrame"))
  data$Variant_Classification[idx]="Frame_shift"
  idx=which(data$Variant_Classification %in% c("In_Frame_Del","In_Frame_Ins","De_novo_Start_InFrame"))
  data$Variant_Classification[idx]="In_frame"
  idx=which(data$Variant_Classification %in% "Missense_Mutation")
  data$Variant_Classification[idx]="Missense"
  idx=which(data$Variant_Classification %in% "Nonsense_Mutation")
  data$Variant_Classification[idx]="Nonsense"
  idx=which(data$Variant_Classification %in% "Splice_Site")
  data$Variant_Classification[idx]="Splice_site"
  return(data)
}
name="CH-EAC"
mutectdir="/fh/scratch/delete30/dai_j/henan/mutect5"
segfiles1=paste0(mutectdir,"/",grp2tumors,".mutectmaf.ancotator")
data=formdata(segfiles1)
nrow(data)
#[1] 743
numsnv=rep(0,length(grp2tumors))
for (i in 1:length(grp2tumors))
{
  tmptable=data[data$Tumor_Sample_Barcode==grp2tumors[i],]
  numsnv[i]=nrow(tmptable)
}

mat=readmutation1(data,grp2tumors)
num_mutations=sapply(1:nrow(mat),function(i){
  sum(mat[i,] != " ")
})
#get 20% cutoff
num_row=sum(num_mutations/ncol(mat)>=0.2)
rownames(mat)[1:num_row]
mart=useMart("ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org/biomart/martservice/", dataset="hsapiens_gene_ensembl") #GRCh37
getgenelength=function(gene)
{
  attributes <- c("chromosome_name","start_position","end_position")
  filters <- c("hgnc_symbol")
  values=list(hgnc_symbol=gene)
  res <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
  res1=NA
  if (nrow(res)>0)
  {
    res1=res[1,3]-res[1,2]
  }
  return(res1)
}

genelength=rep(NA,num_row)
for (i in 1:num_row)
{
  genelength[i]=getgenelength(rownames(mat)[i])
}

mutation_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Missense=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "blue", col = NA),just="top")
  },
  Nonsense=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "red", col = NA),just="bottom")
  },
  Splice_site=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#008000", col = NA))
  },
  Inframe=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "black", col = NA),just="centre")
  },
  Frame_shift=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "#F88000", col = NA),just="centre")
  }
)

col = c("Missense" = "blue", "Nonsense"="red", "Splice_site"="#008000", "In_frame"="black", "Frame_shift"="#F88000")
#Figure3
plotfig3=function()
{
  png("fig3.png", width = 6, height = 8, units = 'in', res=300)
  # ha = HeatmapAnnotation(SNV= anno_barplot(numsnv,axis=T,axis_gp = gpar(fontsize = 12),axis_side="right",border=F,ylim=c(0,500)),
  #                        show_annotation_name = T,annotation_name_offset = unit(2, "cm"),gap = unit(3, "mm"))
  oncoPrint(as.matrix(mat[1:num_row,]), get_type = function(x) strsplit(x, ";")[[1]],
            row_order = NULL,column_order = NULL,
            remove_empty_columns = FALSE,
            alter_fun = mutation_fun, col = col, 
            pct_gp=gpar(fontsize=12),
            axis_gp = gpar(fontsize = 12),
            #column_title = name
            column_title = ""
            # bottom_annotation=ha,
            # bottom_annotation_height=unit(3,"cm"),
            # heatmap_legend_param = list(title = "Mutations", at = c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site","Inframe_Indel","De_novo_Start_OutOfFrame"), 
            #                             labels = c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site","Inframe_Indel","De_novo_Start_OutOfFrame"))
  )
  dev.off()
  
}

#figure4
#plot fig4(a)
plotfig4a=function()
{
  postscript(file="fig4a.ps",width = 8,height = 6,paper= "special",horizontal=F)
  par(mar=c(4.1,4.5,4.1,2.1))
  
  idx=match(c("lengain","lenloss","lenloh"),colnames(restable))
  mat=formmatrix(idx,types=c("lengain","lenloss","lenloh"))
  plotbeeswarm(mat,formula="value ~ locations + types",xlabel="Type of CNA",ylabel="Length of CNA/Mb",labels=c("Gain","Loss","LOH"))
  dev.off()
}
plotfig4b=function()
{
  library("imager")
  if (!file.exists("fig4b.png"))
  {
    gisticdir2="/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel1_brlen0.98_broad1"
    img4b1=load.image(paste0(gisticdir2,"/del_qplot.png"))
    img4b2=load.image(paste0(gisticdir2,"/amp_qplot.png"))
    png(filename = "fig4b.png",res=300,width = 8,height = 6,units ="in")
    par(mar=c(0.3,0.5,0.2,0.3))
    par(mfrow=c(1,2))
    plot(img4b1,axes=F)
    plot(img4b2,axes=F)
    dev.off()
  }
}


getwholetable1=function(alllessionsfile,opt="AMP",opt1="peak",opt2="binary",qcutoff=1,delcutoff=-0.9)
{
  
  alllessiontable=read.table(alllessionsfile,header=T,sep="\t",stringsAsFactors=F)
  #remove the last NA column introduced by \t
  idx=which(colSums(is.na(alllessiontable))==nrow(alllessiontable))
  alllessiontable=alllessiontable[,-idx]
  idx=which(colnames(alllessiontable)=="Amplitude.Threshold")
  idx=idx+1 # tumor samples start
  idxkeep=which(rowSums(abs(alllessiontable[1:(nrow(alllessiontable)/2),idx:ncol(alllessiontable)]),na.rm = T)>1)
  alllessionrealdata=alllessiontable[(nrow(alllessiontable)/2+1):nrow(alllessiontable),]
  idxdel=which(grepl("Amplification",alllessionrealdata$Unique.Name)==F)
  for (i in idxdel)
  {
    for (j in idx:ncol(alllessionrealdata))
    {
      if (alllessionrealdata[i,j]<delcutoff)
      {
        alllessionrealdata[i,j]=alllessiontable[i,j]=2
      }
      
    }
  }
  
  if (opt2=="binary")
  {
    alllessiontable=alllessiontable[1:(nrow(alllessiontable)/2),]
  }else
  {
    alllessiontable=alllessiontable[(nrow(alllessiontable)/2+1):nrow(alllessiontable),]
  }
  alllessiontable=alllessiontable[idxkeep,]
  alllessiontable$Descriptor=gsub(" ","",alllessiontable$Descriptor,fixed=T)
  idx1=grepl("Amplification",alllessiontable$Unique.Name)
  if (opt=="AMP")
  {
    alllessiontable=alllessiontable[idx1,]
  }else
  {
    alllessiontable=alllessiontable[!idx1,]
  }
  alllessiontable=alllessiontable[alllessiontable$Residual.q.values.after.removing.segments.shared.with.higher.peaks<qcutoff,]
  res=data.frame(matrix(NA,nrow=nrow(alllessiontable),ncol=5))
  colnames(res)=c("cytoband","chr","start","end","qvalue")
  res$chr=as.character(res$chr)
  res$start=as.integer(res$start)
  res$end=as.integer(res$end)
  res$cytoband=alllessiontable$Descriptor
  res$qvalue=alllessiontable$Residual.q.values.after.removing.segments.shared.with.higher.peaks
  res=cbind.data.frame(res,alllessiontable[,idx:ncol(alllessiontable)])
  for (i in 1:nrow(res))
  {
    #tmp=unlist(strsplit(alllessiontable$Wide.Peak.Limits[i],":",fixed=T))
    if (opt1=="peak")
    {
      tmp=unlist(strsplit(alllessiontable$Peak.Limits[i],":",fixed=T))
    }else
    {
      tmp=unlist(strsplit(alllessiontable$Region.Limits[i],":",fixed=T))
    }
    
    res$chr[i]=tmp[1]
    tmp=unlist(strsplit(tmp[2],"(",fixed=T))
    tmp=unlist(strsplit(tmp[1],"-",fixed=T))
    res$start[i]=as.integer(tmp[1])
    res$end[i]=as.integer(tmp[2])
  }
  gr_res=GRanges(seqnames = res$chr,ranges=IRanges(start=res$start,end=res$end))
  gr_res1=reduce(gr_res)
  if (length(gr_res)>length(gr_res1)) #overlap regions
  {
    res1=NULL
    for (i in 1:length(gr_res1))
    {
      olapidx=NULL
      for (j in 1:length(gr_res))
      {
        olap=intersect(gr_res[j],gr_res1[i])
        if (length(olap)>0)
          olapidx=c(olapidx,j)
      }
      if (length(olapidx)==1)
      {
        res1=rbind.data.frame(res1,res[olapidx,])
      }
      if (length(olapidx)>1)
      {
        widths=res$end[olapidx]-res$start[olapidx]
        if (max(widths)==min(widths))
        {
          qvalues=res$qvalue[olapidx]
          selidx=olapidx[which.min(qvalues)]
        }else
        {
          selidx=olapidx[which.max(widths)]
        }
        res1=rbind.data.frame(res1,res[selidx,])
      }
    }
    res=res1
  }
  #add rownames, deal with duplicated cytobands
  idxdup=duplicated(res$cytoband)
  rownames(res)[!idxdup]=res$cytoband[!idxdup]
  dupsign=0
  if (sum(idxdup)==1)
  {
    rownames(res)[idxdup]=paste0(res$cytoband[idxdup],"(1)")
  }
  if (sum(idxdup)>1)
  {
    idxdup=which(idxdup==T)
    rownames(res)[idxdup[1]]=paste0(res$cytoband[idxdup[1]],"(1)")
    dupsign=1
    for (i in 2:length(idxdup))
    {
      if(res$cytoband[idxdup[i]]==res$cytoband[idxdup[i-1]])
      {
        dupsign=dupsign+1
        rownames(res)[idxdup[i]]=paste0(res$cytoband[idxdup[i]],"(",dupsign,")")
      }else
      {
        dupsign=1
        rownames(res)[idxdup[i]]=paste0(res$cytoband[idxdup[i]],"(1)")
      }
    }
  }
  
  return(res)
}
gisticdir="/fh/scratch/delete30/dai_j/henan/gistic2data/dulak_henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
alllessionsfile=paste0(gisticdir,"/all_lesions.conf_95.txt")

ampsegdata=getwholetable1(alllessionsfile,opt="AMP",opt2="continuous",qcutoff = 0.01)
delsegdata=getwholetable1(alllessionsfile,opt="DEL",opt2="continuous",qcutoff = 0.01)
colnames(ampsegdata)[6:21]=paste0("US-EAC",1:16)
colnames(ampsegdata)[22:ncol(ampsegdata)]=paste0("CH-EAC",1:10)
colnames(delsegdata)=colnames(ampsegdata)
#continuous copyratio -> discreet
callevents=function(segdata)
{
  res=segdata
  for (i in 1:nrow(segdata))
  {
    for (j in 6:ncol(segdata))
    {
      res[i,j]=0
      if (segdata[i,j] !=0 & segdata[i,j]< -0.9)
      {
        res[i,j]=-2
      }else if (segdata[i,j] !=0 & segdata[i,j]< -0.1)
      {
        res[i,j]=-1
      }else if (segdata[i,j] !=0 & segdata[i,j]>0.9)
      {
        res[i,j]=2
      }else if (segdata[i,j] !=0 & segdata[i,j]>0.1)
      {
        res[i,j]=1
      }
    }
  }
  return(res)
}
ampsegdata=callevents(ampsegdata)
delsegdata=callevents(delsegdata)
#only keep 2, set 1 as 0. >2^0.9 or <2^(-0.9)
useampordel1=function(segdata,type="amp")
{
  for (i in 1:nrow(segdata))
  {
    for (j in 6:ncol(segdata))
    {
      if (segdata[i,j]== -1 | segdata[i,j]== 1)
        segdata[i,j]=0
    }
  }
  #remove high frequency regions < 2/3 #samples
  num=apply(segdata[,6:ncol(segdata)],1,function(x){
    sum(x!=0)
  })
  segdata=segdata[num<2/3*(ncol(segdata)-6) & num>1,]
  #remove regions have conflict def
  numamp=apply(segdata[6:ncol(segdata)],1,function(x){
    sum(x>0)
  })
  numdel=apply(segdata[6:ncol(segdata)],1,function(x){
    sum(x<0)
  })
  if (type=="amp")
  {
    idx=which(numamp>=numdel)
  }else
  {
    idx=which(numdel>=numamp)
  }
  segdata=segdata[idx,]
  return(segdata)
}
ampsegdata_=useampordel1(ampsegdata)
delsegdata_=useampordel1(delsegdata,type="del")
allsegdata=rbind.data.frame(ampsegdata_,delsegdata_)

tmp=which(rownames(allsegdata)=="21p11.2(1)")
rownames(allsegdata)[tmp]="21p11.2"
tmp=which(rownames(allsegdata)=="21p11.2(2)")
rownames(allsegdata)[tmp]="21p11.2(1)"
tmp=which(rownames(allsegdata)=="19q13.331")
rownames(allsegdata)[tmp]="19q13.33(1)"
plotfig4c=function()
{
  png(filename = "fig4c.png",res=300,width = 6,height = 4,units ="in")
  Heatmap(allsegdata[,6:ncol(allsegdata)],cluster_rows = F,name="Event",show_column_names=T,show_column_dend=T,
          rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","gray88","red"),
          column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
          heatmap_legend_param = list(at = c(0,2,-2), labels = c("None", "Amplification","Deletion"),title_gp = gpar(fontsize = 12),labels_gp = gpar(fontsize = 12)),
          row_names_gp = gpar(fontsize = 12)
  )
  dev.off()
}

