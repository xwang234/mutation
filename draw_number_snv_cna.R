#!/usr/bin/env Rscript

formdata=function(segfiles)
{
  dataall=NULL
  for (i in 1:length(segfiles))
  {
    tmp=read.table(file=segfiles[i],header=T,sep="\t",stringsAsFactors=F,quote="")
    tmp$normal_var_freq=as.numeric(gsub("%","",tmp$normal_var_freq))/100
    num_normalvar=(tmp$normal_reads1+tmp$normal_reads2)*tmp$normal_var_freq
    tmp=tmp[num_normalvar<=2,] #filter of normal_alt
    tmp$Tumor_Sample_Barcode=rep(wgstumors[i],nrow(tmp))
    #dataall=rbind.data.frame(dataall,tmp)
    #for varscan
    tmp1=cbind.data.frame(tmp[,1:11],tmp[,c("somatic_p_value","normal_reads1","normal_reads2","normal_var_freq","tumor_reads1","tumor_reads2","tumor_var_freq")])
    dataall=rbind.data.frame(dataall,tmp1)
  }
  # includedtypes1=c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
  #                 "Silent","Splice_Site","Translation_Start_Site","Nonstop_Mutation","RNA","Targeted_Region","De_novo_Start_InFrame","De_novo_Start_OutOfFrame")
  includedtypes2=c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                   "Splice_Site","Translation_Start_Site","Nonstop_Mutation","Targeted_Region","De_novo_Start_InFrame","De_novo_Start_OutOfFrame")
  data=dataall[dataall$Variant_Classification %in% includedtypes2,]
}
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

wgstumors=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")
varscandir="/fh/scratch/delete30/dai_j/henan/varscan2"
cnvnumfile="/fh/scratch/delete30/dai_j/henan/freec/henancnvnumcount.txt"


segfiles1=paste0(varscandir,"/",wgstumors,".snp.Somatic.hc.oncotator.reduced")
segfiles2=paste0(varscandir,"/",wgstumors,".indel.Somatic.hc.oncotator.reduced")


data1=formdata(segfiles1)
data2=formdata(segfiles2)
data=rbind.data.frame(data1,data2)

numsnv=rep(0,length(wgstumors))
for (i in 1:length(wgstumors))
{
  tmptable=data[data$Tumor_Sample_Barcode==wgstumors[i],]
  numsnv[i]=nrow(tmptable)
}
cnvnum=countcnvlength(cnvnumfile,wgstumors)
numcna=rowSums(cnvnum[,2:4]*10^6)

#form the matrix
mat=data.frame(matrix(NA,nrow=2*length(wgstumors),ncol=3))
colnames(mat)=c("sample","type","count")
mat$type=as.character(mat$type)
mat$sample=rep(wgstumors,2)
mat$type[1:length(wgstumors)]="snv"
mat$type[(1+length(wgstumors)):(2*length(wgstumors))]="cna"
mat$count[1:length(wgstumors)]=numsnv
mat$count[(1+length(wgstumors)):(2*length(wgstumors))]=numcna
mat$sample=factor(mat$sample,c("3A","11A","13A","15A","17A","25A","29A","33A","37A","41A"))


library(ggplot2)
ggplot(mat, aes(x=sample)) + geom_bar(data=subset(mat,type =="cna"), aes(y = count,fill=type), stat = "identity")
last_plot()+geom_bar(data=subset(mat,type == "snv"),aes(y = -count,fill=type), stat = "identity")+scale_y_continuous(breaks=c(-500,0,1000,2000,3000),labels=c(500,0,1000,2000,3000))
last_plot() + theme(plot.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),legend.title = element_blank(),axis.text=element_text(size=16),legend.text=element_text(size=16))
last_plot() + scale_fill_manual(labels=c("Copy number alterations", "SNV/indels"),values=c("skyblue", "blue"))
#last_plot() + theme(panel.background = element_blank(),axis.line = element_line(size=3, colour = "black"))
last_plot() + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

