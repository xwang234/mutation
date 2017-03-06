#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library("ComplexHeatmap")
readmutation=function(data,dataall)
{
  #uniq_genes=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/MutSigCV_1.4/uniqgenes_full192.txt")
  #uniq_genes=uniq_genes[,1]
  uniq_genes=unique(data$Hugo_Symbol)
  samples=unique(dataall$Tumor_Sample_Barcode)
  res=data.frame(matrix(0,nrow=length(uniq_genes),ncol=length(samples)+1))
  colnames(res)=c(samples,"count")
  rownames(res)=uniq_genes
  for (i in 1:length(uniq_genes))
  {
    gene=uniq_genes[i]
    tmptable=data[data$Hugo_Symbol==gene,]
    uniq_samples=unique(tmptable$Tumor_Sample_Barcode)
    for (j in 1:length(uniq_samples))
    {
      idx=which(colnames(res)==uniq_samples[j])
      res[i,idx]=sum(tmptable$Tumor_Sample_Barcode==uniq_samples[j])
    }
    res$count[i]=sum(res[i,1:length(samples)]>0)
  }
  res1=res[order(res$count,decreasing=T),]
}

readmutation1=function(data,dataall)
{
  #uniq_genes=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/MutSigCV_1.4/uniqgenes_full192.txt")
  #uniq_genes=uniq_genes[,1]
  uniq_genes=unique(data$Hugo_Symbol)
  samples=unique(dataall$Tumor_Sample_Barcode)
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
  
  res1=res[order(res$count,decreasing=T),]
  res1=res1[,1:length(samples)]
}
listtypes=function(mattable)
{
  types=NULL
  for (i in 1:nrow(mattable))
  {
    for (j in 1:ncol(mattable))
    {
      tmp=mattable[i,j]
      if (tmp != " ")
      {
        tmp1=unlist(strsplit(tmp,";",fixed=T))
        types=c(types,tmp1)
        types=unique(types)
      }
    }
  }
  print(types)
  return(types)
}

#for henan data:
wgstumors=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")
#varscandir="/fh/scratch/delete30/dai_j/henan/varscan"
#segfiles=paste0(varscandir,"/",wgstumors,".somatic.snp.Somatic.hc.annotated.reduced")
varscandir="/fh/scratch/delete30/dai_j/henan/varscan1"
segfiles=paste0(varscandir,"/",wgstumors,".snp.Somatic.hc.annotated.reduced")
name="EA CHINA"

#for dulak data:
wgstumors=c("SRR1001842","SRR1002713","SRR999423","SRR1001466","SRR1002670","SRR1001823","SRR999489","SRR1002343","SRR1002722","SRR1002656",
            "SRR1002929","SRR999438","SRR1001915","SRR999594","SRR1001868")
varscandir="/fh/scratch/delete30/dai_j/varscan"
segfiles=paste0(varscandir,"/",wgstumors,".snp.Somatic.hc.annotated.reduced")
varscandir="/fh/scratch/delete30/dai_j/varscan1"
segfiles=paste0(varscandir,"/",wgstumors,".snp.Somatic.hc.annotated.reduced")
name="EA US"

#for escc data:
wgstumors=paste0("T",c(1:6,8:18))
varscandir="/fh/scratch/delete30/dai_j/escc/varscan1"
segfiles=paste0(varscandir,"/",wgstumors,".snp.Somatic.hc.annotated.reduced")
name="ESCC CHINA"

nummutations=rep(0,length(wgstumors))
dataall=NULL
for (i in 1:length(segfiles))
{
  tmp=read.table(file=segfiles[i],header=T,sep="\t",stringsAsFactors=F,quote="")
  tmp$normal_var_freq=as.numeric(gsub("%","",tmp$normal_var_freq))/100
  num_normalvar=(tmp$normal_reads1+tmp$normal_reads2)*tmp$normal_var_freq
  tmp=tmp[num_normalvar<=2,]
  nummutations[i]=nrow(tmp)
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
numsnv=rep(0,length(wgstumors))
for (i in 1:length(wgstumors))
{
  tmptable=data[data$Tumor_Sample_Barcode==wgstumors[i],]
  numsnv[i]=nrow(tmptable)
}

mat=readmutation1(data,dataall)
num_mutations=sapply(1:nrow(mat),function(i){
  sum(mat[i,] != " ")
})
#get 20% cutoff
num_row=sum(num_mutations/ncol(mat)>=0.2)

mutation_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Missense_Mutation=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(1, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  Nonsense_Mutation=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(1, "mm"), gp = gpar(fill = "red", col = NA))
  },
  Splice_Site=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#008000", col = NA))
  },
  De_novo_Start_OutOfFrame = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "#F88000", col = NA))
  }
)

col = c("Missense_Mutation" = "blue", "Nonsense_Mutation"="red",
        "Splice_Site"="#008000", "De_novo_Start_OutOfFrame" = "#F88000")

ha = HeatmapAnnotation(SNV= anno_barplot(numsnv,axis=T,axis_gp = gpar(fontsize = 12),axis_side="right",border=F),
                       show_annotation_name = T,annotation_name_offset = unit(2, "cm"),gap = unit(3, "mm"))
oncoPrint(as.matrix(mat[1:num_row,]), get_type = function(x) strsplit(x, ";")[[1]],
          row_order = NULL,column_order = NULL,
          remove_empty_columns = FALSE,
          alter_fun = mutation_fun, col = col, 
          column_title = name,
          bottom_annotation=ha,
          bottom_annotation_height=unit(3,"cm"),
          heatmap_legend_param = list(title = "Mutations", at = c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site","De_novo_Start_OutOfFrame"), 
                                      labels = c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site","De_novo_Start_OutOfFrame"))
          )

