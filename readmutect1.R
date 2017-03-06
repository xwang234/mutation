#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
output=as.character(args[1])
segfiles=c() #.maf_keep.annotated.reduced
for (i in 2:length(args))
{
  segfiles[i-1]=as.character(args[i])
}

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
wgstumors=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")
# mutsigdir="/fh/scratch/delete30/dai_j/henan/mutsig"
# mutsigdir="."
# segfiles=paste0(mutsigdir,"/",wgstumors,".maf_keep.annotated.reduced")
# output="/fh/scratch/delete30/dai_j/henan/mutect/henan_mutectcount.txt"
# output1="./henan_mutectcount_with_RNA.txt"
# output2="./henan_mutectcount_without_RNA.txt"
# output4="./henan_mutectmutation_without_RNA.txt"

varscandir="/fh/scratch/delete30/dai_j/henan/varscan"
segfiles=paste0(varscandir,"/",wgstumors,".somatic.snp.Somatic.annotated.reduced")


#combine the segfile

wgstumors=c("SRR1001842","SRR1002713","SRR999423","SRR1001466","SRR1002670","SRR1001823","SRR999489","SRR1002343","SRR1002722","SRR1002656",
"SRR1002929","SRR999438","SRR1001915","SRR999594","SRR1001868")
mutsigdir="/fh/scratch/delete30/dai_j/mutsig"
segfiles=paste0(mutsigdir,"/",wgstumors,".maf_keep.annotated.reduced")
output="/fh/scratch/delete30/dai_j/mutect1/dulak_mutectcount.txt"

varscandir="/fh/scratch/delete30/dai_j/varscan"
segfiles=paste0(varscandir,"/",wgstumors,".snp.Somatic.annotated.reduced")


dataall=NULL
for (i in 1:length(segfiles))
{
  tmp=read.table(file=segfiles[i],header=T,sep="\t",stringsAsFactors=F,quote="")
  tmp$Tumor_Sample_Barcode=rep(wgstumors[i],nrow(tmp))
  #dataall=rbind.data.frame(dataall,tmp)
  #for varscan
  tmp1=cbind.data.frame(tmp[,1:11],tmp[,c("normal_reads1","normal_reads2","normal_var_freq","tumor_reads1","tumor_reads2","tumor_var_freq")])
  dataall=rbind.data.frame(dataall,tmp1)
}
# includedtypes1=c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
#                 "Silent","Splice_Site","Translation_Start_Site","Nonstop_Mutation","RNA","Targeted_Region","De_novo_Start_InFrame","De_novo_Start_OutOfFrame")
includedtypes2=c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                "Splice_Site","Translation_Start_Site","Nonstop_Mutation","Targeted_Region","De_novo_Start_InFrame","De_novo_Start_OutOfFrame")

# data1=dataall[dataall$Variant_Classification %in% includedtypes1,]
data2=dataall[dataall$Variant_Classification %in% includedtypes2,]
# res1=readmutation(data1,dataall)
# res2=readmutation(data2,dataall)

res4=readmutation1(data2,dataall)
# write.table(res4[1:30,],file=output4,sep="\t",quote=F)

mutation_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Missense_Mutation=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  Nonsense_Mutation=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  Silent=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#008000", col = NA))
  },
  Splice_Site=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#018000", col = NA))
  },
  De_novo_Start_InFrame=function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#118000", col = NA))
  },
  De_novo_Start_OutOfFrame = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)

col = c("Missense_Mutation" = "blue", "Nonsense_Mutation"="red", "Silent" = "#008000",
        "Splice_Site"="#018000", "De_novo_Start_InFrame"="#018000", "De_novo_Start_OutOfFrame" = "#008000")

oncoPrint(tmp, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = mutation_fun, col = col, 
          column_title = "EA in Henan",
          heatmap_legend_param = list(title = "Mutations", at = c("Missense_Mutation", "Nonsense_Mutation", "Silent","Splice_Site","De_novo_Start_InFrame","De_novo_Start_OutOfFrame"), 
                                      labels = c("Missense_Mutation", "Nonsense_Mutation", "Silent","Splice_Site","De_novo_Start_InFrame","De_novo_Start_OutOfFrame")))

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

oncoPrint(as.matrix(res4[1:30,]), get_type = function(x) strsplit(x, ";")[[1]],
          row_order = NULL,
          remove_empty_columns = FALSE,
          alter_fun = mutation_fun, col = col, 
          column_title = "EA in Henan",
          heatmap_legend_param = list(title = "Mutations", at = c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site","De_novo_Start_OutOfFrame"), 
                                      labels = c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site","De_novo_Start_OutOfFrame")))


write.table(res1,file=output1,sep="\t",quote=F)
write.table(res2,file=output2,sep="\t",quote=F)
test=read.table(file="41A.maf_keep.annotated.reduced",header=T,sep="\t",quote="")

mat = read.table(paste0(system.file("extdata", package = "ComplexHeatmap"), 
                        "./tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"), 
                 header = TRUE,stringsAsFactors=FALSE, sep = "\t")
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat=  mat[, -ncol(mat)]
mat = t(as.matrix(mat))
mat[1:3, 1:3]

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)

col = c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling",
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "HOMDEL", "MUT"), 
                                      labels = c("Amplification", "Deep deletion", "Mutation")))

