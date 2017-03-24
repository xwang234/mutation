#!/usr/bin/env Rscript
library(GenomicRanges)

#read results from dulak------------------------
# #table to keep SRR->ESO naming transversion
dulaktable=read.table(file="/fh/fast/dai_j/CancerGenomics/EAC/Dulak_fileinfo.txt",header=T)
esotranstable=data.frame(matrix(NA,ncol=4,nrow=nrow(dulaktable)))
dulaktable[,5]=as.character(dulaktable[,5])
colnames(esotranstable)=c("srr","eso","type","wgs")
esotranstable[,1]=as.character(dulaktable[,3])
esotranstable[,4]=dulaktable[,9] #wgs:1,wes:0
for (i in 1:nrow(dulaktable))
{
  tmp=unlist(strsplit(dulaktable[i,5],'-'))
  esotranstable[i,2]=paste0(tmp[[1]],"-",tmp[[2]])
  esotranstable[i,3]=tmp[[3]]
}

dulakstable2=read.table(file="/fh/fast/dai_j/CancerGenomics/EAC/Dulak_Stable2.txt",header=T)
srridx=rep(0,length(grp1tumors))
for (i in 1:length(grp1tumors))
{
  tmp=which(esotranstable[,1]==grp1tumors[i])
  tmp1=esotranstable[tmp,2]
  srridx[i]=which(dulakstable2[,1]==tmp1)
}
#transfer into grp1's order
srrdulakstable2=dulakstable2[srridx,]

dulakstable3=read.table(file="/fh/fast/dai_j/CancerGenomics/EAC/Dulak_Stable3.txt",header=T)
srridx=rep(0,length(grp1tumors))
for (i in 1:length(grp1tumors))
{
  tmp=which(esotranstable[,1]==grp1tumors[i])
  tmp1=esotranstable[tmp,2]
  srridx[i]=which(dulakstable3[,1]==tmp1)
}
#transfer into grp1's order
srrdulakstable3=dulakstable3[srridx,]
#--------------------------

#mutations
grp1normals=c('SRR1002719','SRR999433','SRR999687','SRR1000378','SRR1002786','SRR999559','SRR1001730','SRR10018461','SRR1002703','SRR999599','SRR1002792','SRR1001839','SRR9994281','SRR1002710','SRR9995631','SRR1021476')
grp1tumors=c('SRR1001842','SRR1002713','SRR999423','SRR1001466','SRR1002670','SRR1001823','SRR999489','SRR1002343','SRR1002722','SRR1002656','SRR1002929','SRR999438','SRR1001915','SRR999594','SRR1001868','SRR1001635')

grp2normals=paste0((c(3,11,13,15,17,25,29,33,37,41)+1),"A")
grp2tumors=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")

grp1dir="/fh/scratch/delete30/dai_j/mutect1"
grp2dir="/fh/scratch/delete30/dai_j/henan/mutect3"

#load results
load("mutation/compare_globals_2grps_mutect.RData")
nsvtypes=c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                 "Splice_Site","Translation_Start_Site","Nonstop_Mutation","Targeted_Region","De_novo_Start_InFrame","De_novo_Start_OutOfFrame")

checknummutation=function(grptumors,grpdir)
{
  res=data.frame(matrix(NA,nrow=length(grptumors),ncol=5))
  colnames(res)=c("total","inter","intro","exon","nonsilent")
  rownames(res)=grptumors
  
  for (i in 1:length(grptumors))
  {
    cat(i,"..")
    mutationtable1=read.table(file=paste0(grpdir,"/",grptumors[i],".mutectmaf.ancotator"),header=T,sep="\t",quote="")
    mutationtable=read.table(file=paste0(grpdir,"/",grptumors[i],".variant_function"),sep="\t",quote="")
    res$total[i]=nrow(mutationtable)
    res$inter[i]=sum(mutationtable$V1=="intergenic")
    res$intro[i]=sum(mutationtable$V1 %in% c("intronic","ncRNA_intronic"))
    res$exon[i]=sum(mutationtable$V1 %in% c("exonic","ncRNA_exonic"))
    res$nonsilent[i]=sum(mutationtable1$Variant_Classification %in% nsvtypes)
  }
  return(res)
}
grp1mutations=checknummutation(grp1tumors,grp1dir)
grp2mutations=checknummutation(grp2tumors,grp2dir)

chrs=paste0("chr",c(1:22,"X","Y"))
intron1=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/refgenes_intron.txt",sep="\t",quote="")
intron1=intron1[intron1$V1 %in% chrs,]
gr_intron1=GRanges(seqnames = intron1$V1,ranges=IRanges(start=intron1$V2,end=intron1$V3))
gr_intron1=reduce(gr_intron1)
sum(as.numeric(width(gr_intron1)))

exon1=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/refgenes_exon.txt",sep="\t",quote="")
exon1=exon1[exon1$V1 %in% chrs,]
gr_exon1=GRanges(seqnames = exon1$V1,ranges=IRanges(start=exon1$V2,end=exon1$V3))
gr_exon1=reduce(gr_exon1)
sum(width(gr_exon1))

coding1=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/refgenes_coding.txt",sep="\t",quote="")
coding1=coding1[coding1$V1 %in% chrs,]
gr_coding1=GRanges(seqnames = coding1$V1,ranges=IRanges(start=coding1$V2,end=coding1$V3))
gr_coding1=reduce(gr_coding1)
sum(width(gr_coding1))

wholegene1=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/refgenes_wholegene.txt",sep="\t",quote="")
wholegene1=wholegene1[wholegene1$V1 %in% chrs,]
gr_wholegene1=GRanges(seqnames = wholegene1$V1,ranges=IRanges(start=wholegene1$V2,end=wholegene1$V3))
gr_wholegene1=reduce(gr_wholegene1)
sum(width(gr_wholegene1))

gr_intron_exon1=intersect(gr_intron1,gr_exon1)
sum(width(gr_intron_exon1))

#exon +#intron =wholegen, exon overlap with intron
countmutationrate=function(grpmuations,grptumors,grpdir)
{
  mutationrate=data.frame(matrix(NA,nrow=length(grptumors),ncol=5))
  colnames(mutationrate)=c("total","total1","inter","intro","exon")
  len=data.frame(matrix(NA,nrow=length(grptumors),ncol=4))
  colnames(len)=c("total","inter","intro","exon")
  
  for (i in 1:length(grptumors))
  {
    cat(i,"..")
    wigfile=paste0(grpdir,"/",grptumors[i],".wig.bed.txt")
    wigtable=read.table(wigfile,sep="\t")
    wigtable=wigtable[wigtable$V1 %in% chrs,]
    gr_wigtable=GRanges(seqnames = wigtable$V1,ranges = IRanges(start=wigtable$V2,end=wigtable$V3))
    len$total[i]=sum(as.numeric(width(gr_wigtable)))
    mutationrate$total[i]=grpmuations$total[i]/len$total[i]*10^6
    mutationrate$total1[i]=srrdulakstable2$GenomicMutations[i]/len$total[i]*10^6 #use the total count from dulak
    
    tmp=intersect(gr_intron,gr_wigtable)
    tmp2=intersect(gr_intron_exon,gr_wigtable)
    len$intro[i]=sum(as.numeric(width(tmp)))-sum(as.numeric(width(tmp2)))
    len$intro[i]=len$intro[i]*0.9
    
    tmp1=intersect(gr_coding,gr_wigtable)
    len$exon[i]=sum(width(tmp1))*1.1
    
    tmp3=intersect(gr_wholegene,gr_wigtable)
    len$inter[i]=len$total[i]-sum(as.numeric(width(tmp3)))
    len$inter[i]=len$inter[i]*0.9
    
    mutationrate$inter[i]=grpmuations$inter[i]/len$inter[i]*10^6
    mutationrate$intro[i]=grpmuations$intro[i]/len$intro[i]*10^6
    mutationrate$exon[i]=grpmuations$exon[i]/len$exon[i]*10^6
    #
  }
  return(result=list(len=len,mutationrate=mutationrate))
}
grp1mutationrate=countmutationrate(grp1mutations,grp1tumors,grpdir="/fh/scratch/delete30/dai_j/mutect3")
grp2mutationrate=countmutationrate(grp2mutations,grp2tumors,grpdir=grp2dir)

median(grp1mutations$total)
range(grp1mutations$total)
median(grp2mutations$total)
range(grp2mutations$total)
median(grp1mutationrate$mutationrate$total)
range(grp1mutationrate$mutationrate$total)
median(grp2mutationrate$mutationrate$total)
range(grp2mutationrate$mutationrate$total)

median(grp1mutationrate$mutationrate$total)
range(grp1mutationrate$mutationrate$total)
median(grp2mutationrate$mutationrate$total)
range(grp2mutationrate$mutationrate$total)

median(grp2mutationrate$mutationrate$inter)
range(grp2mutationrate$mutationrate$inter)
median(grp2mutationrate$mutationrate$intro)
range(grp2mutationrate$mutationrate$intro)
median(grp2mutationrate$mutationrate$exon)
range(grp2mutationrate$mutationrate$exon)
median(grp2mutations$nonsilent/grp2mutationrate$len$nonsilent)*10^6
range((grp2mutations$nonsilent/grp2mutationrate$len$nonsilent)*10^6)

restable=data.frame(matrix(NA,nrow=length(grp1tumors)+length(grp2tumors),ncol=2))
colnames(restable)=c("dataset","sample")
restable[1:length(grp1tumors),1]="US-EA"
restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),1]="CH-EA"
restable[,2]=c(grp1tumors,grp2tumors)
#reorder group
restable$dataset=factor(restable$dataset,c("US-EA","CH-EA"))
color=c("green","blue","skyblue")
restable=cbind(restable,totalrate=c(grp1mutationrate$mutationrate$total,grp2mutationrate$mutationrate$total))
outfig="EA_muationfrequency.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
plot2groups(item="totalrate",datatable=restable,adjustymin=T,main='Mutation frequency',usewilcox=T)
dev.off()

draw2barplot=function(data1,data2,ylable="")
{
  tmp=data.frame(data1=data1,data2=data2)
  tmp=t(tmp)
  barplot(as.matrix(tmp),beside=TRUE,col=c("blue","green"),ylab=ylable)
  legend("topright",legend=c("dulak","replicate"),fill=c("blue","green"))
}

draw2barplot(srrdulakstable2$MutationsperMb,grp1mutationrate$mutationrate$total,ylable ="total mutation rate" )
draw2barplot(srrdulakstable3$IntronMutationRate,grp1mutationrate$mutationrate$intro,ylable ="intron mutation rate" )
draw2barplot(srrdulakstable3$IGRMutationRate,grp1mutationrate$mutationrate$inter,ylable ="intergenic mutation rate" )
draw2barplot(srrdulakstable3$ExonMutationRate,grp1mutationrate$mutationrate$exon,ylable ="exon mutation rate" )
t.test(srrdulakstable2$MutationsperMb,grp1mutationrate$mutationrate$total)$p.value
t.test(srrdulakstable3$IntronMutationRate,grp1mutationrate$mutationrate$intro)$p.value
t.test(srrdulakstable3$IGRMutationRate,grp1mutationrate$mutationrate$inter)$p.value
t.test(srrdulakstable3$ExonMutationRate,grp1mutationrate$mutationrate$exon)$p.value
