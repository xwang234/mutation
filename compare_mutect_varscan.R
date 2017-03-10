#!/usr/bin/env Rscript

wgstumor="13A"
mutectfile=paste0("/fh/scratch/delete30/dai_j/henan/mutect/",wgstumor,".chr1.Mutect_out.txt")
mutectkeepfile=paste0("/fh/scratch/delete30/dai_j/henan/mutect/",wgstumor,".Mutect_out_keep.txt")
varscanfile=paste0("/fh/scratch/delete30/dai_j/henan/varscan1/",wgstumor,"_chr1.snp")
varscankeepfile=paste0("/fh/scratch/delete30/dai_j/henan/varscan1/",wgstumor,".snp.Somatic.hc")

library(GenomicRanges)
compare_mutect_varscan=function(mutectfile,mutectkeepfile,varscanfile,varscankeepfile,chr="chr1")
{
  mutectable=read.table(file=mutectfile,header=F,sep="\t",quote="")
  mutectkeeptable=read.table(file=mutectkeepfile,header=T,sep="\t")
  colnames(mutectable)=colnames(mutectkeeptable)
  leftcols=c("contig","position","ref_allele","alt_allele","dbsnp_site","covered","power","tumor_power","normal_power",
             "total_pairs","map_Q0_reads","init_t_lod","t_lod_fstar","contaminant_lod","t_q20_count","t_ref_count","t_alt_count",
             "tumor_f","t_ins_count","t_del_count","n_ref_count","n_alt_count","normal_f","failure_reasons","judgement")
  idx=which(colnames(mutectable) %in% leftcols)
  mutectable=mutectable[,idx]
  colnames(mutectable)[1]=colnames(mutectkeeptable)[1]="chrom"
  GR_mutectable=GRanges(seqnames = mutectable$chrom,ranges=IRanges(start=mutectable$position,width=1))
  #mutect fails coverage filter
  sum((mutectable$t_ref_count+mutectable$t_alt_count<10 | mutectable$n_ref_count+mutectable$n_alt_count<10 | mutectable$tumor_f<0.1) & mutectable$judgement=="KEEP")
  varscantable=read.table(file=varscanfile,header=T,sep="\t",quote="")
  varscankeeptable=read.table(file=varscankeepfile,header=T,sep="\t")
  varscankeeptable=varscankeeptable[varscankeeptable$chrom==chr,]
  mutectkeeptable=mutectkeeptable[mutectkeeptable$chrom==chr,]
  GR_mutectkeeptable=GRanges(seqnames=mutectkeeptable$chrom,IRanges(start=mutectkeeptable$position,width=1))
  GR_varscankeeptable=GRanges(seqnames = varscankeeptable$chrom,IRanges(start=varscankeeptable$position,width=1))
  num_mutect=nrow(mutectkeeptable)
  num_varscan=nrow(varscankeeptable)
  #check the overlap
  idx=which(GR_varscankeeptable %in% GR_mutectkeeptable)
  num_intersect=length(idx)
  varscan_intesect=varscankeeptable[idx,]
  #what are the reasons varscan mutations not found in mutect
  idx=which(! GR_varscankeeptable %in% GR_mutectkeeptable)
  GR_uniq_varscankeeptable=GR_varscankeeptable[idx,]
  uniq_varscankeeptable=varscankeeptable[idx,]
  idx=which(!GR_uniq_varscankeeptable %in% GR_mutectable)
  GR_uniq_varscankeep_notin_mutect=GR_uniq_varscankeeptable[idx,]
  uniq_varscankeep_notin_mutect=uniq_varscankeeptable[idx,]
  #varscan mutations not found in mutect_out table
  num_varscannotfoundinmutect=nrow(uniq_varscankeep_notin_mutect)
  idx=which(GR_mutectable %in% GR_uniq_varscankeeptable)
  #if can be found in mutect_out table, what are the reasons they are filtered out
  reasons=data.frame(matrix(NA,nrow=length(idx),ncol=20))
  for (i in 1:length(idx))
  {
    tmp=as.character(mutectable$failure_reasons[idx[i]])
    tmp=unlist(strsplit(tmp,",",fixed=T))
    tmp=tmp[order(tmp)]
    reasons[i,1:length(tmp)]=tmp
  }
  
  reasonscount=data.frame(matrix(NA,nrow=1,ncol=length(unique(unlist(reasons)))-1))
  colnames(reasonscount)=names(table(unlist(reasons)))
  reasonscount[1,]=table(unlist(reasons))
  reasonscount=cbind.data.frame(num_varscan=num_varscan,num_inmutect=num_intersect,num_notfoundinmutectout=num_varscannotfoundinmutect,
                                num_filteredout=num_varscan-num_intersect-num_varscannotfoundinmutect,reasonscount)
  
}
res_13A_chr1=compare_mutect_varscan(mutectfile,mutectkeepfile,varscanfile,varscankeepfile,chr="chr1")

wgstumor="SRR1001466"
mutectfile=paste0("/fh/scratch/delete30/dai_j/mutect1/",wgstumor,".chr1.Mutect_out.txt")
mutectkeepfile=paste0("/fh/scratch/delete30/dai_j/mutect1/",wgstumor,".Mutect_out_keep.txt")
varscanfile=paste0("/fh/scratch/delete30/dai_j/varscan2/",wgstumor,"_chr1.snp")
varscankeepfile=paste0("/fh/scratch/delete30/dai_j/varscan2/",wgstumor,".snp.Somatic.hc")
res_SRR1001466_chr1=compare_mutect_varscan(mutectfile,mutectkeepfile,varscanfile,varscankeepfile,chr="chr1")


