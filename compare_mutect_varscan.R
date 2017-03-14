#!/usr/bin/env Rscript

wgstumor="3A"
mutectfile=paste0("/fh/scratch/delete30/dai_j/henan/mutect/",wgstumor,".chr1.Mutect_out.txt") #with title,generate_mutect_out_chr.sh
varscanfile=paste0("/fh/scratch/delete30/dai_j/henan/varscan1/",wgstumor,"_chr1.snp")
varscankeepfile=paste0("/fh/scratch/delete30/dai_j/henan/varscan1/",wgstumor,".snp.Somatic.hc")

library(GenomicRanges)
compare_mutect_varscan=function(mutectfile,varscanfile,varscankeepfile,chr="chr1")
{
  mutectable=read.table(file=mutectfile,header=T,sep="\t",quote="")
  
  leftcols=c("contig","position","ref_allele","alt_allele","dbsnp_site","covered","power","tumor_power","normal_power",
             "total_pairs","map_Q0_reads","init_t_lod","t_lod_fstar","contaminant_lod","t_q20_count","t_ref_count","t_alt_count",
             "tumor_f","t_ref_sum","t_alt_sum","t_ins_count","t_del_count","n_ref_count","n_alt_count","t_ref_sum","t_alt_sum","normal_f","failure_reasons","judgement")
  idx=which(colnames(mutectable) %in% leftcols)
  mutectable=mutectable[,idx]
  colnames(mutectable)[1]="chrom"
  GR_mutectable=GRanges(seqnames = mutectable$chrom,ranges=IRanges(start=mutectable$position,width=1))
  #mutect fails coverage filter
  sum((mutectable$t_ref_count+mutectable$t_alt_count<10 | mutectable$n_ref_count+mutectable$n_alt_count<10 | mutectable$tumor_f<0.1) & mutectable$judgement=="KEEP")
  varscantable=read.table(file=varscanfile,header=T,sep="\t",quote="")
  varscankeeptable=read.table(file=varscankeepfile,header=T,sep="\t")
  varscankeeptable=varscankeeptable[varscankeeptable$chrom==chr,]
  mutectkeeptable=mutectable[mutectable$judgement=="KEEP",]
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
  if ("failure_reasons" %in% colnames(mutectable))
  {
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
  
  #check the readdepth
  failedmutect=mutectable[idx,]
  GR_failedmutect=GRanges(seqnames=failedmutect$chrom,ranges=IRanges(start=failedmutect$position,width=1))
  idx1=which(GR_varscankeeptable %in% GR_failedmutect)
  failedvarscan=varscankeeptable[idx1,]
  readdepth=data.frame(matrix(NA,nrow=length(idx),ncol=10))
  colnames(readdepth)=c("chrom","position","mutect_n_ref","mutect_n_alt","mutect_t_ref","mutect_t_alt",
                        "varscan_n_ref","varscan_n_alt","varscan_t_ref","varscan_t_alt")
  readdepth$chrom=failedmutect$chrom
  readdepth$position=failedmutect$position
  for (i in 1:length(idx))
  {
    readdepth$mutect_n_ref[i]=failedmutect$n_ref_count[i]
    readdepth$mutect_n_alt[i]=failedmutect$n_alt_count[i]
    readdepth$mutect_t_ref[i]=failedmutect$t_ref_count[i]
    readdepth$mutect_t_alt[i]=failedmutect$t_alt_count[i]
    readdepth$varscan_n_ref[i]=failedvarscan$normal_reads1[i]
    readdepth$varscan_n_alt[i]=failedvarscan$normal_reads2[i]
    readdepth$varscan_t_ref[i]=failedvarscan$tumor_reads1[i]
    readdepth$varscan_t_alt[i]=failedvarscan$tumor_reads2[i]
  }
  if ("failure_reasons" %in% colnames(mutectable))
  {
    results=list(failedmutect=failedmutect,reasonscount=reasonscount,readdepth=readdepth,mutectkeeptable=mutectkeeptable,
                 varscankeeptable=varscankeeptable,varscan_intesect=varscan_intesect)
  }else
  {
    results=list(failedmutect=failedmutect,readdepth=readdepth,mutectkeeptable=mutectkeeptable,
                 varscankeeptable=varscankeeptable,varscan_intesect=varscan_intesect)
  }
  return(results)
}
res_3A_chr11=compare_mutect_varscan(mutectfile,varscanfile,varscankeepfile,chr="chr1")

wgstumor="SRR1001466"
mutectfile=paste0("/fh/scratch/delete30/dai_j/mutect1/",wgstumor,".chr1.Mutect_out.txt")
varscanfile=paste0("/fh/scratch/delete30/dai_j/varscan2/",wgstumor,"_chr1.snp")
varscankeepfile=paste0("/fh/scratch/delete30/dai_j/varscan2/",wgstumor,".snp.Somatic.hc")
res_SRR1001466_chr1=compare_mutect_varscan(mutectfile,varscanfile,varscankeepfile,chr="chr1")

check_varscan_normalvar=function(varscankeepfile)
{
  varscankeeptable=read.table(file=varscankeepfile,header = T,sep="\t")
  n_alt_count=varscankeeptable$normal_reads2
  normal_f=as.numeric(gsub("%","",varscankeeptable$normal_var_freq,fixed = T))
  result=data.frame(n_alt_count=n_alt_count,normal_f=normal_f)
}

wgstumor="13A"
varscankeepfile=paste0("/fh/scratch/delete30/dai_j/henan/varscan1/",wgstumor,".snp.Somatic.hc")
normal_13A=check_varscan_normalvar(varscankeepfile)

wgstumor="SRR1001466"
varscankeepfile=paste0("/fh/scratch/delete30/dai_j/varscan2/",wgstumor,".snp.Somatic.hc")
normal_SRR1001466=check_varscan_normalvar(varscankeepfile)
sum(normal_13A$n_alt_count<=2)/nrow(normal_13A)
#[1] 0.9021018
sum(normal_SRR1001466$n_alt_count<=2)/nrow(normal_SRR1001466)
#[1] 0.9720614
sum(normal_13A$normal_f<=3)/nrow(normal_13A)
#[1] 0.8549866
sum(normal_SRR1001466$normal_f<=3)/nrow(normal_SRR1001466)
#[1] 0.9163772

#check mapping 0 
m_1001466=read.table(file="/fh/scratch/delete30/dai_j/dulak/test_SRR1001466.mpileup.txt",sep="\t",quote="",fill = T) #all pos
m1_1001466=read.table(file="/fh/scratch/delete30/dai_j/dulak/testq1_SRR1001466.mpileup.txt",sep="\t",quote="",fill = T) #pos with q>=1
m_13A=read.table(file="/fh/scratch/delete30/dai_j/henan/test_13A.mpileup.txt",sep="\t",quote="",fill=T)
m1_13A=read.table(file="/fh/scratch/delete30/dai_j/henan/testq1_13A.mpileup.txt",sep="\t",quote="",fill=T)

myintersect=function(table1,table2)
{
  GR_table1=GRanges(seqnames = table1[,1],ranges = IRanges(start=table1[,2],width = 1))
  GR_table2=GRanges(seqnames = table2[,1],ranges = IRanges(start=table2[,2],width = 1))
  idx1=which(GR_table1 %in% GR_table2)
  intersect_table=table1[idx1,]
  GR_intersect=GRanges(seqnames = intersect_table[,1],ranges = IRanges(start=intersect_table[,2],width=1))
  idx1=which(GR_table1 %in% GR_intersect)
  idx2=which(GR_table2 %in% GR_intersect)
  return(results=list(idx1=idx1,idx2=idx2))
}

idx_1001466=myintersect(table1 = m_1001466,table2 = m1_1001466)
dif_1001466=m_1001466[idx_1001466$idx1,4]-m1_1001466[idx_1001466$idx2,4]
idx_13A=myintersect(table1 = m_13A,table2 = m1_13A)
dif_13A=m_13A[idx_13A$idx1,4]-m1_13A[idx_13A$idx2,4]
