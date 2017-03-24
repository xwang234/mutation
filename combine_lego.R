#!/usr/bin/env Rscript

#henan
#mutationdir="/fh/scratch/delete30/dai_j/henan/varscan1"
mutationdir="/fh/scratch/delete30/dai_j/henan/varscan2"
mutationdir="/fh/scratch/delete30/dai_j/henan/mutect3"
output="henan_golden4_varscan_lego.txt"
output="henan_mutect3_lego.txt"
wgstumors=c("3A","11A","13A","15A")
wgstumors=c("3A","11A","13A","15A","17A","25A","29A","33A","37A","41A")

#escc
mutationdir="/fh/scratch/delete30/dai_j/escc/varscan2"
output="escc_varscan_lego.txt"
wgstumors=c("T1","T2","T3","T4","T5","T6","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17","T18")

#dulak
mutationdir="/fh/scratch/delete30/dai_j/varscan2"
output="dulak_varscan_lego.txt"
wgstumors=c("SRR1001842","SRR1002713","SRR999423","SRR1001466","SRR1002670","SRR1001823","SRR999489","SRR1002343","SRR1002722","SRR1002656",
            "SRR1002929","SRR999438","SRR1001915","SRR999594","SRR1001868","SRR1001635")

mutationdir="/fh/scratch/delete30/dai_j/mutect1"
output="dulak_mutect1_lego.txt"

for (i in 1:length(wgstumors))
{
  #legofile=paste0(mutationdir,"/",wgstumors[i],".snp.Somatic.hc.lego.txt")
  legofile=paste0(mutationdir,"/",wgstumors[i],".lego.txt")
  tmp=read.table(legofile,header=T,sep="\t")
  if (i==1)
  {
    res=tmp
  }else
  {
    res=res+tmp
  }
}
write.table(res,file=paste0(mutationdir,"/",output),row.names = F,sep="\t",quote=F)
