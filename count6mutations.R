#! /usr/bin/env Rscript
#SBATCH -t 0-10
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

#work on a single table (for one sample)
countmutations=function(maffile,header=F,colchr,colstart,colref,colalt,pref)
{
  library(GenomicRanges)
  library(BSgenome)
  library("BSgenome.Hsapiens.UCSC.hg19")
  maftable=read.table(maffile,sep="\t",header=header)
  chrs=paste0("chr",c(1:22,"X","Y"))

  if (!grepl("chr",maftable[1,colchr]))
  {
    maftable[,colchr]=paste0("chr",maftable[,colchr])
    #mychr<-array(paste0("chr",maftable[,colchr]))
  }
  
  idx=maftable[,colchr] %in% chrs
  maftable=maftable[idx,]
  
  trinuc=c("T","C","A","G")
  ref=maftable[,colref]
  alt=maftable[,colalt]
  idx=ref %in% trinuc & alt %in% trinuc #only works on snv
  maftable=maftable[idx,]
  
  mychr<-array(maftable[,colchr])
  mystart<-array(maftable[,colstart])
  mystart=as.integer(mystart)
  myref<-array(maftable[,colref])
  myalt<-array(maftable[,colalt])
  myseq=array(NA,length(mystart))
  myrange=GRanges(seqnames = mychr,ranges=IRanges(start=mystart-1,end=mystart+1))
  genome <- BSgenome.Hsapiens.UCSC.hg19
  tmp=getSeq(genome, myrange)
  
  for (i in 1:length(mystart))
  {
    myseq[i]=toString(tmp[[i]])
  }

  
  ref1=c("C","C","C","A","A","A")
  alt1=c("T","A","G","G","C","T")
  ref2=c("G","G","G","T","T","T")
  alt2=c("A","T","C","C","G","A")
  #use the mutation format used in lego.R
  types=c("C.T.G.A","C.A.G.T","C.G.G.C","T.C.A.G","T.G.A.C","T.A.A.T")
  #count the mutations of AA to AC
  res=data.frame(matrix(0,nrow=4*6,ncol=4))
  #save lego data
  lego=data.frame(matrix(NA,ncol=96,nrow=1))
  row=0
  transtable=data.frame(matrix(NA,ncol=3,nrow=96)) #class,chr,start
  count1=0
  for (j in 1:6) {
    
    for (m in 1:length(trinuc)) {
      pre1=trinuc[m] #the nucleotide before mutation site
      for (n in 1:length(trinuc)) {
        post1=trinuc[n] #after mutation site
        count=0
        for (i in 1:length(mychr)) {
          if (substr(myseq[i],1,1) == pre1 & substr(myseq[i],3,3) == post1 & myref[i] == ref1[j] & myalt[i] == alt1[j]) 
          { 
            count <-count+1
            count1=count1+1
            transtable[count1,1]=paste0(types[j],"_",trinuc[m],"x",trinuc[n])
            transtable[count1,2]=mychr[i]
            transtable[count1,3]=mystart[i]
          }
          if (substr(myseq[i],1,1) == pre1 & substr(myseq[i],3,3) == post1 & myref[i] == ref2[j] & myalt[i] == alt2[j]) 
          { 
            count <-count+1
            count1=count1+1
            transtable[count1,1]=paste0(types[j],"_",trinuc[m],"x",trinuc[n])
            transtable[count1,2]=mychr[i]
            transtable[count1,3]=mystart[i]
          }
        }
        res[row+m,n]=count
        lego[1,(j-1)*16+(m-1)*4+n]=count
        colnames(lego)[(j-1)*16+(m-1)*4+n]=paste0(types[j],"_",trinuc[m],"x",trinuc[n])
      }
      rownames(res)[row+m]=paste0(trinuc[m],j)
    }
    row=row+4
  }
  colnames(res)=trinuc
  colnames(transtable)=c("type","chr","pos")
  c2t=sum(res[1:4,])
  c2t_rate=c2t/sum(res)
  c2a=sum(res[5:8,])
  c2a_rate=c2a/sum(res)
  c2g=sum(res[9:12,])
  c2g_rate=c2g/sum(res)
  a2g=sum(res[13:16,])
  a2g_rate=a2g/sum(res)
  a2c=sum(res[17:20,])
  a2c_rate=a2c/sum(res)
  a2t=sum(res[21:24,])
  a2t_rate=a2t/sum(res)
  aa2c=sum(res[19,])
  aa2c_rate=aa2c/sum(res)
  resrate=data.frame(matrix(0,ncol=1,nrow=7))
  resrate[1,]=c2t_rate
  resrate[2,]=c2a_rate
  resrate[3,]=c2g_rate
  resrate[4,]=a2g_rate
  resrate[5,]=a2c_rate
  resrate[6,]=a2t_rate
  resrate[7,]=aa2c_rate
  rownames(resrate)=c('c2t','c2a','c2g','a2g','a2c','a2t','aa2c')
  colnames(resrate)=c("rate")
  results=list(convmat=res,covrate=resrate,lego=lego,transtable=transtable)
  
  output1=paste0(pref,".transmat.txt")
  output2=paste0(pref,".transrate.txt")
  output3=paste0(pref,".lego.txt")
  output4=paste0(pref,".transtable.txt")
  write.table(results$convmat,file=output1,na="NA",sep="\t",row.names=T,col.names=T,quote=F)
  write.table(results$covrate,file=output2,na="NA",sep="\t",row.names=T,col.names=T,quote=F)
  write.table(results$lego,file=output3,na="NA",sep="\t",row.names=F,col.names=T,quote=F)
  write.table(results$transtable,file=output4,sep="\t",row.names=F,col.names=T,quote=F)
  
  return(results)
  
}

#work on a combined data (all the samples in a table)
countmutations2=function(maffile,header=F,colchr,colstart,colref,colalt,outfolder,colsample=1)
{
  print("countmutations2")
  print(maffile)
  library(GenomicRanges)
  library(BSgenome)
  library("BSgenome.Hsapiens.UCSC.hg19")
  allmaftable=read.table(maffile,sep="\t",header=header)
  chrs=paste0("chr",c(1:22,"X","Y"))
  if (!grepl("chr",allmaftable[1,colchr]))
  {
    allmaftable[,colchr]=paste0("chr",allmaftable[,colchr])
    #mychr<-array(paste0("chr",maftable[,colchr]))
  }
  
  idx=allmaftable[,colchr] %in% chrs
  allmaftable=allmaftable[idx,]
  
  trinuc=c("T","C","A","G")
  ref1=c("C","C","C","A","A","A")
  alt1=c("T","A","G","G","C","T")
  ref2=c("G","G","G","T","T","T")
  alt2=c("A","T","C","C","G","A")
  #use the mutation format used in lego.R
  types=c("C.T.G.A","C.A.G.T","C.G.G.C","T.C.A.G","T.G.A.C","T.A.A.T")
  #count the mutations of AA to AC
  
  samples=unique(allmaftable[,colsample])
  write.table(samples,file=paste0(outfolder,'/','samplenames.txt'),sep="\t",row.names = F,col.names = F,quote=F)
  genome <- BSgenome.Hsapiens.UCSC.hg19
  for (numpair in 1:length(samples))
  {
    idx=which(allmaftable[,colsample]==samples[numpair])
    maftable=allmaftable[idx,]
    
    mychr<-array(maftable[,colchr])
    mystart<-array(maftable[,colstart])
    mystart=as.integer(mystart)
    myref<-array(maftable[,colref])
    myalt<-array(maftable[,colalt])
    myseq=array(NA,length(mystart))
    myrange=GRanges(seqnames = mychr,ranges=IRanges(start=mystart-1,end=mystart+1))
    
    tmp=getSeq(genome, myrange)
    
    for (i in 1:length(mystart))
    {
      myseq[i]=toString(tmp[[i]])
    }
    res=data.frame(matrix(0,nrow=4*6,ncol=4))
    #save lego data
    lego=data.frame(matrix(NA,ncol=96,nrow=1))
    row=0
    transtable=data.frame(matrix(NA,ncol=3,nrow=96)) #class,chr,start
    count1=0
    for (j in 1:6) {
      
      for (m in 1:length(trinuc)) {
        pre1=trinuc[m] #the nucleotide before mutation site
        for (n in 1:length(trinuc)) {
          post1=trinuc[n] #after mutation site
          count=0
          for (i in 1:length(mychr)) {
            if (substr(myseq[i],1,1) == pre1 & substr(myseq[i],3,3) == post1 & myref[i] == ref1[j] & myalt[i] == alt1[j]) 
            { 
              count <-count+1
              count1=count1+1
              transtable[count1,1]=paste0(types[j],"_",trinuc[m],"x",trinuc[n])
              transtable[count1,2]=mychr[i]
              transtable[count1,3]=mystart[i]
            }
            if (substr(myseq[i],1,1) == pre1 & substr(myseq[i],3,3) == post1 & myref[i] == ref2[j] & myalt[i] == alt2[j]) 
            { 
              count <-count+1
              count1=count1+1
              transtable[count1,1]=paste0(types[j],"_",trinuc[m],"x",trinuc[n])
              transtable[count1,2]=mychr[i]
              transtable[count1,3]=mystart[i]
            }
          }
          res[row+m,n]=count
          lego[1,(j-1)*16+(m-1)*4+n]=count
          colnames(lego)[(j-1)*16+(m-1)*4+n]=paste0(types[j],"_",trinuc[m],"x",trinuc[n])
        }
        rownames(res)[row+m]=paste0(trinuc[m],j)
      }
      row=row+4
    }
    colnames(res)=trinuc
    colnames(transtable)=c("type","chr","pos")
    c2t=sum(res[1:4,])
    c2t_rate=c2t/sum(res)
    c2a=sum(res[5:8,])
    c2a_rate=c2a/sum(res)
    c2g=sum(res[9:12,])
    c2g_rate=c2g/sum(res)
    a2g=sum(res[13:16,])
    a2g_rate=a2g/sum(res)
    a2c=sum(res[17:20,])
    a2c_rate=a2c/sum(res)
    a2t=sum(res[21:24,])
    a2t_rate=a2t/sum(res)
    aa2c=sum(res[19,])
    aa2c_rate=aa2c/sum(res)
    resrate=data.frame(matrix(0,ncol=1,nrow=7))
    resrate[1,]=c2t_rate
    resrate[2,]=c2a_rate
    resrate[3,]=c2g_rate
    resrate[4,]=a2g_rate
    resrate[5,]=a2c_rate
    resrate[6,]=a2t_rate
    resrate[7,]=aa2c_rate
    rownames(resrate)=c('c2t','c2a','c2g','a2g','a2c','a2t','aa2c')
    colnames(resrate)=c("rate")
    results=list(convmat=res,covrate=resrate,lego=lego,transtable=transtable)
    output1=paste0(outfolder,'/',samples[numpair],".transmat.txt")
    output2=paste0(outfolder,'/',samples[numpair],".transrate.txt")
    output3=paste0(outfolder,'/',samples[numpair],".lego.txt")
    output4=paste0(outfolder,'/',samples[numpair],".transtable.txt")
    write.table(results$convmat,file=output1,na="NA",sep="\t",row.names=T,col.names=T,quote=F)
    write.table(results$covrate,file=output2,na="NA",sep="\t",row.names=T,col.names=T,quote=F)
    write.table(results$lego,file=output3,na="NA",sep="\t",row.names=F,col.names=T,quote=F)
    write.table(results$transtable,file=output4,sep="\t",row.names=F,col.names=T,quote=F)
    
  }
}




args <- commandArgs(trailingOnly = TRUE)
maffile=as.character(args[1])
colchr=as.integer(args[2])
colstart=as.integer(args[3])
colref=as.integer(args[4])
colalt=as.integer(args[5])
header=as.logical(args[6])
pref=as.character(args[7])
opt=as.character(args[8])

#colchr=4
#colstart=5
#colref=7
#colalt=8
if (opt==1)
{
  countmutations2(maffile,header=header,colchr=colchr,colstart=colstart,colref=colref,colalt=colalt,outfolder=pref)
}else
{
  results=countmutations(maffile,header=header,colchr=colchr,colstart=colstart,colref=colref,colalt=colalt,pref=pref)
}

#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/TCGA/Somatic_Mutations/BCGSC__IlluminaHiSeq_DNASeq_automated/Level_2/bcgsc.ca__IlluminaHiSeq_automated_DNA_sequencing_level2.maf 5 6 11 13 T BCGSC

#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/TCGA/Somatic_Mutations/BCM__IlluminaGA_DNASeq_automated/Level_2/hgsc.bcm.edu__IlluminaGA_automated_DNA_sequencing_level2.maf 5 6 11 13 T BCM

#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/TCGA/Somatic_Mutations/BI__IlluminaGA_DNASeq_automated/Level_2/broad.mit.edu__IlluminaGA_automated_DNA_sequencing_level2.maf 5 6 11 13 T BI

#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/TCGA/Somatic_Mutations/UCSC__IlluminaGA_DNASeq_automated/Level_2/ucsc.edu__IlluminaGA_automated_DNA_sequencing_level2.maf 5 6 11 13 T UCSC

#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/TCGA/Somatic_Mutations/WUSM__IlluminaHiSeq_DNASeq_automated/Level_2/genome.wustl.edu__IlluminaHiSeq_automated_DNA_sequencing_level2.maf 5 6 11 13 T WUSM

#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/1A.variant_function 3 4 6 7 F 1A

#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/1A.exonic_variant_function 4 5 7 8 F 1A.exonic

#12/5
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/17A.variant_function 3 4 6 7 F 17A 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/41A.variant_function 3 4 6 7 F 41A 0

#12/9
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/13A.variant_function 3 4 6 7 F 13A 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/15A.variant_function 3 4 6 7 F 15A 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/23A.variant_function 3 4 6 7 F 23A 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/25A.variant_function 3 4 6 7 F 25A 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/29A.variant_function 3 4 6 7 F 29A 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/33A.variant_function 3 4 6 7 F 33A 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/35A.variant_function 3 4 6 7 F 35A 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/37A.variant_function 3 4 6 7 F 37A 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/mutect/41A.variant_function 3 4 6 7 F 41A 0

#12/16
#process varscan somatic mutation result
varscandir="/fh/scratch/delete30/dai_j/henan/varscan"
tumors=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")
chrs=paste0("chr",c(1:22,"X","Y"))
process_varscan_somatic=function(varscandir,tumors)
{
  for (numpair in 1:length(tumors))
  {
    tumor=tumors[numpair]
    maffile=paste0(varscandir,'/',tumor,".somatic.snp")
    tmp=read.table(maffile,header = T)
    tmp=tmp[tmp$chrom %in% chrs,]
    tmp=tmp[tmp$normal_reads1+tmp$normal_reads2>=10 & tmp$tumor_reads1+tmp$tumor_reads2>=10 & tmp$somatic_status=="Somatic",]
    tmp$normal_var_freq=as.character(tmp$normal_var_freq)
    tmp$normal_var_freq=as.numeric(gsub("%","",tmp$normal_var_freq,fixed=T))
    tmp$tumor_var_freq=as.character(tmp$tumor_var_freq)
    tmp$tumor_var_freq=as.numeric(gsub("%","",tmp$tumor_var_freq,fixed=T))
    tmp=tmp[tmp$normal_var_freq<=5,] #less than 5%
    write.table(tmp,file=paste0(varscandir,"/",tumor,".somatic_keep.txt"),row.names=F,col.names = T,sep="\t",quote=F ) 
  }
  
}
#process_varscan_somatic(varscandir,tumors)
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan/3A.somatic_keep.txt 1 2 3 4 T 3A.varscan_somatic 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan/11A.somatic_keep.txt 1 2 3 4 T 11A.varscan_somatic 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan/13A.somatic_keep.txt 1 2 3 4 T 13A.varscan_somatic 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan/15A.somatic_keep.txt 1 2 3 4 T 15A.varscan_somatic 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan/17A.somatic_keep.txt 1 2 3 4 T 17A.varscan_somatic 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan/23A.somatic_keep.txt 1 2 3 4 T 23A.varscan_somatic 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan/25A.somatic_keep.txt 1 2 3 4 T 25A.varscan_somatic 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan/33A.somatic_keep.txt 1 2 3 4 T 33A.varscan_somatic 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan/35A.somatic_keep.txt 1 2 3 4 T 35A.varscan_somatic 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan/37A.somatic_keep.txt 1 2 3 4 T 37A.varscan_somatic 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan/41A.somatic_keep.txt 1 2 3 4 T 41A.varscan_somatic 0


#1/18
#process_varscan_somatic(varscandir,tumors)
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan1/3A.snp.Somatic.hc 1 2 3 4 T 3A.snp.Somatic.hc 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan1/11A.snp.Somatic.hc 1 2 3 4 T 11A.snp.Somatic.hc 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan1/13A.snp.Somatic.hc 1 2 3 4 T 13A.snp.Somatic.hc 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan1/15A.snp.Somatic.hc 1 2 3 4 T 15A.snp.Somatic.hc 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan1/17A.snp.Somatic.hc 1 2 3 4 T 17A.snp.Somatic.hc 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan1/25A.snp.Somatic.hc 1 2 3 4 T 25A.snp.Somatic.hc 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan1/33A.snp.Somatic.hc 1 2 3 4 T 33A.snp.Somatic.hc 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan1/37A.snp.Somatic.hc 1 2 3 4 T 37A.snp.Somatic.hc 0
#./count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan1/41A.snp.Somatic.hc 1 2 3 4 T 41A.snp.Somatic.hc 0

#1/23
#varscandir="/fh/scratch/delete30/dai_j/henan/varscan2"
#cd $varscandir 
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan2/3A.snp.Somatic.hc 1 2 3 4 T 3A.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan2/11A.snp.Somatic.hc 1 2 3 4 T 11A.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan2/13A.snp.Somatic.hc 1 2 3 4 T 13A.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan2/15A.snp.Somatic.hc 1 2 3 4 T 15A.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan2/17A.snp.Somatic.hc 1 2 3 4 T 17A.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan2/25A.snp.Somatic.hc 1 2 3 4 T 25A.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan2/29A.snp.Somatic.hc 1 2 3 4 T 29A.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan2/33A.snp.Somatic.hc 1 2 3 4 T 33A.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan2/37A.snp.Somatic.hc 1 2 3 4 T 37A.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/henan/varscan2/41A.snp.Somatic.hc 1 2 3 4 T 41A.snp.Somatic.hc 0

#varscandir="/fh/scratch/delete30/dai_j/escc/varscan2"
#cd $varscandir 
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T1.snp.Somatic.hc 1 2 3 4 T T1.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T2.snp.Somatic.hc 1 2 3 4 T T2.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T3.snp.Somatic.hc 1 2 3 4 T T3.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T4.snp.Somatic.hc 1 2 3 4 T T4.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T5.snp.Somatic.hc 1 2 3 4 T T5.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T6.snp.Somatic.hc 1 2 3 4 T T6.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T8.snp.Somatic.hc 1 2 3 4 T T8.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T9.snp.Somatic.hc 1 2 3 4 T T9.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T10.snp.Somatic.hc 1 2 3 4 T T10.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T11.snp.Somatic.hc 1 2 3 4 T T11.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T12.snp.Somatic.hc 1 2 3 4 T T12.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T13.snp.Somatic.hc 1 2 3 4 T T13.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T14.snp.Somatic.hc 1 2 3 4 T T14.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T15.snp.Somatic.hc 1 2 3 4 T T15.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T16.snp.Somatic.hc 1 2 3 4 T T16.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T17.snp.Somatic.hc 1 2 3 4 T T17.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/escc/varscan2/T18.snp.Somatic.hc 1 2 3 4 T T18.snp.Somatic.hc 0

#1/26
wgstumors=c("SRR1001842","SRR1002713","SRR999423","SRR1001466","SRR1002670","SRR1001823","SRR999489","SRR1002343","SRR1002722","SRR1002656",
            "SRR1002929","SRR999438","SRR1001915","SRR999594","SRR1001868","SRR1001635")
#varscandir="/fh/scratch/delete30/dai_j/varscan2"
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1001842.snp.Somatic.hc 1 2 3 4 T SRR1001842.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1002713.snp.Somatic.hc 1 2 3 4 T SRR1002713.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR999423.snp.Somatic.hc 1 2 3 4 T SRR999423.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1001466.snp.Somatic.hc 1 2 3 4 T SRR1001466.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1002670.snp.Somatic.hc 1 2 3 4 T SRR1002670.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1001823.snp.Somatic.hc 1 2 3 4 T SRR1001823.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR999489.snp.Somatic.hc 1 2 3 4 T SRR999489.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1002343.snp.Somatic.hc 1 2 3 4 T SRR1002343.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1002722.snp.Somatic.hc 1 2 3 4 T SRR1002722.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1002656.snp.Somatic.hc 1 2 3 4 T SRR1002656.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1002929.snp.Somatic.hc 1 2 3 4 T SRR1002929.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR999438.snp.Somatic.hc 1 2 3 4 T SRR999438.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1001915.snp.Somatic.hc 1 2 3 4 T SRR1001915.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR999594.snp.Somatic.hc 1 2 3 4 T SRR999594.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1001868.snp.Somatic.hc 1 2 3 4 T SRR1001868.snp.Somatic.hc 0
#/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/count6mutations.R /fh/scratch/delete30/dai_j/varscan2/SRR1001635.snp.Somatic.hc 1 2 3 4 T SRR1001635.snp.Somatic.hc 0














