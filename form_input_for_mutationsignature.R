#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
#variant_function 3,4,6,7,9;exonic_variant_function: 4,5,7,8,10
mutectdir=as.character(args[1]) #the current folder
colchr=as.integer(args[2])
colstart=as.integer(args[3])
colref=as.integer(args[4])
colalt=as.integer(args[5])
pref=as.character(args[6])

print(mutectdir)
print(pref)
mutationsignaturedir="/fh/fast/dai_j/CancerGenomics/Tools/mutationsignature"
output_mutationmatrix=paste0(mutectdir,"/",pref,"_mutationsignaturematrix.txt")
output_samplename=paste0(mutectdir,"/",pref,"_samplenames.txt")
output_inputmfile=paste0(mutationsignaturedir,"/input/",pref,"_inputformutationsignature.m")
output_inputmatfile=paste0(mutationsignaturedir,"/input/",pref,"_inputformutationsignature.mat")
print(output_mutationmatrix)
print(output_samplename)
print(output_inputmfile)
print(output_inputmatfile)

numsamples=(length(args)-6)/2
print(numsamples)
tumors=c()
for (i in 7:(6+numsamples))
{
  tumors[i-6]=as.character(args[i])
} 
maffiles=c()
for (i in (7+numsamples):length(args))
{
  maffiles[i-6-numsamples]=as.character(args[i])
  print(maffiles[i-6-numsamples])
}

countmutations=function(maffile,header=F,colchr,colstart,colref,colalt,samplename=NULL,pref=NULL,write2files=0)
{
  library(GenomicRanges)
  library(BSgenome)
  library("BSgenome.Hsapiens.UCSC.hg19")
  maftable=read.table(maffile,sep="\t",header=header)
  #samplename=as.character(maftable[1,colsamplename])
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
  results=list(convmat=res,covrate=resrate,lego=lego,transtable=transtable,samplename=samplename)
  
  if (write2files==1)
  {
    output1=paste0(pref,".transmat.txt")
    output2=paste0(pref,".transrate.txt")
    output3=paste0(pref,".lego.txt")
    output4=paste0(pref,".transtable.txt")
    write.table(results$convmat,file=output1,na="NA",sep="\t",row.names=T,col.names=T,quote=F)
    write.table(results$covrate,file=output2,na="NA",sep="\t",row.names=T,col.names=T,quote=F)
    write.table(results$lego,file=output3,na="NA",sep="\t",row.names=F,col.names=T,quote=F)
    write.table(results$transtable,file=output4,sep="\t",row.names=F,col.names=T,quote=F)
  }
  
  return(results)
  
}


types=read.table("/fh/fast/dai_j/CancerGenomics/Tools/mutationsignature/input/mutationsignature_types.txt")
subtypes=read.table("/fh/fast/dai_j/CancerGenomics/Tools/mutationsignature/input/mutationsignature_subtypes.txt")
uniqtypes=unique(types)
uniqtypesid=c(2,3,1,6,4,5) #block id used in count6mutations.R
trinuc=c("T","C","A","G")

formdata=function(tranmatfile)
{
  if (class(tranmatfile)=="character")
  {
    origdata=read.table(transmatfile,header=T)
  }
  if (class(tranmatfile)=="data.frame")
  {
    origdata=tranmatfile
  }
  
  res=data.frame(matrix(NA,ncol=1,nrow=nrow(types)))
  colnames(res)='mutation'
  for (i in 1:nrow(subtypes))
  {
    subtype=as.character(subtypes[i,1])
    pre=which(trinuc==substr(subtype,1,1))
    post=which(trinuc==substr(subtype,3,3))
    type=as.character(types[i,1])
    blockid=uniqtypesid[which(uniqtypes==type)]
    blocktable=origdata[((blockid-1)*4+1):(blockid*4),]
    res[i,1]=blocktable[pre,post]
    rownames(res)[i]=paste0(substr(subtype,1,1),'_',type,'_',substr(subtype,3,3))
  }
  return(res)
}

mutationmatrix=data.frame(matrix(NA,ncol=length(maffiles),nrow=nrow(types)))

samplenames=data.frame(matrix(NA,ncol=1,nrow=length(maffiles)))
samplenames[,1]=as.character(samplenames[,1])
for (i in 1:length(maffiles))
{
  resultlist=countmutations(maffiles[i],header=F,colchr,colstart,colref,colalt,samplename=tumors[i],pref,write2files=0)
  mutationmatrix[,i]=formdata(resultlist$convmat)
  samplenames[i,1]=resultlist$samplename
  print(samplenames[i,1])
}
#write to files
write.table(mutationmatrix,file=output_mutationmatrix,sep="\t",quote=F,row.names=F,col.names=F)
write.table(samplenames,file=output_samplename,quote=F,row.names=F,col.names=F)
filecon=file(output_inputmfile,"w")
print(output_inputmfile)
#load subtype and type
writeLines("load(\'21_WTSI_BRCA_whole_genome_substitutions.mat\')",con=filecon)
writeLines(paste0("cancerType=\'",pref,"\'"),con=filecon)
writeLines(paste0("originalGenomes=textread(\'",output_mutationmatrix,"\')"),con=filecon)
writeLines(paste0("sampleNames=textread(\'",output_samplename,"\',\'%s\')"),con=filecon)
writeLines(paste0("save(\'",output_inputmatfile,"\')"),con=filecon)
close(filecon)

