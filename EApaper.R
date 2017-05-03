#!/usr/bin/env Rscript
library(GenomicRanges)
grp1normals=c('SRR1002719','SRR999433','SRR999687','SRR1000378','SRR1002786','SRR999559','SRR1001730','SRR10018461','SRR1002703','SRR999599','SRR1002792','SRR1001839','SRR9994281','SRR1002710','SRR9995631','SRR1021476')
grp1tumors=c('SRR1001842','SRR1002713','SRR999423','SRR1001466','SRR1002670','SRR1001823','SRR999489','SRR1002343','SRR1002722','SRR1002656','SRR1002929','SRR999438','SRR1001915','SRR999594','SRR1001868','SRR1001635')

grp2normals=paste0((c(3,11,13,15,17,25,29,33,37,41)+1),"A")
grp2tumors=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")

grp1dir="/fh/scratch/delete30/dai_j/mutect4"
grp2dir="/fh/scratch/delete30/dai_j/henan/mutect5"

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

#mutation frequency
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

#mutation subtype proportion
dulakstable6=read.table(file="/fh/fast/dai_j/CancerGenomics/EAC/Dulak_Stable6.txt",header=T)
srridx1=rep(0,length(grp1tumors))
for (i in 1:length(grp1tumors))
{
  tmp=which(esotranstable[,1]==grp1tumors[i])
  tmp1=esotranstable[tmp,2]
  srridx1[i]=which(dulakstable6[,1]==tmp1)
}
#transfer into grp1's order
srrdulakstable6=dulakstable6[srridx1,]


#--------------------------

#mutations-----------------------------

#compute mutation proportions
countmutations_fromtable=function(maftable,colchr,colstart,colref,colalt,colsamplename=NULL,pref=NULL,write2files=0)
{
  library(GenomicRanges)
  library(BSgenome)
  library("BSgenome.Hsapiens.UCSC.hg19")
  if (! is.null(colsamplename))
  {
    samplename=as.character(maftable[1,colsamplename])
  }else
  {
    samplename="NA"
  }
  
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
  transcount=resrate=data.frame(matrix(0,ncol=1,nrow=7))
  resrate[1,]=c2t_rate
  resrate[2,]=c2a_rate
  resrate[3,]=c2g_rate
  resrate[4,]=a2g_rate
  resrate[5,]=a2c_rate
  resrate[6,]=a2t_rate
  resrate[7,]=aa2c_rate
  transcount[1,]=c2t
  transcount[2,]=c2a
  transcount[3,]=c2g
  transcount[4,]=a2g
  transcount[5,]=a2c
  transcount[6,]=a2t
  transcount[7,]=aa2c
  rownames(transcount)=rownames(resrate)=c('c2t','c2a','c2g','a2g','a2c','a2t','aa2c')
  colnames(resrate)=c("rate")
  results=list(convmat=res,transrate=resrate,transcount=transcount,lego=lego,transtable=transtable,samplename=samplename)
  
  if (write2files==1)
  {
    output1=paste0(pref,".transmat.txt")
    output2=paste0(pref,".transrate.txt")
    output3=paste0(pref,".lego.txt")
    output4=paste0(pref,".transtable.txt")
    write.table(results$convmat,file=output1,na="NA",sep="\t",row.names=T,col.names=T,quote=F)
    write.table(results$transrate,file=output2,na="NA",sep="\t",row.names=T,col.names=T,quote=F)
    write.table(results$lego,file=output3,na="NA",sep="\t",row.names=F,col.names=T,quote=F)
    write.table(results$transtable,file=output4,sep="\t",row.names=F,col.names=T,quote=F)
  }
  
  return(results)
}

computeall=function(grpdir,grptumors)
{
  transcount=transratemat=data.frame(matrix(NA,nrow=7,ncol=length(grptumors)))
  rownames(transcount)=rownames(transratemat)=c('c2t','c2a','c2g','a2g','a2c','a2t','aa2c')
  colnames(transcount)=colnames(transratemat)=grptumors
  legomat=data.frame(matrix(NA,nrow=length(grptumors),ncol=96))
  rownames(legomat)=grptumors
  nummutation=data.frame(matrix(0,nrow=length(grptumors),ncol=1))
  rownames(nummutation)=grptumors
  
  for (i in 1:length(grptumors))
  {
    print(grptumors[i])
    maffile=paste0(grpdir,"/",grptumors[i],".mutectmaf.ancotator")
    if ( ! file.exists(maffile))
    {
      stop(paste0("run runmutect_mutsig_usefull192coverage.sh first!"))
    }
    maftable=read.table(maffile,header=T,sep="\t",quote="")
    maftable$Variant_Classification=as.character(maftable$Variant_Classification)
    colchr=which(colnames(maftable)=="Chromosome")
    colstart=which(colnames(maftable)=="Start_position")
    colref=which(colnames(maftable)=="Reference_Allele")
    colalt=which(colnames(maftable)=="Tumor_Seq_Allele2")
    result=countmutations_fromtable(maftable,colchr,colstart,colref,colalt)
    transratemat[,i]=result$transrate
    transcount[,i]=result$transcount
    legomat[i,]=result$lego
    colnames(legomat)=names(result$lego)
    nummutation[i,1]=nrow(maftable)
  }
  
  #used for lego plot
  alllegofile=paste0(grpdir,"/all_lego.txt")
  alllego=data.frame(matrix(0,nrow=1,ncol=ncol(legomat)))
  colnames(alllego)=colnames(legomat)
  for (i in 1:ncol(legomat))
  {
    alllego[1,i]=sum(legomat[,i])
  }
  write.table(alllego,file=alllegofile,row.names = F,sep="\t",quote=F)
  #write lego files
  for (i in 1:length(grptumors))
  {
    file1=paste0(grpdir,"/",grptumors[i],".lego.txt")
    write.table(legomat[i,],file=file1,row.names = F,sep="\t",quote=F)
  }
  
  #form mutationmatrix for signaure
  mutationmatrixfile=paste0(grpdir,"/mutationsignaturematrix.txt")
  mutationmatrix=data.frame(matrix(0,nrow=96,ncol=length(grptumors)))
  subtypestable=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/mutationsignature/mutationsubtypes.txt",stringsAsFactors = F)
  for (i in 1:nrow(subtypestable))
  {
    tmp=unlist(strsplit(subtypestable[i,1],"[",fixed = T))
    b=tmp[1]
    tmp1=unlist(strsplit(tmp[2],"]",fixed=T))
    m=tmp1[1]
    e=tmp1[2]
    if (m=="C>A") legotxt=paste0("C.A.G.T_",b,"x",e)
    if (m=="C>G") legotxt=paste0("C.G.G.C_",b,"x",e)
    if (m=="C>T") legotxt=paste0("C.T.G.A_",b,"x",e)
    if (m=="T>A") legotxt=paste0("T.A.A.T_",b,"x",e)
    if (m=="T>C") legotxt=paste0("T.C.A.G_",b,"x",e)
    if (m=="T>G") legotxt=paste0("T.G.A.C_",b,"x",e)
    idxrow=which(colnames(legomat)==legotxt)
    for (j in 1:length(grptumors))
    {
      mutationmatrix[i,j]=legomat[j,idxrow]
    }
  }
  write.table(mutationmatrix,file=mutationmatrixfile,row.names = F,col.names = F,sep="\t",quote=F)
  
  results=list(transratemat=transratemat,legomat=legomat,nummutation=nummutation)
}

dulak_res=computeall(grp1dir,grp1tumors)
henan_res=computeall(grp2dir,grp2tumors)
#save(dulak_res,henan_res,file="compare_globals_2grps_mutect4.RData")

#load results
load("compare_globals_2grps_mutect4.RData")
nsvtypes=c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                 "Splice_Site","Translation_Start_Site","Nonstop_Mutation","Targeted_Region","De_novo_Start_InFrame","De_novo_Start_OutOfFrame")

#stratify mutations
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
  mutationrate=data.frame(matrix(NA,nrow=length(grptumors),ncol=6))
  colnames(mutationrate)=c("total","total1","inter","intro","exon","nonsilent")
  len=data.frame(matrix(NA,nrow=length(grptumors),ncol=5))
  colnames(len)=c("total","inter","intro","exon","nonsilent")
  
  for (i in 1:length(grptumors))
  {
    cat(i,"..")
    wigfile=paste0(grpdir,"/",grptumors[i],".wig.bed.txt")
    if ( ! file.exists(wigfile)) stop("run run_wig2bed.sh!")
    wigtable=read.table(wigfile,sep="\t")
    wigtable=wigtable[wigtable$V1 %in% chrs,]
    gr_wigtable=GRanges(seqnames = wigtable$V1,ranges = IRanges(start=wigtable$V2,end=wigtable$V3))
    len$total[i]=sum(as.numeric(width(gr_wigtable)))
    mutationrate$total[i]=grpmuations$total[i]/len$total[i]*10^6
    mutationrate$total1[i]=srrdulakstable2$GenomicMutations[i]/len$total[i]*10^6 #use the total count from dulak
    
    tmp=intersect(gr_intron1,gr_wigtable)
    tmp2=intersect(gr_intron_exon1,gr_wigtable)
    len$intro[i]=sum(as.numeric(width(tmp)))-sum(as.numeric(width(tmp2)))
    len$intro[i]=len$intro[i]*0.9
    
    tmp1=intersect(gr_coding1,gr_wigtable)
    len$exon[i]=sum(width(tmp1))*1.1
    
    tmp3=intersect(gr_wholegene1,gr_wigtable)
    len$inter[i]=len$total[i]-sum(as.numeric(width(tmp3)))
    len$inter[i]=len$inter[i]*0.9
    
    tmp4=intersect(gr_coding1,gr_wigtable)
    len$nonsilent[i]=sum(as.numeric(width(tmp4)))
    
    mutationrate$inter[i]=grpmuations$inter[i]/len$inter[i]*10^6
    mutationrate$intro[i]=grpmuations$intro[i]/len$intro[i]*10^6
    mutationrate$exon[i]=grpmuations$exon[i]/len$exon[i]*10^6
    mutationrate$nonsilent[i]=grpmuations$nonsilent[i]/len$nonsilent[i]*10^6
    #
  }
  return(result=list(len=len,mutationrate=mutationrate))
}


grp1mutationrate=countmutationrate(grp1mutations,grp1tumors,grpdir=grp1dir)
grp2mutationrate=countmutationrate(grp2mutations,grp2tumors,grpdir=grp2dir)

#save(grp1mutations,grp2mutations,grp1mutationrate,grp2mutationrate,file="mutationrate.RData")
load("mutationrate.RData")
median(grp1mutations$total)
range(grp1mutations$total)
median(grp2mutations$total)
#[1] 7047.5
range(grp2mutations$total)
#[1]  2070 53570
median(grp1mutationrate$mutationrate$total)
range(grp1mutationrate$mutationrate$total)
median(grp2mutationrate$mutationrate$total)
#[1] 2.562489
range(grp2mutationrate$mutationrate$total)
#[1]  0.7535181 19.6544562

median(grp1mutationrate$mutationrate$total)
range(grp1mutationrate$mutationrate$total)
median(grp1mutationrate$mutationrate$inter)
median(grp1mutationrate$mutationrate$intro)
median(grp1mutationrate$mutationrate$exon)

median(grp2mutationrate$mutationrate$total)
#[1] 2.562489
range(grp2mutationrate$mutationrate$total)
#[1]  0.7535181 19.6544562
median(grp2mutationrate$mutationrate$inter)
#[1] 3.446513
range(grp2mutationrate$mutationrate$inter)
median(grp2mutationrate$mutationrate$intro)
#[1] 2.122243
range(grp2mutationrate$mutationrate$intro)
median(grp2mutationrate$mutationrate$exon)
#[1] 2.081696
range(grp2mutationrate$mutationrate$nonsilent)
median(grp2mutations$nonsilent/grp2mutationrate$len$nonsilent)*10^6
#[1] 0.9159641
range((grp2mutations$nonsilent/grp2mutationrate$len$nonsilent)*10^6)
#[1]  0.09154503 13.55586453

median(grp1mutations$nonsilent/grp1mutationrate$len$nonsilent)*10^6
#[1] 4.184271
range((grp1mutations$nonsilent/grp1mutationrate$len$nonsilent)*10^6)
#[1] 2.953471 8.717601

computep=function(x,y)
{
  
  print(paste0("meanx:",mean(x)))
  print(paste0("meany:",mean(y)))
  print(paste0("medianx:",median(x)))
  print(paste0("mediany:",median(y)))
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
  print(paste0("meanx:",mean(x)))
  print(paste0("meany:",mean(y)))
  print(paste0("medianx:",median(x)))
  print(paste0("mediany:",median(y)))
  print(paste0("rangex:",range(x)))
  print(paste0("rangey:",range(y)))
  p=wilcox.test(x,y)$p.value
  res=list(min1=min(x,na.rm=T),max1=max(x,na.rm=T),min2=min(y,na.rm=T),max2=max(y,na.rm=T),
           mean1=mean(x),mean2=mean(y),median1=median(x),median2=median(y),pvalue=p)
}
plot2groups=function(item,datatable,main="",xlabel="",ylabel="",ymax=NULL,yat=NULL,adjustymin=T,colorgolden4=F,usewilcox=T)
{
  if (usewilcox==F)
  {
    res1=computep(datatable[1:length(grp1tumors),item],
                  datatable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),item])
  }else
  {
    res1=computep_wilcox(datatable[1:length(grp1tumors),item],
                         datatable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),item])
  }
  
  pvalue1=format(res1$pvalue,digits=2)
  idx=which(colnames(datatable)==item)
  if (is.null(ymax))
  {
    ymax=max(datatable[,idx]) #mannualy set ymax
    if (colorgolden4==T) #keep space for legend
      ymax=1.2*ymax
  }
  ymin=min(datatable[,idx])
  if (adjustymin==T)
  {
    ymin1=ymin-0.15*ymax
    if (ymin1<0) ymin1=0
  }else
  {
    ymin1=0
  }
  ymax1=max(datatable[datatable$sample %in% grp1tumors,idx])
  ymax2=max(datatable[datatable$sample %in% grp2tumors,idx])
  # tmp=boxplot(as.formula(paste0(item," ~ dataset")), data=datatable, 
  #         outline = FALSE,     ## avoid double-plotting outliers, if any
  #         cex.axis=1.3,
  #         cex.main=1.3,
  #         ylim=c(ymin1,ymax*1.25),
  #         plot=F,
  #         main=main)
  if (colorgolden4==T & "col" %in% colnames(datatable))
  {
    beeswarm(as.formula(paste0(item," ~ dataset")), data = datatable, 
             pwcol = datatable$col, pch = 19,add = TRUE, cex=1.2)
    legend("topright",legend=c("US-EAC","CH-EAC","CH-EAC-Golden4"),col=c("green","blue","firebrick"),pch=16,horiz=T,text.width=0.7,x.intersp=0.1,xjust=0.5,yjust=0)
  }else
  {
    if (is.null(yat))
    {
      beeswarm(as.formula(paste0(item," ~ dataset")), data = datatable, 
               col = color, pch = 19,add = FALSE, cex=1.2,ylim=c(ymin1,ymax*1.25),main=main,
               xlab=xlabel,ylab=ylabel,xlim=c(0.45,2.55),cex.axis=1.2,cex.main=1.2,cex.lab=1.2)
    }else
    {
      beeswarm(as.formula(paste0(item," ~ dataset")), data = datatable, yaxt="n",
               col = color, pch = 19,add = FALSE, cex=1.5,ylim=c(ymin1,ymax*1.25),main=main,
               xlab=xlabel,ylab=ylabel,xlim=c(0.45,2.55),cex.axis=1.3,cex.main=1.3,cex.lab=1.3)
      axis(2,at=yat,labels=as.character(yat),cex.axis=1.2)
    }
  }
  ymax12=max(ymax1,ymax2)
  segments(1,ymax12*1.05+ymax*0.05,2,ymax12*1.05+ymax*0.05,lwd=2)
  segments(1,ymax12*1.05+ymax*0.05,1,ymax12*1.05,lwd=2)
  segments(2,ymax12*1.05+ymax*0.05,2,ymax12*1.05,lwd=2)
  text(1.5,ymax12*1.05+ymax*0.13,paste0("p=",pvalue1),cex=1.2)
  median1=median(datatable[datatable$sample %in% grp1tumors,idx])
  median2=median(datatable[datatable$sample %in% grp2tumors,idx])
  segments(0.55,median1,1.45,median1,col=color[1],lwd=2,lty=3)
  segments(1.55,median2,2.45,median2,col=color[2],lwd=2,lty=3)
}


#draw boxplot
library("beeswarm")
restable=data.frame(matrix(NA,nrow=length(grp1tumors)+length(grp2tumors),ncol=2))
colnames(restable)=c("dataset","sample")
restable[1:length(grp1tumors),1]="US-EAC"
restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),1]="CH-EAC"
restable[,2]=c(grp1tumors,grp2tumors)
#reorder group
restable$dataset=factor(restable$dataset,c("US-EAC","CH-EAC"))
color=c("green","blue")
#total mutation frequency
restable=cbind(restable,totalrate=c(grp1mutationrate$mutationrate$total,grp2mutationrate$mutationrate$total))
#figure1
#par(mar=c(5.1,4.1,4.1,2.1))
par(mar=c(3.5,5.1,3.1,2.1))
par(mfrow=c(1,1))
outfig="EA_muationfrequency.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
plot2groups(item="totalrate",datatable=restable,adjustymin=T,main='',ylab="Mutations/Mb",usewilcox=T)
dev.off()

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
median(grp2a2c[,2])
#[1] 0.07514751
range((grp2a2c[,2]))
median(grp1a2c[,2])
#[1] 0.3030532

#aa2c
grp1aa2ac=counttransrate(dulak_res$transratemat,grp1tumors,opt="aa2c")
grp2aa2ac=counttransrate(henan_res$transratemat,grp2tumors,opt="aa2c")

restable=cbind(restable,c(grp1aa2ac[,2],grp2aa2ac[,2]))
colnames(restable)[ncol(restable)]="aa2ac"
#aa2ac/a2c
grp1aa2acprop=dulak_res$transratemat["aa2c",]/dulak_res$transratemat["a2c",]
range(grp1aa2acprop)
median(as.numeric(grp1aa2acprop[1,]))
#[1] 0.7044438
grp2aa2acprop=henan_res$transratemat["aa2c",]/henan_res$transratemat["a2c",]
range(grp2aa2acprop)
median(as.numeric(grp2aa2acprop[1,]))
#[1] 0.3115777
median(as.numeric(henan_res$transratemat["aa2c",]))
#[1] 0.02400435
median(as.numeric(dulak_res$transratemat["aa2c",]))
#[1] 0.2221981
tmp=computep_wilcox(as.numeric(dulak_res$transratemat["aa2c",]),as.numeric(henan_res$transratemat["aa2c",]))
tmp$pvalue
#[1] 3.765248e-07
#a2g
grp1a2g=counttransrate(dulak_res$transratemat,grp1tumors,opt="a2g")
grp2a2g=counttransrate(henan_res$transratemat,grp2tumors,opt="a2g")
restable=cbind(restable,c(grp1a2g[,2],grp2a2g[,2]))
colnames(restable)[ncol(restable)]="a2g"
#a2t
grp1a2t=counttransrate(dulak_res$transratemat,grp1tumors,opt="a2t")
grp2a2t=counttransrate(henan_res$transratemat,grp2tumors,opt="a2t")
restable=cbind(restable,c(grp1a2t[,2],grp2a2t[,2]))
colnames(restable)[ncol(restable)]="a2t"
median(grp2c2t[,2])
range(grp2c2t[,2])
#c2a
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
  text(mean(c(x1,x2)),ymax12*1.05+ymax*0.1,paste0("p=",pvalue),cex=1.2)
}

plotbeeswarm=function(mat,formula,colors=c("skyblue", "black", "red", "gray", "forestgreen", "tan3"),xlabel,ylabel,labels)
{
  ymax=max(mat$value,na.rm=T)
  numgroups=length(unique(mat[,2]))
  colors1=rep(colors,each=2)
  ats=NULL
  for (i in 1:numgroups)
  {
    at1=i-0.2
    at2=i+0.2
    ats=c(ats,at1,at2)
  }
  beeswarm(as.formula(formula), data = mat, xaxt='n',
                     col = colors1, pch = 19,add = FALSE, cex=1.2,ylim=c(0,ymax*1.25),main="",
                     xlab=xlabel,ylab=ylabel,cex.axis=1.2,cex.main=1.2,cex.lab=1.2,at=ats)
  axis(side=1, at=seq(1,numgroups,1), labels=labels,
       cex.axis=1.2,cex.lab=1.2,line=F,lwd.ticks=1.2)
  #add median
  for (i in 1:numgroups)
  {
    groupdata1=mat[mat[,1]==unique(mat[,1])[1] & mat[,2]==unique(mat[,2])[i] ,3]
    groupdata2=mat[mat[,1]==unique(mat[,1])[2] & mat[,2]==unique(mat[,2])[i] ,3]
    median1=median(groupdata1,na.rm=T)
    median2=median(groupdata2,na.rm=T)
    segments(i-0.38,median1,i-0.02,median1,col=colors1[(i-1)*2+1],lwd=3,lty=3)
    segments(i+0.02,median2,i+0.38,median2,col=colors1[(i-1)*2+2],lwd=3,lty=3)
    addpvlue(groupdata1,groupdata2,ymax,i-0.2,i+0.2)
  }
}



ymax=max(c(grp1c2a[,2],grp2c2a[,2],grp1c2g[,2],grp2c2g[,2],grp1c2t[,2],grp2c2t[,2],
           grp1a2g[,2],grp2a2g[,2],grp1a2c[,2],grp2a2c[,2],grp1a2t[,2],grp2a2t[,2]),na.rm=T)

#plot figure1(a)
plotfig1a=function()
{
  par(par(mar=c(4.1,5.1,2.1,2.1)))
  par(mfrow=c(1,1))
  plot2groups(item="totalrate",datatable=restable,adjustymin=T,main='',ylab="Mutations/Mb",usewilcox=T)
  #mtext("(a)",side=1,line=3.5,cex=1.6)
}
postscript(file="fig1a.ps",width = 6,height = 6,paper= "special",horizontal=F)
plotfig1a()
dev.off()

plotfig1b=function()
{
  #draw boxplot
  #form matrix
  mat=data.frame(matrix(NA,nrow=12,ncol=length(grp1tumors)))
  mat[1,]=grp1c2a[,2]
  mat[2,1:length(grp2tumors)]=grp2c2a[,2]
  mat[3,]=grp1c2g[,2]
  mat[4,1:length(grp2tumors)]=grp2c2g[,2]
  mat[5,]=grp1c2t[,2]
  mat[6,1:length(grp2tumors)]=grp2c2t[,2]
  mat[7,]=grp1a2g[,2]
  mat[8,1:length(grp2tumors)]=grp2a2g[,2]
  mat[9,]=grp1a2c[,2]
  mat[10,1:length(grp2tumors)]=grp2a2c[,2]
  mat[11,]=grp1a2t[,2]
  mat[12,1:length(grp2tumors)]=grp2a2t[,2]
  mutationtypes=c("C>A","C>G","C>T","A>G","A>C","A>T")
  mutationtypes=rep(mutationtypes,each=2)
  locations=rep(c("US","CH"),6)
  mat$mutationtypes=mutationtypes
  mat$locations=locations
  mat$locations=factor(mat$locations,c("US","CH"))
  mat$mutationtypes=factor(mat$mutationtypes,c("C>A","C>G","C>T","A>G","A>C","A>T"))
  library(reshape2)
  mat1=melt(mat, id = c('locations', 'mutationtypes'))
  mat1=mat1[,-3]
  mat1$value=as.numeric(mat1$value)
  boxplots=boxplot(value~locations + mutationtypes,data = mat1,boxwex=0.2,
                   at = c(0.5,0.8,1.5,1.8,2.5,2.8,3.5,3.8,4.5,4.8,5.5,5.8),
                   xaxt='n',ylim=c(0,1.3*max(mat1$value,na.rm=T)),
                   xlim=c(0,6),
                   border=rep(c("skyblue", "black", "red", "gray", "forestgreen", "tan3"),each=2),
                   cex=1,cex.lab=1.2,cex.axis=1.2,boxlwd=2,
                   frame=F,
                   xlab="Type of mutation",ylab="Fraction of mutation"
  )
  
  axis(side=1, at=c(0.65,1.65,2.65,3.65,4.65,5.65), labels=c("C>A","C>G","C>T","A>G","A>C","A>T"),
       cex.axis=1.2,cex.lab=1.2,line=F,lwd.ticks=1.2)
  segments(0,par("usr")[3],6,par("usr")[3],lwd=1.2)
  
  pc2a=computep_wilcox(grp1c2a[,2],grp2c2a[,2])$pvalue
  if (pc2a<0.05)
  {
    addpvlue(grp1c2a[,2],grp2c2a[,2],ymax,0.5,0.8)
  }
  pc2g=computep_wilcox(grp1c2g[,2],grp2c2g[,2])$pvalue
  if (pc2g<0.05)
  {
    addpvlue(grp1c2g[,2],grp2c2g[,2],ymax,1.5,1.8)
  }
  pc2t=computep_wilcox(grp1c2t[,2],grp2c2t[,2])$pvalue
  if (pc2t<0.05)
  {
    addpvlue(grp1c2t[,2],grp2c2t[,2],ymax,2.5,2.8)
  }
  pa2g=computep_wilcox(grp1a2g[,2],grp2a2g[,2])$pvalue
  if (pa2g<0.05)
  {
    addpvlue(grp1a2g[,2],grp2a2g[,2],ymax,3.5,3.8)
  }
  pa2c=computep_wilcox(grp1a2c[,2],grp2a2c[,2])$pvalue
  if (pa2c<0.05)
  {
    addpvlue(grp1a2c[,2],grp2a2c[,2],ymax,4.5,4.8)
  }
  pa2t=computep_wilcox(grp1a2t[,2],grp2a2t[,2])$pvalue
  if (pa2t<0.05)
  {
    addpvlue(grp1a2t[,2],grp2a2t[,2],ymax,5.5,5.8)
  }
  mtext("(b)",side=1,line=4.5,cex=1.6)
}

plotfig1c=function()
{
  # library(png)
  # img1c = readPNG("figure1c.png")
  # plot.new()
  # usr<-par("usr")
  # rasterImage(img1c, usr[1], usr[3], usr[2], usr[4],interpolate = F)
  library("imager")
  img1c=load.image("figure1c.png")
  plot(img1c,axes=F)
  mtext("(c)",side=1,line=1,cex=1.6)
}

#plot figure1
outfig="EA_figure1.png"
png(outfig, width = 8, height = 8, units = 'in', res=300)
layout(mat=matrix(c(1,2,3,3), ncol=2, byrow=TRUE),heights=c(1,1.5))
par(mar=c(5.5,5.1,3.1,2.1))
plotfig1a()
plotfig1b()
plotfig1c()
dev.off()



#draw plot figure2
outfig="EA_muationproportion.png"

png(outfig, width = 6, height = 9, units = 'in', res=300)
par(mar=c(4.5,5.1,3.1,2.1))
par(mfrow=c(3,2))
plot2groups(item="a2c",datatable=restable,main='A>C',ylabel="Proportion",xlabel="(a)",yat=c(0.1,0.2,0.3,0.4,0.5,0.6))
plot2groups(item="a2g",datatable=restable,main='A>G',ylabel="Proportion",xlabel="(b)",yat=c(0.15,0.2,0.25,0.3))
plot2groups(item="a2t",datatable=restable,main='A>T',adjustymin=T,ylabel="Proportion",xlabel="(c)",yat=c(0.05,0.1,0.15,0.2))
plot2groups(item="c2a",datatable=restable,main='C>A',ylabel="Proportion",xlabel="(d)",yat=c(0.05,0.1,0.15,0.2,0.25,0.3))
plot2groups(item="c2g",datatable=restable,main='C>G',ylabel="Proportion",xlabel="(e)",yat=c(0.05,0.1,0.15,0.2))
plot2groups(item="c2t",datatable=restable,main='C>T',ylabel="Proportion",xlabel="(f)",yat=c(0.1,0.2,0.3,0.4,0.5))
dev.off()

draw2barplot=function(data1,data2,ylable="",legends=c("dulak","replicate"),main="")
{
  tmp=data.frame(data1=data1,data2=data2)
  tmp=t(tmp)
  ymax=max(c(data1,data2))
  barplot(as.matrix(tmp),beside=TRUE,col=c("blue","green"),ylab=ylable,ylim=c(0,ymax*1.3),main=main)
  legend("topright",legend=legends,fill=c("blue","green"))
  computep(data1,data2)
}
par(mfrow=c(1,1))
par(mar=c(5.1,5.1,3.1,2.1))
draw2barplot(srrdulakstable2$MutationsperMb,grp1mutationrate$mutationrate$total,ylable ="total mutation rate" )
draw2barplot(srrdulakstable3$IntronMutationRate,grp1mutationrate$mutationrate$intro,ylable ="intron mutation rate" )
draw2barplot(srrdulakstable3$IGRMutationRate,grp1mutationrate$mutationrate$inter,ylable ="intergenic mutation rate" )
draw2barplot(srrdulakstable3$ExonMutationRate,grp1mutationrate$mutationrate$exon,ylable ="exon mutation rate" )
t.test(srrdulakstable2$MutationsperMb,grp1mutationrate$mutationrate$total)$p.value
t.test(srrdulakstable3$IntronMutationRate,grp1mutationrate$mutationrate$intro)$p.value
t.test(srrdulakstable3$IGRMutationRate,grp1mutationrate$mutationrate$inter)$p.value
t.test(srrdulakstable3$ExonMutationRate,grp1mutationrate$mutationrate$exon)$p.value

#mutation heatmap
library("ComplexHeatmap")
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
mutectdir="/fh/scratch/delete30/dai_j/henan/mutect3"
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
outfig="EA_muationheatmap.png"
png(outfig, width = 6, height = 8, units = 'in', res=300)
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

name="US-EAC"
mutectdir="/fh/scratch/delete30/dai_j/mutect4"
segfiles1=paste0(mutectdir,"/",grp1tumors,".mutectmaf.ancotator")
data=formdata(segfiles1)
nrow(data)
#[1] 2257
numsnv=rep(0,length(grp1tumors))
for (i in 1:length(grp1tumors))
{
  tmptable=data[data$Tumor_Sample_Barcode==grp1tumors[i],]
  numsnv[i]=nrow(tmptable)
}

mat=readmutation1(data,grp1tumors)
num_mutations=sapply(1:nrow(mat),function(i){
  sum(mat[i,] != " ")
})
#get 20% cutoff
num_row=sum(num_mutations/ncol(mat)>=0.2)
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
#on average coverage
sum(grp2cnvlen[,2:4])/length(grp2tumors)/3200
#gain
median(grp1cnvlen[,2])
median(grp2cnvlen[,2])
#loss
median(grp1cnvlen[,3])
median(grp2cnvlen[,3])
#loh
median(grp1cnvlen[,4])
median(grp2cnvlen[,4])
#Figure4
outfig="EA_CNAlength.png"
png(outfig, width = 9, height = 4, units = 'in', res=300)
par(mfrow=c(1,3))
plot2groups(item="lengain",datatable=restable,main='Amplifications',xlabel="(a)",ylabel="Length of CNA (Mb)")
plot2groups(item="lenloss",datatable=restable,main='Deletions',xlabel="(b)",ylabel="Length of CNA (Mb)")
plot2groups(item="lenloh",datatable=restable,main='LOH',xlabel="(c)",ylabel="Length of CNA (Mb)")
dev.off()

#figure4(a)
plotfig4a=function()
{
  #par(mar=c(6.5,5.1,1.1,2.1))
  mat=data.frame(matrix(NA,nrow=6,ncol=length(grp1tumors)))
  mat[1,]=grp1cnvlen$lengain
  mat[2,1:length(grp2tumors)]=grp2cnvlen$lengain
  mat[3,]=grp1cnvlen$lenloss
  mat[4,1:length(grp2tumors)]=grp2cnvlen$lenloss
  mat[5,]=grp1cnvlen$lenloh
  mat[6,1:length(grp2tumors)]=grp2cnvlen$lenloh
  cntypes=c("gain","loss","loh")
  cntypes=rep(cntypes,each=2)
  locations=rep(c("US","CH"),3)
  mat$cntypes=cntypes
  mat$locations=locations
  mat$locations=factor(mat$locations,c("US","CH"))
  mat$cntypes=factor(mat$cntypes,c("gain","loss","loh"))
  library(reshape2)
  mat1=melt(mat, id = c('locations', 'cntypes'))
  mat1=mat1[,-3]
  mat1$value=as.numeric(mat1$value)
  boxplots=boxplot(value~locations + cntypes,data = mat1,boxwex=0.2,
                   at = c(0.5,0.8,1.5,1.8,2.5,2.8),
                   xaxt='n',ylim=c(0,1.3*max(mat1$value,na.rm=T)),
                   xlim=c(0,3),
                   border=rep(c("red", "black", "skyblue"),each=2),
                   cex=1,cex.lab=1.2,cex.axis=1.2,boxlwd=2,
                   frame=F,
                   xlab="Type of CNA",ylab="Length of CNA (Mb)"
  )
  
  axis(side=1, at=c(0.65,1.65,2.65), labels=c("Gain","Loss","LOH"),
       cex.axis=1.2,cex.lab=1.2,line=F,lwd.ticks=1.2)
  segments(0,par("usr")[3],3,par("usr")[3],lwd=1.2)
  ymax=max(mat1$value,na.rm=T)
  pgain=computep_wilcox(grp1cnvlen$lengain,grp2cnvlen$lengain)$pvalue
  addpvlue(grp1cnvlen$lengain,grp2cnvlen$lengain,ymax,0.5,0.8)
  ploss=computep_wilcox(grp1cnvlen$lenloss,grp2cnvlen$lenloss)$pvalue
  addpvlue(grp1cnvlen$lenloss,grp2cnvlen$lenloss,ymax,1.5,1.8)
  ploh=computep_wilcox(grp1cnvlen$lenloh,grp2cnvlen$lenloh)$pvalue
  addpvlue(grp1cnvlen$lenloh,grp2cnvlen$lenloh,ymax,2.5,2.8)
  mtext("(a)",side=1,line=4.5,cex=1.6)
}

#plot figure4(b)
plotfig4b=function()
{
  library("imager")
  if (!file.exists("fig4b.png"))
  {
    gisticdir2="/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel1_brlen0.98_broad1"
    img4b1=load.image(paste0(gisticdir2,"/del_qplot.png"))
    img4b2=load.image(paste0(gisticdir2,"/amp_qplot.png"))
    png(filename = "fig4b.png",res=300,width = 6,height = 3,units ="in")
    par(mar=c(0.3,0.5,0.2,0.3))
    par(mfrow=c(1,2))
    plot(img4b1,axes=F)
    plot(img4b2,axes=F)
    dev.off()
  }
  
  img4b=load.image("fig4b.png")
  par(mar=c(0.3,0.5,0.2,0.3))
  plot(img4b,axes=F)
  mtext("(b)",side=1,line=4.5,cex=1.6)
}

#gistic result
library(GenomicRanges)
library("biomaRt")
library(RCircos)
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

countoverlap2segs=function(segdata1,segdata2)
{
  cytobands=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt',header=T)
  gr_cytobands=GRanges(seqnames=cytobands$chrom,ranges=IRanges(start=cytobands$start,end=cytobands$end),name=cytobands$name)
  segdata=rbind.data.frame(segdata1[,1:4],segdata2[,1:4])
  chrs=paste0("chr",c(1:22,"X","Y"))
  segdata_=NULL
  for (i in 1:length(chrs))
  {
    tmp=segdata[segdata$chr==chrs[i],]
    tmp=tmp[order(tmp$start),]
    segdata_=rbind(segdata_,tmp)
  }
  gr_segdata1=GRanges(seqnames=segdata1$chr,ranges=IRanges(start=segdata1$start,end=segdata1$end),cytoband=segdata1$cytoband)
  gr_segdata1=reduce(gr_segdata1)
  gr_segdata2=GRanges(seqnames=segdata2$chr,ranges=IRanges(start=segdata2$start,end=segdata2$end),cytoband=segdata2$cytoband)
  gr_segdata2=reduce(gr_segdata2)
  gr_segdata=GRanges(seqnames=segdata_$chr,ranges=IRanges(start=segdata_$start,end=segdata_$end),cytoband=segdata_$cytoband)
  gr_segdata=reduce(gr_segdata)
  overlapseg=data.frame(matrix(NA,nrow=0,ncol=4))
  colnames(overlapseg)=c("cytoband","chr","start","end")
  overlapseg$cytoband=as.character(overlapseg$cytoband)
  overlapseg$chr=as.character(overlapseg$chr)
  for (i in 1:length(gr_segdata1))
  {
    olap=intersect(gr_segdata1[i,],gr_segdata2)
    if (length(olap)>0)
    {
      for (j in 1:length(olap))
      {
        chr=as.character(seqnames(olap[j,]))
        olap1=subsetByOverlaps(gr_cytobands,olap[j,])
        if (length(olap1)>0)
        {
          idx1=which.max(width(olap1))
          tmpcytoband=paste0(substr(chr,4,nchar(chr)),as.character(mcols(olap1)["name"][,1])[idx1])
          overlapseg=rbind(overlapseg,data.frame(cytoband=tmpcytoband,chr=chr,start=start(olap[j,]),end=end(olap[j,])))
        }
      }
    }
  }
  num_segdata1=num_segdata2=num_overlapseg=0
  print("overlapsegs:")
  for (i in 1:length(gr_segdata))
  {
    olap1=intersect(gr_segdata[i,],gr_segdata1)
    olap2=intersect(gr_segdata[i,],gr_segdata2)
    if (length(olap1)>0) num_segdata1=num_segdata1+1
    if (length(olap2)>0) num_segdata2=num_segdata2+1
    if (length(olap1)>0 & length(olap2)>0)
    {
      num_overlapseg=num_overlapseg+1
      print(gr_segdata[i,])
    }
  }
  
  width_segdata=as.numeric(format(sum(width(gr_segdata))/10^6,digits = 3,nsmall=1))
  width_segdata1=as.numeric(format(sum(width(gr_segdata1))/10^6,digits = 3,nsmall=1))
  width_segdata2=as.numeric(format(sum(width(gr_segdata2))/10^6,digits = 3,nsmall=1))
  gr_overlapseg=GRanges(seqnames=overlapseg$chr,ranges=IRanges(start=overlapseg$start,end=overlapseg$end),cytoband=overlapseg$cytoband)
  width_overlapseg=as.numeric(format(sum(width(reduce(gr_overlapseg)))/10^6,digits = 3,nsmall=1))
  result=list(allsegs=gr_segdata,overlapseg=overlapseg,width_segdata=width_segdata,width_segdata1=width_segdata1,width_segdata2=width_segdata2,width_overlapseg=width_overlapseg,
              num_segdata1=num_segdata1,num_segdata2=num_segdata2,num_overlapseg=num_overlapseg)
}

detectcancergenes=function(genes)
{
  cancerfile="/fh/fast/dai_j/CancerGenomics/Tools/GISTIC/cancer_gene_census.csv"
  cancergenetable=read.table(cancerfile,header=T,sep=',',stringsAsFactors=F,quote="\"")
  cancergenes=cancergenetable[,1]
  if (class(genes)=="character" & length(genes)==1) genes=unlist(strsplit(genes,",",fixed=T))
  if (is.data.frame(genes)) genes=as.character(genes[,1])
  res=intersect(genes,cancergenes)
  if (length(res)==0)
    res=""
  cancerfile1="/fh/fast/dai_j/CancerGenomics/Tools/database/other/Cancercandidategenes.txt"
  cancergenes1=read.table(cancerfile1,fill=T,stringsAsFactors = F)
  cancergenes1=unlist(cancergenes1)
  cancergenes1=cancergenes1[cancergenes1!=""]
  res1=intersect(genes,cancergenes1)
  if (length(res1)==0) res1=""
  oncogenetable=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/ongene_human.txt",header=T,sep="\t",quote="",stringsAsFactors = F)
  oncogenes=oncogenetable$OncogeneName
  res2=intersect(genes,oncogenes)
  if (length(res2)==0) res2=""
  library(gdata)
  suppressortable=read.xls("/fh/fast/dai_j/CancerGenomics/Tools/database/other/Human_716_TSGs.txt.xlsx",skip=1)
  suppressorgenes=suppressortable$Gene_symbol
  res3=intersect(genes,suppressorgenes)
  if (length(res3)==0) res3=""
  res=list(cosmic=res,candidate=res1,oncogene=res2,suppressorgene=res3)
}

printsegdata=function(segdata,output=NULL)
{
  if (! "gene" %in% colnames(segdata))
  {
    if (ncol(segdata)>4)
    {
      segdata=cbind.data.frame(segdata[,1:5],gene=rep("",nrow(segdata)),cosmicgene=rep("",nrow(segdata)),candidategene=rep("",nrow(segdata)))
    }else
    {
      segdata=cbind.data.frame(segdata,gene=rep("",nrow(segdata)),cosmicgene=rep("",nrow(segdata)),candidategene=rep("",nrow(segdata)))
    }
    
    segdata$gene=as.character(segdata$gene)
    for (i in 1:nrow(segdata))
    {
      segdata$gene[i]=findgenes(chr=segdata$chr[i],start=segdata$start[i],end=segdata$end[i])
    }
    for (i in 1:nrow(segdata))
    {
      tmp=detectcancergenes(segdata$gene[i])
      segdata$cosmicgene[i]=paste0(tmp$cosmic,collapse = ",")
      segdata$candidategene[i]=paste0(tmp$candidate,collapse = ",")
      segdata$oncogene[i]=paste0(tmp$oncogenes,collapse = ",")
      segdata$suppressorgene[i]=paste0(tmp$suppressorgene,collapse = ",")
    }
  }
  
  segdata$start[segdata$start<1000]=0
  segdata$start=format(segdata$start/10^6,nsmall=2,digits = 2)
  segdata$end[segdata$end<1000]=0
  segdata$end=format(segdata$end/10^6,nsmall=2,digits = 2)
  if (is.null(output))
  {
    print(segdata)
  }else
  {
    write.table(segdata,file=output,sep="\t",row.names=F,quote=F)  
  }
}

mart=useMart("ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org/biomart/martservice/", dataset="hsapiens_gene_ensembl") #GRCh37
findgenes=function(chr,start,end)
{
  attributes <- c("hgnc_symbol")
  filters <- c("chromosome_name","start","end")
  if (grepl('chr',chr))
    chr=substr(chr,4,nchar(as.character(chr)))
  
  values=list(chromosome=chr,start=start,end=end)
  genes <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
  genes=unique(genes[,1])
  idxkeep=which(genes!="")
  genes=genes[idxkeep]
  
  genes1=NULL
  if (length(genes)>0)
  {
    if (length(genes)==1)
    {
      genes1=genes
    }else
    {
      for (j in 1:length(genes))
      {
        if (j==1)
        {
          genes1=genes[j]
        }
        else
        {
          genes1=paste(genes1,genes[j],sep=",")
        }
      }
    }
  }
  if (is.null(genes1)) genes1=NA
  return(genes1)
}


segonlyinseg1=function(segdata1,segdata2)
{
  gr_segdata1=GRanges(seqnames=segdata1$chr,ranges=IRanges(start=segdata1$start,end=segdata1$end),cytoband=segdata1$cytoband)
  #gr_segdata1=reduce(gr_segdata1)
  gr_segdata2=GRanges(seqnames=segdata2$chr,ranges=IRanges(start=segdata2$start,end=segdata2$end),cytoband=segdata2$cytoband)
  #gr_segdata2=reduce(gr_segdata2)
  seg1only=data.frame(matrix(NA,nrow=0,ncol=7))
  colnames(seg1only)=c("cytoband","chr","start","end","gene","cosmicgene","candidategene")
  for (i in 1:length(gr_segdata1))
  {
    olap=subsetByOverlaps(gr_segdata1[i,],gr_segdata2,maxgap=3000000)
    if (length(olap)==0)
    {
      #chr=as.character(seqnames(gr_segdata1[i,]))
      #start=start(gr_segdata1[i,])
      #end=end(gr_segdata1[i,])
      #tmp=data.frame(cytoband=findcytoband(chr=chr,start=start,end=end),chr=chr,start=start,end=end)
      tmp=data.frame(cytoband=segdata1$cytoband[i],chr=segdata1$chr[i],start=segdata1$start[i],end=segdata1$end[i])
      seg1only=rbind.data.frame(seg1only,tmp)
    }
  }
  for (i in 1:nrow(seg1only))
  {
    seg1only$gene[i]=findgenes(chr=seg1only$chr[i],start=seg1only$start[i],end=seg1only$end[i])
  }
  for (i in 1:nrow(seg1only))
  {
    tmp=detectcancergenes(seg1only$gene[i])
    seg1only$cosmicgene[i]=paste0(tmp$cosmic,collapse = ",")
    seg1only$candidategene[i]=paste0(tmp$candidate,collapse = ",")
    seg1only$oncogene[i]=paste0(tmp$oncogenes,collapse = ",")
    seg1only$suppressorgene[i]=paste0(tmp$suppressorgene,collapse = ",")
  }
  return(seg1only)
}

RCircos.Get.Plot.Data1=function (genomic.data, plot.type) 
{
  genomic.data <- RCircos.Validate.Genomic.Data(genomic.data, plot.type)
  #data.points <- rep(0, nrow(genomic.data))
  data.points=data.frame(matrix(0,nrow=nrow(genomic.data),ncol=2))
  for (a.row in 1:nrow(genomic.data)) {
    chromosome <- as.character(genomic.data[a.row, "chr"])
    data.points[a.row,1] <- RCircos.Data.Point(chromosome, genomic.data[a.row,"start"]) #for start
    data.points[a.row,2] <- RCircos.Data.Point(chromosome, genomic.data[a.row,"end"])
  }
  genomic.data["Location_start"] <- data.points[,1]
  genomic.data["Location_end"] <- data.points[,2]
  genomic.data <- genomic.data[order(genomic.data$Location_start), 
                               ]
  return(genomic.data)
}
#use the actual segment width
RCircos.Heatmap.Plot2=function(heatmap.data, track.num, side="in") 
{
  heatmap.data=heatmap.data[,c("chr","start","end","col")]
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  heatmap.data <- RCircos.Get.Plot.Data1(heatmap.data, "plot")
  #ColorRamp <- RCircos.Get.Heatmap.ColorScales(RCircos.Par$heatmap.color)
  #columns <- 5:(ncol(heatmap.data) - 1)
  #min.value <- min(as.matrix(heatmap.data[, columns]))
  #max.value <- max(as.matrix(heatmap.data[, columns]))
  #min.value=max.value=1
  #ColorLevel <- seq(min.value, max.value, length = length(ColorRamp))
  #heatmap.locations <- as.numeric(heatmap.data[, ncol(heatmap.data)])
  start=heatmap.data$Location_start
  end=heatmap.data$Location_end
  for (i in 1:length(start))
  {
    if (start[i]==end[i])
    {
      start[i]=start[i]-RCircos.Par$heatmap.width/2
      end[i]=end[i]+RCircos.Par$heatmap.width/2
    }
  }
  #start <- heatmap.locations - RCircos.Par$heatmap.width/2
  #end <- heatmap.locations + RCircos.Par$heatmap.width/2
  data.chroms <- as.character(heatmap.data[, 1])
  chromosomes <- unique(data.chroms)
  cyto.chroms <- as.character(RCircos.Cyto$Chromosome)
  for (a.chr in 1:length(chromosomes)) {
    cyto.rows <- which(cyto.chroms == chromosomes[a.chr])
    locations <- as.numeric(RCircos.Cyto$Location[cyto.rows])
    chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]]
    chr.end <- max(locations)
    data.rows <- which(data.chroms == chromosomes[a.chr])
    start[data.rows[start[data.rows] < chr.start]] <- chr.start
    end[data.rows[end[data.rows] > chr.end]] <- chr.end
  }
  locations <- RCircos.Track.Positions(side, track.num)
  out.pos <- locations[1]
  in.pos <- locations[2]
  chroms <- unique(RCircos.Cyto$Chromosome)
  for (a.chr in 1:length(chroms)) {
    the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr], 
                            ]
    the.start <- the.chr$Location[1] - the.chr$Unit[1] + 
      1
    the.end <- the.chr$Location[nrow(the.chr)]
    polygon.x <- c(RCircos.Pos[the.start:the.end, 1] * out.pos, 
                   RCircos.Pos[the.end:the.start, 1] * in.pos)
    polygon.y <- c(RCircos.Pos[the.start:the.end, 2] * out.pos, 
                   RCircos.Pos[the.end:the.start, 2] * in.pos)
    polygon(polygon.x, polygon.y, col = "white")
  }
  #heatmap.value <- as.numeric(heatmap.data[, data.col])
  for (a.point in 1:nrow(heatmap.data)) {
    #the.level <- which(ColorLevel >= heatmap.value[a.point])
    cell.color <- heatmap.data$col[a.point]
    the.start <- start[a.point]
    the.end <- end[a.point]
    polygon.x <- c(RCircos.Pos[the.start:the.end, 1] * out.pos, 
                   RCircos.Pos[the.end:the.start, 1] * in.pos)
    polygon.y <- c(RCircos.Pos[the.start:the.end, 2] * out.pos, 
                   RCircos.Pos[the.end:the.start, 2] * in.pos)
    polygon(polygon.x, polygon.y, col = cell.color, border = cell.color,lwd=2)
  }
}

data(UCSC.HG19.Human.CytoBandIdeogram)
#chr.exclude <- NULL
chr.exclude <- c("chrX", "chrY")
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
#how many tracks inside:
tracks.inside <- 2
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
rcircos.params <- RCircos.Get.Plot.Parameters()
#rcircos.params$base.per.unit <- 30000
rcircos.params$track.in.start = 2
RCircos.Reset.Plot.Parameters(rcircos.params)

#processed gistic data
gisticdir1="/fh/scratch/delete30/dai_j/gistic/dulak_ploid2degree3force0_cnv_rx1_conf0.95_armpeel1_brlen0.98_broad1"
alllessionsfile1=paste0(gisticdir1,"/all_lesions.conf_95.txt")
wgstumors1=grp1tumors
wgstumors2=grp2tumors

gisticdir2="/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel1_brlen0.98_broad1"
alllessionsfile2=paste0(gisticdir2,"/all_lesions.conf_95.txt")

ampsegdata1=getwholetable1(alllessionsfile1,opt="AMP")
#printsegdata(ampsegdata1,output="/fh/scratch/delete30/dai_j/gistic/dulak_gistic_ampseg.txt")
nrow(ampsegdata1)
#[1] 22
sum(ampsegdata1$end-ampsegdata1$start)/10^6
#[1] 81.34437
ampsegdata2=getwholetable1(alllessionsfile2,opt="AMP")
#printsegdata(ampsegdata2,output="/fh/scratch/delete30/dai_j/henan/gistic/henan_gistic_ampseg.txt")
nrow(ampsegdata2)
#[1] 17
sum(ampsegdata2$end-ampsegdata2$start)/10^6
#[1] 44.94492
delsegdata1=getwholetable1(alllessionsfile1,opt="DEL")
#printsegdata(delsegdata1,output="/fh/scratch/delete30/dai_j/gistic/dulak_gistic_delseg.txt")
nrow(delsegdata1)
#[1] 22
sum(delsegdata1$end-delsegdata1$start)/10^6
#[1] 44.84071
delsegdata2=getwholetable1(alllessionsfile2,opt="DEL")
#printsegdata(delsegdata2,output="/fh/scratch/delete30/dai_j/henan/gistic/henan_gistic_delseg.txt")
nrow(delsegdata2)
#[1] 17
sum(delsegdata2$end-delsegdata2$start)/10^6
#[1] 25.13484

uniqampsegdata1=segonlyinseg1(ampsegdata1,ampsegdata2)
#printsegdata(uniqampsegdata1,output="/fh/scratch/delete30/dai_j/gistic/dulak_gistic_uniqampseg.txt")
uniqampsegdata2=segonlyinseg1(ampsegdata2,ampsegdata1)
#printsegdata(uniqampsegdata2,output="/fh/scratch/delete30/dai_j/henan/gistic/henan_gistic_uniqampseg.txt")
uniqdelsegdata1=segonlyinseg1(delsegdata1,delsegdata2)
#printsegdata(uniqdelsegdata1,output="/fh/scratch/delete30/dai_j/gistic/dulak_gistic_uniqdelseg.txt")
uniqdelsegdata2=segonlyinseg1(delsegdata2,delsegdata1)
#printsegdata(uniqdelsegdata2,output="/fh/scratch/delete30/dai_j/henan/gistic/henan_gistic_uniqdelseg.txt")



ampsegdata1=cbind.data.frame(ampsegdata1,col=rep("red",nrow(ampsegdata1)))
ampsegdata1$col=as.character(ampsegdata1$col)
ampsegdata1=ampsegdata1[,-1]
ampsegdata2=cbind.data.frame(ampsegdata2,col=rep("red",nrow(ampsegdata2)))
ampsegdata2$col=as.character(ampsegdata2$col)
ampsegdata2=ampsegdata2[,-1]
delsegdata1=cbind.data.frame(delsegdata1,col=rep("blue",nrow(delsegdata1)))
delsegdata1$col=as.character(delsegdata1$col)
delsegdata1=delsegdata1[,-1]
delsegdata2=cbind.data.frame(delsegdata2,col=rep("blue",nrow(delsegdata2)))
delsegdata2$col=as.character(delsegdata2$col)
delsegdata2=delsegdata2[,-1]
#Figure5
#circos plot
outfig="EA_gistic_circosplot.png"
png(outfig, width = 8, height = 4, units = 'in', res=300)
par(mfrow=c(1,2))
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Heatmap.Plot2(ampsegdata1, track.num=1)
RCircos.Heatmap.Plot2(ampsegdata2, track.num=2)
mtext("Amplification",side=3)
mtext("(a)",side=1)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
RCircos.Heatmap.Plot2(delsegdata1, track.num=1)
RCircos.Heatmap.Plot2(delsegdata2, track.num=2)
mtext("Deletion",side=3)
mtext("(b)",side=1)
dev.off()

#check armlevel
readarmlevel=function(armfile)
{
  tmp=read.table(armfile,header=T,sep="\t")
  idxamp=which(tmp$Amp.q.value<0.1)
  idxdel=which(tmp$Del.q.value<0.1)
  res=data.frame(matrix(NA,nrow=0,ncol=2))
  colnames(res)=c("arm","type")
  
  if (length(idxamp)>0)
  {
    tmp1=data.frame(arm=tmp$Arm[idxamp],type="amp")
    res=rbind.data.frame(res,tmp1)
  }
  if (length(idxdel)>0)
  {
    tmp1=data.frame(arm=tmp$Arm[idxdel],type="del")
    res=rbind.data.frame(res,tmp1)
  }
  return(res)
}
grp1arm=readarmlevel(armfile="/fh/scratch/delete30/dai_j/gistic/dulak_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1/broad_significance_results.txt")
grp2arm=readarmlevel(armfile="/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1/broad_significance_results.txt")

grp1arm=readarmlevel(armfile="/fh/scratch/delete30/dai_j/gistic/dulak_ploid2degree3force0_cnv_rx1_conf0.95_armpeel1_brlen0.98_broad1/broad_significance_results.txt")
grp2arm=readarmlevel(armfile="/fh/scratch/delete30/dai_j/henan/gistic/henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel1_brlen0.98_broad1/broad_significance_results.txt")

ampolapseg12=countoverlap2segs(ampsegdata1,ampsegdata2)
nrow(ampolapseg12$overlapseg)
unique(as.character(ampolapseg12$overlapseg$cytoband))
ampolapseg12$width_overlapseg/ampolapseg12$width_segdata
delolapseg12=countoverlap2segs(delsegdata1,delsegdata2)
nrow(delolapseg12$overlapseg)
unique(as.character(delolapseg12$overlapseg$cytoband))
delolapseg12$width_overlapseg/delolapseg12$width_segdata
delolapseg12$width_segdata1
delolapseg12$width_segdata2

#combined 2 datasets
gisticdir="/fh/scratch/delete30/dai_j/henan/gistic2data/dulak_henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel0_brlen0.98_broad1"
#gisticdir="/fh/scratch/delete30/dai_j/henan/gistic2data/dulak_henan_ploid2degree3force0_cnv_rx1_conf0.95_armpeel1_brlen0.98_broad1"

alllessionsfile=paste0(gisticdir,"/all_lesions.conf_95.txt")

ampsegdata=getwholetable1(alllessionsfile,opt="AMP",qcutoff = 0.01)
delsegdata=getwholetable1(alllessionsfile,opt="DEL",qcutoff = 0.01)
colnames(ampsegdata)[6:21]=paste0("US-EAC",1:16)
colnames(ampsegdata)[22:ncol(ampsegdata)]=paste0("CH-EAC",1:10)
colnames(delsegdata)=colnames(ampsegdata)
#only keep 2, set 1 as 0. >2^0.9 or <2^(-0.9)
useampordel=function(segdata)
{
  for (i in 1:nrow(segdata))
  {
    for (j in 6:ncol(segdata))
    {
      if (segdata[i,j]<2)
        segdata[i,j]=0
    }
  }
  #remove high frequency regions < 2/3 #samples
  num=apply(segdata[,6:ncol(segdata)],1,function(x){
    sum(x>1)
  })
  segdata=segdata[num<2/3*(ncol(segdata)-6) & num>0,]
  return(segdata)
}
library("ComplexHeatmap")
library(circlize)
Heatmap(ampsegdata[,6:ncol(ampsegdata)],cluster_rows = T,name="amp event",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","orange","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(0,1,2), labels = c("none", "gain","amplification")),
        row_names_gp = gpar(fontsize = 10)
)
Heatmap(delsegdata[,6:ncol(delsegdata)],cluster_rows = T,name="del event",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","orange","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(0,1,2), labels = c("none", "loss","deletion")),
        row_names_gp = gpar(fontsize = 10)
)
ampsegdata_=useampordel(ampsegdata)
Heatmap(ampsegdata_[,6:ncol(ampsegdata_)],cluster_rows = T,name="amp event",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","orange","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(0,1), labels = c("none", "amplification")),
        row_names_gp = gpar(fontsize = 10)
)
delsegdata_=useampordel(delsegdata)
Heatmap(delsegdata_[,6:ncol(delsegdata_)],cluster_rows = T,name="del event",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","orange","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(0,1), labels = c("none", "deletion")),
        row_names_gp = gpar(fontsize = 10)
)
allsegdata=rbind.data.frame(ampsegdata_,delsegdata_)
for (i in (nrow(ampsegdata_)+1):nrow(allsegdata))
{
  for (j in 6:ncol(allsegdata))
  {
    if (allsegdata[i,j]==2)
    {
      allsegdata[i,j]=-2
    }
    if (allsegdata[i,j]==1)
    {
      allsegdata[i,j]=-1
    }  
  }
}
tmp=which(rownames(allsegdata)=="21p11.2(1)")
rownames(allsegdata)[tmp]="21p11.2"
tmp=which(rownames(allsegdata)=="21p11.2(2)")
rownames(allsegdata)[tmp]="21p11.2(1)"

Heatmap(allsegdata[,6:ncol(allsegdata)],cluster_rows = T,name="Event",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","gray88","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(0,2,-2), labels = c("NONE", "AMP","DEL"),title_gp = gpar(fontsize = 14),labels_gp = gpar(fontsize = 12)),
        row_names_gp = gpar(fontsize = 10)
)

#compare 1d with 2d
uniq_ampsegdata=segonlyinseg1(ampsegdata,rbind.data.frame(ampsegdata1[,1:4],ampsegdata2[,1:4]))
uniq_ampsegdata_=segonlyinseg1(ampsegdata_,rbind.data.frame(ampsegdata1[,1:4],ampsegdata2[,1:4]))
uniq_delsegdata=segonlyinseg1(delsegdata,rbind.data.frame(delsegdata1[,1:4],delsegdata2[,1:4]))
uniq_delsegdata_=segonlyinseg1(delsegdata_,rbind.data.frame(delsegdata1[,1:4],delsegdata2[,1:4]))


#consider all the events
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

Heatmap(ampsegdata[,6:ncol(ampsegdata)],cluster_rows = T,name="amp event",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","orange","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(-2,-1,0,1,2), labels = c("loss","deletion","none","gain","amplification")),
        row_names_gp = gpar(fontsize = 10)
)
Heatmap(delsegdata[,6:ncol(delsegdata)],cluster_rows = T,name="del event",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),#col=c("blue","orange","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(-2,-1,0,1,2), labels = c("loss","deletion","none","gain","amplification")),
        row_names_gp = gpar(fontsize = 10)
)
ampsegdata_=useampordel1(ampsegdata)
Heatmap(ampsegdata_[,6:ncol(ampsegdata_)],cluster_rows = T,name="amp event",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","orange","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(0,-2,2), labels = c("none","deletion", "amplification")),
        row_names_gp = gpar(fontsize = 10)
)
delsegdata_=useampordel1(delsegdata,type="del")
Heatmap(delsegdata_[,6:ncol(delsegdata_)],cluster_rows = T,name="del event",show_column_names=T,show_column_dend=T,
        rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","orange","red"),
        column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
        heatmap_legend_param = list(at = c(0,-2,2), labels = c("none", "deletion","amplification")),
        row_names_gp = gpar(fontsize = 10)
)
allsegdata=rbind.data.frame(ampsegdata_,delsegdata_)

tmp=which(rownames(allsegdata)=="21p11.2(1)")
rownames(allsegdata)[tmp]="21p11.2"
tmp=which(rownames(allsegdata)=="21p11.2(2)")
rownames(allsegdata)[tmp]="21p11.2(1)"
tmp=which(rownames(allsegdata)=="19q13.331")
rownames(allsegdata)[tmp]="19q13.33(1)"

#Figure6
outfig="EA_gistic_heatmap.png"
png(outfig, width = 8, height = 8, units = 'in', res=300)
#figure4(c)
plotfig4c=function()
{
  Heatmap(allsegdata[,6:ncol(allsegdata)],cluster_rows = F,name="Event",show_column_names=T,show_column_dend=T,
          rect_gp = gpar(col = "white", lty = 1, lwd = 1),col=c("blue","gray88","red"),
          column_dend_height = unit(30, "mm"),row_dend_width  = unit(30, "mm"),
          heatmap_legend_param = list(at = c(0,2,-2), labels = c("None", "Amplification","Deletion"),title_gp = gpar(fontsize = 12),labels_gp = gpar(fontsize = 12)),
          row_names_gp = gpar(fontsize = 12)
  )
}

dev.off()

#plot figure4
outfig="EA_figure4.ps"
postscript(outfig, width = 12, height = 20)
layout(mat=matrix(c(1,2,3,3), ncol=2, byrow=TRUE))
par(mar=c(0.5,1.1,1.1,1.1))
plotfig4a()
plotfig4b()
plotfig4c()
dev.off()

#compare 1d with 2d
uniq_ampsegdata=segonlyinseg1(ampsegdata,rbind.data.frame(ampsegdata1[,1:4],ampsegdata2[,1:4]))
uniq_ampsegdata_=segonlyinseg1(ampsegdata_,rbind.data.frame(ampsegdata1[,1:4],ampsegdata2[,1:4]))
uniq_delsegdata=segonlyinseg1(delsegdata,rbind.data.frame(delsegdata1[,1:4],delsegdata2[,1:4]))
uniq_delsegdata_=segonlyinseg1(delsegdata_,rbind.data.frame(delsegdata1[,1:4],delsegdata2[,1:4]))
nrow(ampsegdata)
nrow(delsegdata)
nrow(ampsegdata_)
nrow(delsegdata_)

#SV
#generated by read_manta.R, save(mat1,mat2,grp1table,grp2table,genes1,genes2,grp1tumors,grp2tumors,file="mantares.Rdata")
load("/fh/fast/dai_j/CancerGenomics/Tools/wang/SV/mantares1.Rdata")
tmp=detectcancergenes(rownames(mat2)[1:11])
tmp=detectcancergenes(rownames(mat1)[1:7])
mat=mat1 #US
mat=mat2 #CH
num_sv=sapply(1:nrow(mat),function(i){
  sum(mat[i,] != " ")
})


CFSgenes=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/CFSgenes.txt",stringsAsFactors = F)
idx_CFS=rownames(mat) %in% CFSgenes[,1]
#CFS genes
mat_CFS=mat[idx_CFS,]
num_sv=sapply(1:nrow(mat_CFS),function(i){
  sum(mat_CFS[i,] != " ")
})
idx_rowCFS=num_sv/ncol(mat_CFS)>=0.2
idx_rowCFS=num_sv/ncol(mat_CFS)>=0.5 #mat1
sum(idx_rowCFS)
tmp=data.frame(gene=rownames(mat_CFS)[idx_rowCFS],prop=num_sv[idx_rowCFS]/ncol(mat_CFS))
tmp1=detectcancergenes(tmp)

mat_noneCFS=mat[!idx_CFS,]
num_sv=sapply(1:nrow(mat_noneCFS),function(i){
  sum(mat_noneCFS[i,] != " ")
})
idx_nonerowCFS=num_sv/ncol(mat_noneCFS)>=0.2
sum(idx_nonerowCFS)
tmp=data.frame(gene=rownames(mat_noneCFS)[idx_nonerowCFS],prop=num_sv[idx_nonerowCFS]/ncol(mat_noneCFS))
tmp1=detectcancergenes(tmp)

#plot figure5a
library(ComplexHeatmap)
mutation_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  DEL=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "blue", col = NA),just="top")
  },
  DUP=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "red", col = NA),just="bottom")
  },
  INV=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  },
  INS=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "black", col = NA),just="centre")
  },
  BND=function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "#F88000", col = NA),just="centre")
  }
)

col = c("DEL" = "blue", "DUP"="red", "INV"="#008000", "INS"="black", "BND"="#F88000")

num_sv=sapply(1:nrow(mat),function(i){
  sum(mat[i,] != " ")
})
# ha = HeatmapAnnotation(SNV= anno_barplot(num_sv,axis=T,axis_gp = gpar(fontsize = 12),axis_side="right",border=F,ylim=c(0,500)),
#                        show_annotation_name = T,annotation_name_offset = unit(2, "cm"),gap = unit(3, "mm"))
png("fig5a.png", width = 4, height = 6, units = 'in', res=300)
oncoPrint(as.matrix(rbind.data.frame(mat_noneCFS[idx_nonerowCFS,],mat_CFS[idx_rowCFS,])), get_type = function(x) strsplit(x, ";")[[1]],
          row_order = NULL,column_order = NULL,
          remove_empty_columns = FALSE,
          alter_fun = mutation_fun, col = col,
          pct_gp=gpar(fontsize=12),
          axis_gp = gpar(fontsize = 12),
          show_row_barplot=F,
          # column_title = name,
          # bottom_annotation=ha,
          # bottom_annotation_height=unit(3,"cm"),
          heatmap_legend_param = list(title = "SV", at = c("DEL", "DUP", "INV","INS","BND"), 
                                      labels = c("DEL", "DUP", "INV","INS","BND"))
)
dev.off()

png("fig5b.png", width = 4, height = 6, units = 'in', res=300)
oncoPrint(as.matrix(rbind.data.frame(mat_noneCFS[idx_nonerowCFS,],mat_CFS[idx_rowCFS,])), get_type = function(x) strsplit(x, ";")[[1]],
          row_order = NULL,column_order = NULL,
          remove_empty_columns = FALSE,
          alter_fun = mutation_fun, col = col,
          pct_gp=gpar(fontsize=12),
          axis_gp = gpar(fontsize = 12),
          show_row_barplot=F,
          # column_title = name,
          # bottom_annotation=ha,
          # bottom_annotation_height=unit(3,"cm"),
          heatmap_legend_param = list(title = "SV", at = c("DEL", "DUP", "INV","INS","BND"), 
                                      labels = c("DEL", "DUP", "INV","INS","BND"))
)
dev.off()

png("fig5b.png", width = 4, height = 6, units = 'in', res=300)
oncoPrint(as.matrix(mat_noneCFS[idx_nonerowCFS,]), get_type = function(x) strsplit(x, ";")[[1]],
          row_order = NULL,column_order = NULL,
          remove_empty_columns = FALSE,
          alter_fun = mutation_fun, col = col,
          pct_gp=gpar(fontsize=12),
          axis_gp = gpar(fontsize = 12),
          show_row_barplot=F,
          # column_title = name,
          # bottom_annotation=ha,
          # bottom_annotation_height=unit(3,"cm"),
          heatmap_legend_param = list(title = "SV", at = c("DEL", "DUP", "INV","INS","BND"), 
                                      labels = c("DEL", "DUP", "INV","INS","BND"))
)
dev.off()

png("fig5c.png", width = 4, height = 6, units = 'in', res=300)
oncoPrint(as.matrix(mat_CFS[idx_rowCFS,]), get_type = function(x) strsplit(x, ";")[[1]],
          row_order = NULL,column_order = NULL,
          remove_empty_columns = FALSE,
          alter_fun = mutation_fun, col = col,
          pct_gp=gpar(fontsize=12),
          axis_gp = gpar(fontsize = 12),
          show_row_barplot=F,
          # column_title = name,
          # bottom_annotation=ha,
          # bottom_annotation_height=unit(3,"cm"),
          heatmap_legend_param = list(title = "SV", at = c("DEL", "DUP", "INV","INS","BND"), 
                                      labels = c("DEL", "DUP", "INV","INS","BND"))
)
dev.off()

restable=cbind(restable,numsv=c(grp1numsv,grp2numsv))
par(mar=c(3.5,5.1,3.1,2.1))
par(mfrow=c(1,1))
outfig="EA_numsv.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
plot2groups(item="numsv",datatable=restable,adjustymin=T,main='',ylab="",usewilcox=T)
dev.off()
library(gdata)
dulakstable4=read.table(file="/fh/fast/dai_j/CancerGenomics/EAC/Dulak_Stable4.txt",header=T,sep="\t",quote="",fill=T)
dulakstable4$Sample=as.character(dulakstable4$Sample)
dulakstable4$Sample=paste0("ESO-",dulakstable4$Sample)
numsv_dulak=rep(0,length(grp1tumors))
for (i in 1:length(grp1tumors))
{
  numsv_dulak[i]=sum(dulakstable4$Sample==as.character(srrdulakstable2$Individual[i]))
}
draw2barplot(data1=numsv_dulak,data2=grp1numsv)

#draw manta result
drawmanta=function(mantatable=grp1table)
{
  chrs=paste0("chr",1:22)
  mantatable=mantatable[mantatable$chr1 %in% chrs & mantatable$chr2 %in% chrs,]
  deltable=data.frame(matrix(NA,nrow=sum(mantatable$svtype=="DEL"),ncol=3))
  colnames(deltable)=c("chr","start","end")
  idxdel=which(mantatable$svtype=="DEL")
  deltable$chr=mantatable$chr1[idxdel]
  deltable$start=mantatable$pos[idxdel]
  deltable$end=mantatable$end[idxdel]
  deltable$col="black"
  duptable=data.frame(matrix(NA,nrow=sum(mantatable$svtype=="DUP"),ncol=3))
  colnames(duptable)=c("chr","start","end")
  idxdup=which(mantatable$svtype=="DUP")
  duptable$chr=mantatable$chr1[idxdup]
  duptable$start=mantatable$pos[idxdup]
  duptable$end=mantatable$end[idxdup]
  duptable$col="red"
  instable=data.frame(matrix(NA,nrow=sum(mantatable$svtype=="INS"),ncol=3))
  colnames(instable)=c("chr","start","end")
  idxins=which(mantatable$svtype=="INS")
  instable$chr=mantatable$chr1[idxins]
  instable$start=mantatable$pos[idxins]
  instable$end=mantatable$end[idxins]
  instable$col="pink"
  invtable=data.frame(matrix(NA,nrow=sum(mantatable$svtype=="INV"),ncol=3))
  colnames(invtable)=c("chr","start","end")
  idxinv=which(mantatable$svtype=="INV")
  invtable$chr=mantatable$chr1[idxinv]
  invtable$start=mantatable$pos[idxinv]
  invtable$end=mantatable$end[idxinv]
  invtable$col="orange"
  bndtable=data.frame(matrix(NA,nrow=sum(mantatable$svtype=="BND"),ncol=6))
  colnames(bndtable)=c("chr1","start1","end1","chr2","start2","end2")
  idxbnd=which(mantatable$svtype=="BND")
  bndtable$chr1=mantatable$chr1[idxbnd]
  bndtable$start1=bndtable$end1=mantatable$pos[idxbnd]
  bndtable$chr2=mantatable$chr2[idxbnd]
  bndtable$start2=bndtable$end2=mantatable$end[idxbnd]
  bndtable$col="blue"
  RCircos.Set.Plot.Area();
  RCircos.Chromosome.Ideogram.Plot()
  RCircos.Heatmap.Plot2(deltable, track.num=1)
  RCircos.Heatmap.Plot2(duptable, track.num=2)
  RCircos.Heatmap.Plot2(invtable, track.num=3)
  RCircos.Heatmap.Plot2(instable, track.num=4)
  RCircos.Link.Plot(bndtable, track.num=5)
  bndchrmatrix=data.frame(matrix(nrow=length(chrs),ncol=length(chrs)))
  rownames(bndchrmatrix)=colnames(bndchrmatrix)=chrs
  for (i in 1:length(chrs))
  {
    for (j in 1:length(chrs))
    {
      bndchrmatrix[i,j]=sum(bndtable$chr1==chrs[i] & bndtable$chr2==chrs[j])+sum(bndtable$chr1==chrs[j] & bndtable$chr2==chrs[i])
    }
  }
}
drawmanta(mantatable=grp1table)
drawmanta(mantatable=grp2table)


##lego plot
library(rgl)
stackplot.3d=function(x,y,z,alpha=1,topcol="#078E53",sidecol="#aaaaaa"){
  
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to draw the sides and ends of the column separately
  #save on.exit(par3d(save))
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  
  ## Determine the coordinates of each surface of the column and its edges
  x1=c(rep(c(x[1],x[2],x[2],x[1]),3),rep(x[1],4),rep(x[2],4))
  z1=c(rep(0,4),rep(c(0,0,z,z),4))
  y1=c(y[1],y[1],y[2],y[2],rep(y[1],4),rep(y[2],4),rep(c(y[1],y[2],y[2],y[1]),2))
  x2=c(rep(c(x[1],x[1],x[2],x[2]),2),rep(c(x[1],x[2],rep(x[1],3),rep(x[2],3)),2))
  z2=c(rep(c(0,z),4),rep(0,8),rep(z,8) )
  y2=c(rep(y[1],4),rep(y[2],4),rep(c(rep(y[1],3),rep(y[2],3),y[1],y[2]),2) )
  
  ## These lines create the sides of the column and its coloured top surface
  rgl.quads(x1,z1,y1,col=rep(sidecol,each=4),alpha=alpha,lit=F)
  rgl.quads(c(x[1],x[2],x[2],x[1]),rep(z,4),c(y[1],y[1],y[2],y[2]),
            col=rep(topcol,each=4),alpha=1,lit=F)
  ## This line adds black edges to the column
  rgl.lines(x2,z2,y2,col="#000000",lit=F)
}
context3d=function(z,alpha=1,scalexy=10,scalez=1,gap=0.2,filename="graph"){
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to add each column sequentially
  #save on.exit(par3d(save))
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  ## Recreate Broad order
  types=c("C.G.G.C","T.A.A.T","C.A.G.T","T.G.A.C","C.T.G.A","T.C.A.G")
  contexts=c("TxT","CxT","AxT","GxT","TxC","CxC","AxC","GxC",
             "TxA","CxA","AxA","GxA","TxG","CxG","AxG","GxG")
  typeorder=c()
  for(type in types){
    typeorder=c(typeorder,paste(type,contexts,sep="_"))
  }
  z=z[typeorder]
  
  ## Reorder data into 6 regions
  set1=c(1:4,17:20,5:8,21:24,9:12,25:28,13:16,29:32)
  set2=set1+32
  set3=set1+64
  neworder=c(set1,set2,set3)
  
  ## Define dimensions of the plot
  dimensions=c(12,8)
  
  ## Scale column area and the gap between columns
  y=seq(1,dimensions[1])*scalexy
  x=seq(1,dimensions[2])*scalexy
  gap=gap*scalexy
  
  ## Scale z coordinate
  z=z*scalez
  
  ## Set up colour palette
  broadcolors=c("#805D3F","#72549A","#5EAFB2","#3F4F9D","#F2EC3C","#74B655")
  colors=as.vector(sapply(broadcolors,rep,16))
  
  ## Plot each of the columns
  for(i in 1:dimensions[1]){
    for(j in 1:dimensions[2]){
      it=(i-1)*dimensions[2]+j # Variable to work out which column to plot; counts from 1:96
      stackplot.3d(c(gap+x[j],x[j]+scalexy),
                   c(-gap-y[i],-y[i]-scalexy),
                   z[neworder[it]],
                   alpha=alpha,
                   topcol=colors[neworder[it]],
                   sidecol=colors[neworder[it]])
    }
  }
  ## Set the viewpoint and add axes and labels
  rgl.viewpoint(theta=50,phi=40,fov=0,scale=c(1,0.9,1))
  #axes3d("y-+",labels=TRUE,at=c(0,20,40,60))
  ats=seq(0,ceiling(max(counts)/10)*10,10)
  axes3d("y-+",labels=TRUE,at=ats)
  filename=paste0(filename,".ps")
  rgl.postscript(filename, fmt="ps")
}
library(BSgenome.Hsapiens.UCSC.hg19)
genome=getSeq(BSgenome.Hsapiens.UCSC.hg19)
trinuctable=trinucleotideFrequency(genome)
compute_trinucleotide_freq1=function(counts0,trinuctable,numsample)
{
  res=counts0
  for (i in 1:length(counts0))
  {
    myname=names(counts0)[i]
    tmp=unlist(strsplit(myname,"_"))
    tmp1=unlist(strsplit(tmp,".",fixed=T))
    tmp2=unlist(strsplit(tmp1[5],"x"))
    tri1=paste0(tmp2[1],tmp1[1],tmp2[2])
    tri2=paste0(tmp2[1],tmp1[3],tmp2[2])
    count1=sum(trinuctable[,colnames(trinuctable)==tri1])
    count2=sum(trinuctable[,colnames(trinuctable)==tri2])
    count=count1+count2
    res[i]=counts0[i]/count/numsample*10^6
  }
  return(res)
}


plotlego=function(rawdata,numsample,outputprefix)
{
  counts0=as.numeric(rawdata)
  names(counts0)=colnames(rawdata)
  counts=compute_trinucleotide_freq1(counts0,trinuctable,numsample)
  context3d(z=counts,filename=outputprefix)
}
rawdata1=read.table("/fh/scratch/delete30/dai_j/mutect4/all_lego.txt",header=T)
plotlego(rawdata=rawdata1,numsample=length(grp1tumors),outputprefix="lego_dulak")
rawdata2=read.table("/fh/scratch/delete30/dai_j/henan/mutect5/all_lego.txt",header=T)
plotlego(rawdata=rawdata2,numsample=length(grp2tumors),outputprefix="lego_henan")

#check normal snp variations
mutectdata=read.table(file="/fh/scratch/delete30/dai_j/henan/mutect5/11A.Mutect_out.txt",header=T,sep="\t",quote="")
tumorvar=mutectdata$tumor_f
normalvar=mutectdata$normal_f
par(mfrow=c(2,1))
hist(tumorvar,probability = T,xlab="Allele  fraction  of  variant  in tumor",main="CH-EA")
hist(normalvar,probability = T,xlab="Allele  fraction  of  variant  in normal",main="CH-EA")

mutectdata1=read.table(file="/fh/scratch/delete30/dai_j/mutect4/SRR1002929.Mutect_out.txt",header=T,sep="\t",quote="")
tumorvar1=mutectdata1$tumor_f
normalvar1=mutectdata1$normal_f
hist(tumorvar1,probability = T,xlab="Allele  fraction  of  variant  in tumor",main="US-EA")
hist(normalvar1,probability = T,xlab="Allele  fraction  of  variant  in normal",main="US-EA")

sum(normalvar<0.1)/length(normalvar)
sum(normalvar1<0.1)/length(normalvar1)

sum(normalvar>0.9)/length(normalvar)
sum(normalvar1>0.9)/length(normalvar1)
idx=which(normalvar==1)
View(mutectdata[idx,])

checkvafinnormals=function(grptumors=grp1tumors,grpdir=grp1dir)
{
  res=data.frame(matrix(NA,nrow=length(grptumors),ncol=8))
  colnames(res)=c("n_numvar","n_lowprop","n_lowprop1","n_highprop","t_numvar","t_lowprop","t_lowprop1","t_highprop")
  for (i in 1:length(grptumors))
  {
    cat(i,"..")
    mutectdata=read.table(file=paste0(grpdir,"/",grptumors[i],".Mutect_out.txt"),header=T,sep="\t",quote="")
    mutectdata=mutectdata[mutectdata$t_q20_count>=10 & mutectdata$t_q20_count<=100 & mutectdata$n_q20_count>=10 & mutectdata$n_q20_count<=100,]
    res$n_numvar[i]=sum(mutectdata$normal_f>0)
    res$n_lowprop[i]=sum(mutectdata$normal_f<0.1 & mutectdata$normal_f>0)/sum(mutectdata$normal_f>0)
    res$n_lowprop1[i]=sum(mutectdata$normal_f<0.01 & mutectdata$normal_f>0)/sum(mutectdata$normal_f>0)
    res$n_highprop[i]=sum(mutectdata$normal_f>0.9 & mutectdata$normal_f>0)/sum(mutectdata$normal_f>0)
    res$t_numvar[i]=sum(mutectdata$tumor_f>0)
    res$t_lowprop[i]=sum(mutectdata$tumor_f<0.1 & mutectdata$tumor_f>0)/sum(mutectdata$tumor_f>0)
    res$t_lowprop1[i]=sum(mutectdata$tumor_f<0.01 & mutectdata$tumor_f>0)/sum(mutectdata$tumor_f>0)
    res$t_highprop[i]=sum(mutectdata$tumor_f>0.9 & mutectdata$tumor_f>0)/sum(mutectdata$tumor_f>0)
    
    # file=paste0("./tmp/",grptumors[i],".vaf.hist.png")
    # png(file,width = 800, height=800)
    # par(mfrow=c(2,1))
    # hist(mutectdata$tumor_f[mutectdata$tumor_f>0],probability = F,xlab="Allele  fraction  of  variant  in tumor",main="",axes=F)
    # axis(1)
    # axis(2,at=seq(0,sum(mutectdata$tumor_f>0),by=100000))
    # hist(mutectdata$normal_f[mutectdata$normal_f>0],probability = F,xlab="Allele  fraction  of  variant  in normal",main="",axes=F)
    # axis(1)
    # axis(2,at=seq(0,sum(mutectdata$tumor_f>0),by=100000))
    # dev.off()
  }
  return(res)
}

grp1vaf=checkvafinnormals(grptumors=grp1tumors,grpdir=grp1dir)
grp2vaf=checkvafinnormals(grptumors=grp2tumors,grpdir=grp2dir)

restable=cbind(restable,n_numvar=c(grp1vaf$n_numvar,grp2vaf$n_numvar))
#figure1
#par(mar=c(5.1,4.1,4.1,2.1))
par(mar=c(3.5,5.1,3.1,2.1))
par(mfrow=c(1,1))
outfig="EA_n_numvar.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
plot2groups(item="n_numvar",datatable=restable,adjustymin=T,main='Variants in normals',ylab="Number of mutations",usewilcox=T)
dev.off()

restable=cbind(restable,t_numvar=c(grp1vaf$t_numvar,grp2vaf$t_numvar))
plot2groups(item="t_numvar",datatable=restable,adjustymin=T,main='Variants in tumors',ylab="Number of mutations",usewilcox=T)

draw2barplot(data1=grp1vaf$n_numvar,data2=grp1vaf$t_numvar,ylable="Number of mutations",legends=c("normals","tumors"),main="US-EAC")
draw2barplot(data1=grp2vaf$n_numvar,data2=grp2vaf$t_numvar,ylable="Number of mutations",legends=c("normals","tumors"),main="CH-EAC")

restable=cbind(restable,n_lowprop=c(grp1vaf$n_lowprop,grp2vaf$n_lowprop))
plot2groups(item="n_lowprop",datatable=restable,adjustymin=T,main='Variants in normals, VAF<0.1',ylab="Proportion",usewilcox=T)

restable=cbind(restable,n_low=c(grp1vaf$n_lowprop*grp1vaf$n_numvar,grp2vaf$n_lowprop*grp2vaf$n_numvar))
plot2groups(item="n_low",datatable=restable,adjustymin=T,main='Variants in normals, VAF<0.1',ylab="Number of mutations",usewilcox=T)

restable=cbind(restable,n_lowprop1=c(grp1vaf$n_lowprop1,grp2vaf$n_lowprop1))
plot2groups(item="n_lowprop1",datatable=restable,adjustymin=T,main='Variants in normals, VAF<0.01',ylab="Proportion",usewilcox=T)

restable=cbind(restable,n_low1=c(grp1vaf$n_lowprop1*grp1vaf$n_numvar,grp2vaf$n_lowprop1*grp2vaf$n_numvar))
plot2groups(item="n_low1",datatable=restable,adjustymin=T,main='Variants in normals, VAF<0.01',ylab="Number of mutations",usewilcox=T)


# MSI results
msidir="/fh/scratch/delete30/dai_j/henan/msisensor"
readmsi=function(msidir,grptumors=grp2tumors)
{
  res=data.frame(matrix(NA,ncol=2,nrow=length(grptumors)))
  colnames(res)=c("sample","score")
  res$sample=grptumors
  for (i in 1:length(grptumors))
  {
    file=paste0(msidir,"/",grptumors[i])
    tmp=read.table(file,header = T)
    res$score[i]=tmp[1,3]
  }
  return(res)
}
msiscore=readmsi(msidir)
write.table(msiscore,file=paste0(msidir,"/msiscore.txt"),row.names = F,col.names = T,sep="\t",quote=F)

#ConTest results
contestdir="/fh/scratch/delete30/dai_j/henan"
readcontest=function(contestdir,grptumors=grp2tumors)
{
  res=data.frame(matrix(NA,ncol=2,nrow=length(grptumors)))
  colnames(res)=c("sample","contamination")
  res$sample=grptumors
  for (i in 1:length(grptumors))
  {
    file=paste0(contestdir,"/",grptumors[i],".contest.txt")
    tmp=read.table(file,header = T)
    res$contamination[i]=tmp[1,4]
  }
  return(res)
}
contestres=readcontest(contestdir)
write.table(contestres,file=paste0(contestdir,"/contestresult.txt"),row.names = F,col.names = T,sep="\t",quote=F)

#mutation signature:
formatContexts=function (contexts) 
{
  for_plotting = reshape2::melt(contexts)
  colnames(for_plotting)[1] = "sample.id"
  colnames(for_plotting)[2] = "full_context"
  colnames(for_plotting)[3] = "fraction"
  for_plotting$mutation = substr(for_plotting$full_context, 
                                 3, 5)
  for_plotting$full_context = as.factor(for_plotting$full_context)
  for_plotting$full_context = factor(for_plotting$full_context, 
                                     levels = for_plotting$full_context, ordered = T)
  return(for_plotting)
}



plotsignature2=function (signature, name, y_limit=NULL) 
{
  #sort signature subtypes
  mutationsubtypes=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/mutationsignature/mutationsubtypes.txt")
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
  #signature: matrix 1 by 96, rowname: title
  #      A[C>A]A     A[C>A]C      A[C>A]G     A[C>A]T   ...
  # 1 0.001314202 0.001961651 9.687434e-05 0.002529105  ...
 op <- graphics::par()
  signature <- as.matrix(signature)
  if (is.null(y_limit))
  y_limit <- 1.2 * max(signature)
  tumor_plotting <- formatContexts(signature)
  #name <- unique(tumor_plotting$sample.id)
  mutationcols=c("skyblue", "black", "red", "gray", "forestgreen", "tan3")
  # grDevices::palette(c("skyblue", "black", "red", "gray", 
  #                      "forestgreen", "tan3"))
  graphics::barplot(tumor_plotting$fraction, names.arg = "", 
                    cex.names = 0.7, las = 2, col = NA, ylim = c(0, y_limit), 
                    border = NA, xaxt = "n", ann = FALSE, yaxt = "n", space = 0.3)
  graphics::box()
  x = graphics::par("usr")
  # graphics::abline(h = seq(from = 0, to = y_limit, by = 0.01), 
  #                  col = "#d3d3d350", lty = 1)
  # graphics::abline(v = seq(from = x[1], to = x[2], by = 1), 
  #                  col = "#d3d3d350", lty = 1)
  cols=NULL
  for (i in 1:length(mutationtypes))
  {
    cols=c(cols,rep(mutationcols[i],16))
  }
  graphics::barplot(tumor_plotting$fraction, names.arg ="", 
                    cex.names = 0.7, cex.lab=1.5,las = 2, col = cols, 
                    ylim = c(0, y_limit), border = NA, space = 0.3, main = name, 
                    ylab = "fraction", add = TRUE,cex.axis=1.2)
  #   graphics::legend("topright", legend = unique(tumor_plotting$mutation), 
  #                    col = c("skyblue", "black", "red", "gray", "forestgreen", 
  #                            "tan3"), bty = "n", ncol = 6, inset = c(-0, 0), 
  #                    pch = 15, xpd = TRUE, pt.cex = 4)
  mutationtypes=c('C>A','C>G','C>T','A>T','A>G','A>C')
  
  ats=c()
  xs=c(1,3,5,7,9,11)
  len=(par("usr")[2]-par("usr")[1])/(6+2/3)
  for (i in 1:6)
  {
    ats[i]=par("usr")[1]+1/3*len+xs[i]*len/2
  }
  for (j in 1:length(mutationtypes))
  {
    mtext(mutationtypes[j],side=1,line=1,at=ats[j],col=mutationcols[j],cex=1.2)
  }
  on.exit(suppressWarnings(graphics::par(op)))
  return(signature)
}

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
    y_limit <- 1.3 * max(signature)
  
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
   text(x[i],-0.002,substr(trinucleotide[i],1,1),col="black",srt=90,cex=0.5)
   text(x[i],-0.004,substr(trinucleotide[i],2,2),col=cols[i],srt=90,cex=0.5)
   text(x[i],-0.006,substr(trinucleotide[i],3,3),col="black",srt=90,cex=0.5)
 }
 axis(side=2, at=seq(0,y_limit,0.01),cex.axis=1.2,cex.lab=1.2,line=F,lwd.ticks=1.1)
 for (i in seq(from = 0, to = y_limit-0.01, by = 0.01))
 {
   segments(par("usr")[1],i,125,i,col="gray66")
 }
 # abline(h = seq(from = 0, to = y_limit-0.01, by = 0.01), col="gray66")
 y1=max(signature)+0.005+y_limit*0.03
 # y2=y1+(y_limit-max(signature))/4
 y2=y1+0.004
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

signaturemat1=read.table(file="/fh/fast/dai_j/CancerGenomics/Tools/mutationsignature/henan_mutect5_4_signatures.txt")

#plot figure2
outfig="EA_figure2.ps"
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

#read absolute results
absolutedir="/fh/scratch/delete30/dai_j/henan/absolute"

readabsolute=function(absolutedir,tumors)
{
  res=data.frame(matrix(NA,nrow=length(tumors),ncol=2))
  colnames(res)=c("purity","ploidy")
  rownames(res)=tumors
  for (i in 1:length(tumors))
  {
    load(paste0(absolutedir,"/",tumors[i],"/",tumors[i],".ABSOLUTE.RData"))
    tmp=as.data.frame(seg.dat$mode.res$mode.tab)
    tmp1=tmp[tmp$alpha>0.7 & tmp$alpha<0.9 & tmp$tau <2.5 & tmp$tau >1.5,]
    if (nrow(tmp1)>0)
    {
      res$purity[i]=tmp1$alpha[1]
      res$ploidy[i]=tmp1$tau[1]
    }else
    {
      tmp1=tmp[tmp$alpha>0.6 & tmp$alpha<0.9 & tmp$tau <3 & tmp$tau >1,]
      if (nrow(tmp1)>0)
      {
        idx=which.min(abs(tmp1$tau-2))
        res$purity[i]=tmp1$alpha[idx]
        res$ploidy[i]=tmp1$tau[idx]
      }
    }
  }
  return(res)
}

absoluteres=readabsolute(absolutedir,tumors=grp2tumors)
write.table(absoluteres,file="absolute.txt",sep="\t",quote=F)
