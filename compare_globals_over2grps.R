#!/usr/bin/env Rscript
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

includedtypes2=c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
                 "Splice_Site","Translation_Start_Site","Nonstop_Mutation","Targeted_Region","De_novo_Start_InFrame","De_novo_Start_OutOfFrame")


grp1normals=c('SRR1002719','SRR999433','SRR999687','SRR1000378','SRR1002786','SRR999559','SRR1001730','SRR10018461','SRR1002703','SRR999599','SRR1002792','SRR1001839','SRR9994281','SRR1002710','SRR9995631')
grp1tumors=c('SRR1001842','SRR1002713','SRR999423','SRR1001466','SRR1002670','SRR1001823','SRR999489','SRR1002343','SRR1002722','SRR1002656','SRR1002929','SRR999438','SRR1001915','SRR999594','SRR1001868')

#grp2normals=c('2A','4A','6A','8A','10A','12A')
#grp2tumors=c('1A','3A','5A','7A','9A','11A')
grp2normals=paste0((c(3,11,13,15,17,25,29,33,37,41)+1),"A")
grp2tumors=paste0(c(3,11,13,15,17,25,29,33,37,41),"A")

grp1dir="/fh/scratch/delete30/dai_j/varscan"
grp2dir="/fh/scratch/delete30/dai_j/henan/varscan"
#compute transrate and lego
computeall=function(grpdir,grptumors)
{
  transcount=transcount_snv=transratemat=transrate_snv_mat=data.frame(matrix(NA,nrow=7,ncol=length(grptumors)))
  rownames(transcount)=rownames(transcount_snv)=rownames(transratemat)=rownames(transrate_snv_mat)=c('c2t','c2a','c2g','a2g','a2c','a2t','aa2c')
  colnames(transcount)=colnames(transcount_snv)=colnames(transratemat)=colnames(transrate_snv_mat)=grptumors
  legomat=lego_snv_mat=data.frame(matrix(NA,nrow=length(grptumors),ncol=96))
  rownames(legomat)=rownames(lego_snv_mat)=grptumors
  nummutation=nummutation_snv=data.frame(matrix(0,nrow=length(grptumors),ncol=1))
  rownames(nummutation)=rownames(nummutation_snv)=grptumors
  
  for (i in 1:length(grptumors))
  {
    print(grptumors[i])
    maffile=paste0(grpdir,"/",grptumors[i],".snp.Somatic.hc.annotated")
    if (! file.exists(maffile))
    {
      maffile=paste0(grpdir,"/",grptumors[i],".somatic.snp.Somatic.annotated")
    }
    maftable=read.table(maffile,header=T,sep="\t",quote="")
    maftable$Variant_Classification=as.character(maftable$Variant_Classification)
    colchr=5
    colstart=6
    colref=11
    colalt=13
    result=countmutations_fromtable(maftable,colchr,colstart,colref,colalt)
    transratemat[,i]=result$transrate
    transcount[,i]=result$transcount
    legomat[i,]=result$lego
    colnames(legomat)=colnames(lego_snv_mat)=names(result$lego)
    nummutation[i,1]=nrow(maftable)
    maftable_snv=maftable[maftable$Variant_Classification %in% includedtypes2,]
    result_snv=countmutations_fromtable(maftable_snv,colchr,colstart,colref,colalt)
    transrate_snv_mat[,i]=result_snv$transrate
    transcount_snv[,i]=result_snv$transcount
    lego_snv_mat[i,]=result_snv$lego
    nummutation_snv[i,1]=nrow(maftable_snv)
  }
  results=list(transratemat=transratemat,transcount=transcount,legomat=legomat,nummutation=nummutation,transrate_snv_mat=transrate_snv_mat,transcount_snv=transcount_snv,lego_snv_mat=lego_snv_mat,nummutation_snv=nummutation_snv)
}

dulak_res=computeall(grp1dir,grp1tumors)
henan_res=computeall(grp2dir,grp2tumors)
#save(dulak_res,henan_res,file="tmp.RData")

# #table to keep SRR->ESO naming transversion
# dulaktable=read.table(file="/fh/fast/dai_j/CancerGenomics/EAC/Dulak_fileinfo.txt",header=T)
# esotranstable=data.frame(matrix(NA,ncol=4,nrow=nrow(dulaktable)))
# dulaktable[,5]=as.character(dulaktable[,5])
# colnames(esotranstable)=c("srr","eso","type","wgs")
# esotranstable[,1]=as.character(dulaktable[,3])
# esotranstable[,4]=dulaktable[,9] #wgs:1,wes:0
# for (i in 1:nrow(dulaktable))
# {
#  tmp=unlist(strsplit(dulaktable[i,5],'-'))
#  esotranstable[i,2]=paste0(tmp[[1]],"-",tmp[[2]])
#  esotranstable[i,3]=tmp[[3]]
# }
# 
# dulakstable2=read.table(file="/fh/fast/dai_j/CancerGenomics/EAC/Dulak_Stable2.txt",header=T)
# srridx=rep(0,length(grp1tumors))
# for (i in 1:length(grp1tumors))
# {
#   tmp=which(esotranstable[,1]==grp1tumors[i])
#   tmp1=esotranstable[tmp,2]
#   srridx[i]=which(dulakstable2[,1]==tmp1)
# }
# #transfer into grp1's order
# srrdulakstable2=dulakstable2[srridx,]
# 
# dulakstable6=read.table(file="/fh/fast/dai_j/CancerGenomics/EAC/Dulak_Stable6.txt",header=T)
# srridx1=rep(0,length(grp1tumors))
# for (i in 1:length(grp1tumors))
# {
#   tmp=which(esotranstable[,1]==grp1tumors[i])
#   tmp1=esotranstable[tmp,2]
#   srridx1[i]=which(dulakstable6[,1]==tmp1)
# }
# #transfer into grp1's order
# srrdulakstable6=dulakstable6[srridx1,]


computep=function(x,y)
{
  
  print(mean(x))
  print(mean(y))
  p=2*(1-pnorm(abs((mean(x)-mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))))
  res=list(min1=min(x,na.rm=T),max1=max(x,na.rm=T),min2=min(y,na.rm=T),max2=max(y,na.rm=T),
           mean1=mean(x),mean2=mean(y),median1=median(x),median2=median(y),pvalue=p)
  print(p)
  return(res)
}

library("beeswarm")
library(GenomicRanges)
#the overall table, col1:dataset,col2:samplename
restable=data.frame(matrix(NA,nrow=length(grp1tumors)+length(grp2tumors),ncol=2))
colnames(restable)=c("dataset","sample")
restable[1:length(grp1tumors),1]="EA US"
restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),1]="EA CHINA"
restable[,2]=c(grp1tumors,grp2tumors)
color=c("green","blue")

#
ploidy=c(1.91,2,34,2.97,1.94,2.46,1.95,2.55,2.78,2.29,1.77,2.6,2.26,2.61,2.13,1.75,2.81,2.04,1.89,2.11,2.61,1.98,1.96,1.95,2.09)
contamination=c(0.32,0.25,0.06,0.29,0.26,0.12,0.16,0.21,0.22,0.18,0.07,0.09,0.2,0.19,0.25,0.22,0.32,0.25,0.35,0.13,0.15,0.21,0.3,0.14,0.17)
restable=cbind(restable,ploidy=ploidy)
outfig="ploidy_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(ploidy ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Ploidy')
beeswarm(ploidy ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

restable=cbind(restable,contamination=contamination)
outfig="contamination_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2c ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Contamination')
beeswarm(contamination ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#count num of mutations
#using .Somatic.annotated maf files
countallmutations=function(nummutation,tumors)
{
  res=data.frame(matrix(NA,ncol=2,nrow=length(tumors)))
  colnames(res)=c("tumor","num_mutations")
  for (numpair in 1:length(tumors))
  {
    res[numpair,1]=tumors[numpair]
    res[numpair,2]=nummutation[numpair,1]
  }
  return(res)
}

grp1mutations=countallmutations(dulak_res$nummutation,grp1tumors)
grp2mutations=countallmutations(henan_res$nummutation,grp2tumors)

restable=cbind(restable,c(grp1mutations[,2],grp2mutations[,2]))
colnames(restable)[ncol(restable)]="mutation"

res_mutation=computep(restable[1:length(grp1tumors),"mutation"],
                      restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"mutation"])

outfig="number_mutations_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(mutation ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        ylim=c(0,62000),
        main = 'Number of mutations')
beeswarm(mutation ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
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

grp1a2c=counttransrate(dulak_res$transratemat,grp1tumors)
grp2a2c=counttransrate(henan_res$transratemat,grp2tumors)
#add mutations
restable=cbind(restable,c(grp1a2c[,2],grp2a2c[,2]))
colnames(restable)[ncol(restable)]="a2c"
res_a2c=computep(restable[1:length(grp1tumors),"a2c"],
                 restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"a2c"])
outfig="a2c_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2c ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of A>C mutations')
beeswarm(a2c ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1aa2ac=counttransrate(dulak_res$transratemat,grp1tumors,opt="aa2c")
grp2aa2ac=counttransrate(henan_res$transratemat,grp2tumors,opt="aa2c")
#add mutations
restable=cbind(restable,c(grp1aa2ac[,2],grp2aa2ac[,2]))
colnames(restable)[ncol(restable)]="aa2ac"
res_aa2ac=computep(restable[1:length(grp1tumors),"aa2ac"],
                   restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"aa2ac"])
outfig="aa2ac_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(aa2ac ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Percentage of A>C at AA sites (%)')
beeswarm(aa2ac ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1a2g=counttransrate(dulak_res$transratemat,grp1tumors,opt="a2g")
grp2a2g=counttransrate(henan_res$transratemat,grp2tumors,opt="a2g")
#add mutations
restable=cbind(restable,c(grp1a2g[,2],grp2a2g[,2]))
colnames(restable)[ncol(restable)]="a2g"
res_a2g=computep(restable[1:length(grp1tumors),"a2g"],
                 restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"a2g"])
outfig="a2g_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2g ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of A>G mutations')
beeswarm(a2g ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1a2t=counttransrate(dulak_res$transratemat,grp1tumors,opt="a2t")
grp2a2t=counttransrate(henan_res$transratemat,grp2tumors,opt="a2t")
#add mutations
restable=cbind(restable,c(grp1a2t[,2],grp2a2t[,2]))
colnames(restable)[ncol(restable)]="a2t"
res_a2t=computep(restable[1:length(grp1tumors),"a2t"],
                 restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"a2t"])
outfig="a2t_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2t ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of A>T mutations')
beeswarm(a2t ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1c2a=counttransrate(dulak_res$transratemat,grp1tumors,opt="c2a")
grp2c2a=counttransrate(henan_res$transratemat,grp2tumors,opt="c2a")
#add mutations
restable=cbind(restable,c(grp1c2a[,2],grp2c2a[,2]))
colnames(restable)[ncol(restable)]="c2a"
res_c2a=computep(restable[1:length(grp1tumors),"c2a"],
                 restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"c2a"])
outfig="c2a_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2a ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of C>A mutations')
beeswarm(c2a ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1c2g=counttransrate(dulak_res$transratemat,grp1tumors,opt="c2g")
grp2c2g=counttransrate(henan_res$transratemat,grp2tumors,opt="c2g")
#add mutations
restable=cbind(restable,c(grp1c2g[,2],grp2c2g[,2]))
colnames(restable)[ncol(restable)]="c2g"
res_c2g=computep(restable[1:length(grp1tumors),"c2g"],
                 restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"c2g"])
outfig="c2g_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2g ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of C>G mutations')
beeswarm(c2g ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1c2t=counttransrate(dulak_res$transratemat,grp1tumors,opt="c2t")
grp2c2t=counttransrate(henan_res$transratemat,grp2tumors,opt="c2t")
#add mutations
restable=cbind(restable,c(grp1c2t[,2],grp2c2t[,2]))
colnames(restable)[ncol(restable)]="c2t"
res_c2t=computep(restable[1:length(grp1tumors),"c2t"],
                 restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"c2t"])
outfig="c2t_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2t ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of C>T mutations')
beeswarm(c2t ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#SNV mutations---------------------------------------------
restable1=data.frame(matrix(NA,nrow=length(grp1tumors)+length(grp2tumors),ncol=2))
colnames(restable1)=c("dataset","sample")
restable1[1:length(grp1tumors),1]="EA US"
restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),1]="EA CHINA"
restable1[,2]=c(grp1tumors,grp2tumors)

grp1mutations_snv=countallmutations(dulak_res$nummutation_snv,grp1tumors)
grp2mutations_snv=countallmutations(henan_res$nummutation_snv,grp2tumors)

restable1=cbind(restable1,c(grp1mutations_snv[,2],grp2mutations_snv[,2]))
colnames(restable1)[ncol(restable1)]="mutation"

res_mutation_snv=computep(restable1[1:length(grp1tumors),"mutation"],
                      restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"mutation"])

outfig="number_mutations_dulakhenan_snv.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(mutation ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        ylim=c(0,300),
        main = 'Number of mutations')
beeswarm(mutation ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1a2c_snv=counttransrate(dulak_res$transrate_snv_mat,grp1tumors)
grp2a2c_snv=counttransrate(henan_res$transrate_snv_mat,grp2tumors)
#add mutations
restable1=cbind(restable1,c(grp1a2c_snv[,2],grp2a2c_snv[,2]))
colnames(restable1)[ncol(restable1)]="a2c"
res_a2c_snv=computep(restable1[1:length(grp1tumors),"a2c"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"a2c"])
outfig="a2c_dulakhenan_snv.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2c ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of A>C mutations')
beeswarm(a2c ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1aa2ac_snv=counttransrate(dulak_res$transrate_snv_mat,grp1tumors,opt="aa2c")
grp2aa2ac_snv=counttransrate(henan_res$transrate_snv_mat,grp2tumors,opt="aa2c")
#add mutations
restable1=cbind(restable1,c(grp1aa2ac_snv[,2],grp2aa2ac_snv[,2]))
colnames(restable1)[ncol(restable1)]="aa2ac"
res_aa2ac_snv=computep(restable1[1:length(grp1tumors),"aa2ac"],
                   restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"aa2ac"])
outfig="aa2ac_dulakhenan_snv.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(aa2ac ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Percentage of A>C at AA sites (%)')
beeswarm(aa2ac ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1a2g_snv=counttransrate(dulak_res$transrate_snv_mat,grp1tumors,opt="a2g")
grp2a2g_snv=counttransrate(henan_res$transrate_snv_mat,grp2tumors,opt="a2g")
#add mutations
restable1=cbind(restable1,c(grp1a2g_snv[,2],grp2a2g_snv[,2]))
colnames(restable1)[ncol(restable1)]="a2g"
res_a2g_snv=computep(restable1[1:length(grp1tumors),"a2g"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"a2g"])
outfig="a2g_dulakhenan_snv.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2g ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of A>G mutations')
beeswarm(a2g ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1a2t_snv=counttransrate(dulak_res$transrate_snv_mat,grp1tumors,opt="a2t")
grp2a2t_snv=counttransrate(henan_res$transrate_snv_mat,grp2tumors,opt="a2t")
#add mutations
restable1=cbind(restable1,c(grp1a2t_snv[,2],grp2a2t_snv[,2]))
colnames(restable1)[ncol(restable1)]="a2t"
res_a2t_snv=computep(restable1[1:length(grp1tumors),"a2t"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"a2t"])
outfig="a2t_dulakhenan_snv.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2t ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of A>T mutations')
beeswarm(a2t ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1c2a_snv=counttransrate(dulak_res$transrate_snv_mat,grp1tumors,opt="c2a")
grp2c2a_snv=counttransrate(henan_res$transrate_snv_mat,grp2tumors,opt="c2a")
#add mutations
restable1=cbind(restable1,c(grp1c2a_snv[,2],grp2c2a_snv[,2]))
colnames(restable1)[ncol(restable1)]="c2a"
res_c2a_snv=computep(restable1[1:length(grp1tumors),"c2a"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"c2a"])
outfig="c2a_dulakhenan_snv.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2a ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of C>A mutations')
beeswarm(c2a ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1c2g_snv=counttransrate(dulak_res$transrate_snv_mat,grp1tumors,opt="c2g")
grp2c2g_snv=counttransrate(henan_res$transrate_snv_mat,grp2tumors,opt="c2g")
#add mutations
restable1=cbind(restable1,c(grp1c2g_snv[,2],grp2c2g_snv[,2]))
colnames(restable1)[ncol(restable1)]="c2g"
res_c2g_snv=computep(restable1[1:length(grp1tumors),"c2g"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"c2g"])
outfig="c2g_dulakhenan_snv.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2g ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of C>G mutations')
beeswarm(c2g ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

grp1c2t_snv=counttransrate(dulak_res$transrate_snv_mat,grp1tumors,opt="c2t")
grp2c2t_snv=counttransrate(henan_res$transrate_snv_mat,grp2tumors,opt="c2t")
#add mutations
restable1=cbind(restable1,c(grp1c2t_snv[,2],grp2c2t_snv[,2]))
colnames(restable1)[ncol(restable1)]="c2t"
res_c2t_snv=computep(restable1[1:length(grp1tumors),"c2t"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"c2t"])
outfig="c2t_dulakhenan_snv.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2t ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of C>T mutations')
beeswarm(c2t ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()




















#compute aa2ac mutation frequency
computeaa2acrate=function(legofiles,tumors)
{
  trinucleotidetable=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/trinucleotide_count.txt",header=T,sep="\t")
  #trinucleotidetable=read.table("/fh/scratch/delete30/dai_j/henan/mutect/trinucleotide_count.txt",header=T,sep="\t")
  #idx=which(grepl("^A",trinucleotidetable[,1],perl=T))
  #idx=which(grepl("^AA",trinucleotidetable[,1],perl=T))
  #len=sum(trinucleotidetable[idx,2]) #length of ^AXX
  idx1=which(grepl("^AT",trinucleotidetable[,1],perl=T))
  idx2=which(grepl("^AA",trinucleotidetable[,1],perl=T))
  len=sum(trinucleotidetable[idx1,2])+sum(trinucleotidetable[idx2,2])
  res=data.frame(matrix(NA,nrow=length(tumors),ncol=2))
  colnames(res)=c("sample",'aa2acfreq')
  res[,1]=tumors
  for (numpair in 1:length(tumors))
  {
    tumor=tumors[numpair]
    rawdata=read.table(file=legofiles[i],header=T)
    #idx1=which(grepl("T.G.A.C_",colnames(rawdata),perl=T))
    idx1=which(grepl("T.G.A.C_Ax",colnames(rawdata),perl=T))
    num=sum(rawdata[idx1])
    res[numpair,2]=num/len*10^6
  }
  return(res)
}


#count aa>ac mutation rate on coding region
countaa2ac_coding=function(transratefiles,tumors)
{
  res=data.frame(matrix(NA,ncol=2,nrow=length(tumors)))
  colnames(res)=c("tumor","mutationrate")
  for (numpair in 1:length(tumors))
  {
    tumor=tumors[numpair]
    tmp=read.table(transratefiles[i],header = T)
    res[numpair,1]=tumor
    res[numpair,2]=tmp[7,1]*100
  }
  return(res)
}

#compute aa2ac mutation frequency on coding region
computeaa2acrate_coding=function(mutectdir,tumors)
{
  trinucleotidetable=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/trinucleotide_count.txt",header=T,sep="\t")
  #trinucleotidetable=read.table("/fh/scratch/delete30/dai_j/henan/mutect/trinucleotide_count.txt",header=T,sep="\t")
  idx=which(grepl("^A",trinucleotidetable[,1],perl=T))
  #idx=which(grepl("^AA",trinucleotidetable[,1],perl=T))
  len=sum(trinucleotidetable[idx,2])/100 #length of ^AXX, 100 divided to approximate length in coding
  
  res=data.frame(matrix(NA,nrow=length(tumors),ncol=2))
  colnames(res)=c("sample",'aa2acfreq')
  res[,1]=tumors
  for (numpair in 1:length(tumors))
  {
    tumor=tumors[numpair]
    rawdata=read.table(file=paste0(mutectdir,"/",tumor,".exonic.lego.txt"),header=T)
    #idx1=which(grepl("T.G.A.C_",colnames(rawdata),perl=T))
    idx1=which(grepl("T.G.A.C_Ax",colnames(rawdata),perl=T))
    num=sum(rawdata[idx1])
    res[numpair,2]=num/len*10^6
  }
  return(res)
}

#count a>c mutation rate
counta2c=function(transratefiles,tumors,opt="a2c")
{
  res=data.frame(matrix(NA,ncol=2,nrow=length(tumors)))
  colnames(res)=c("tumor","mutationrate")
  for (numpair in 1:length(tumors))
  {
    tumor=tumors[numpair]
    tmp=read.table(transratefiles[i],header = T)
    res[numpair,1]=tumor
    idx=which(rownames(tmp)==opt)
    res[numpair,2]=tmp[idx,1]
  }
  return(res)
}

#cnv length count
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





#AA>AC mutation rates
grp1aa2acrate=computeaa2acrate(mutectdir1,grp1tumors)
grp2aa2acrate=computeaa2acrate(mutectdir2,grp2tumors)
#add mutations
restable=cbind(restable,c(grp1aa2acrate[,2],grp2aa2acrate[,2]))
colnames(restable)[5]="aa2acrate"
res_aa2acrate=computep(restable[1:length(grp1tumors),"aa2acrate"],
                   restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"aa2acrate"])
#outfig="aa2acrate_dulakhenan.png"
outfig="aa2acrate_dulakhenan2.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(aa2acrate ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=1.5,
        #ylim=c(2000,52000),
        main = 'Frequency of A>C at AA sites (per MB)')
beeswarm(aa2acrate ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()





#total length of CNV,lengths were generated by compute_cnvlength.R
cnvlenfile1="/fh/scratch/delete30/dai_j/freec/dulakcnvlengthcount.txt"
grp1cnvlen=countcnvlength(cnvlenfile1,grp1tumors)
cnvlenfile2="/fh/scratch/delete30/dai_j/henan/freec/henancnvlengthcount.txt"
grp2cnvlen=countcnvlength(cnvlenfile2,grp2tumors)
restable=cbind(restable,rbind(grp1cnvlen,grp2cnvlen)[,2:4])
res_lengain=computep(restable[1:length(grp1tumors),"lengain"],
                   restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"lengain"])
outfig="lengain_dulakhenan1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(lengain ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Length of Amplifications (Mb)')
beeswarm(lengain ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

res_lenloss=computep(restable[1:length(grp1tumors),"lenloss"],
                     restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"lenloss"])
outfig="lenloss_dulakhenan1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(lenloss ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Length of Deletions (Mb)')
beeswarm(lenloss ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

outfig="lenloh_dulakhenan1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
res_lenloh=computep(restable[1:length(grp1tumors),"lenloh"],
                     restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"lenloh"])
boxplot(lenloh ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Length of LOH (Mb)')
beeswarm(lenloh ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()
#combine lens
cnvtable=data.frame(matrix(NA,ncol=3,nrow=3*nrow(restable)))
colnames(cnvtable)=c("sample","status","len")
cnvtable[,1]=rep(restable[,2],3)
cnvtable[1:length(grp1tumors),2]="Dulak_gain"
start=length(grp1tumors)+1; end=length(grp1tumors)+length(grp2tumors)
cnvtable[start:end,2]="Our_gain"
start=length(grp1tumors)+length(grp2tumors)+1
end=2*length(grp1tumors)+length(grp2tumors)
cnvtable[start:end,2]="Dulak_loss"
start=2*length(grp1tumors)+length(grp2tumors)+1
end=2*length(grp1tumors)+2*length(grp2tumors)
cnvtable[start:end,2]="Our_loss"
start=2*length(grp1tumors)+2*length(grp2tumors)+1
end=3*length(grp1tumors)+2*length(grp2tumors)
cnvtable[start:end,2]="Dulak_LOH"
start=3*length(grp1tumors)+2*length(grp2tumors)+1
end=3*length(grp1tumors)+3*length(grp2tumors)
cnvtable[start:end,2]="Our_LOH"
cnvtable[,3]=c(restable[,"lengain"],restable[,"lenloss"],restable[,"lenloh"])
color1=rep(c("red","green","blue"),2)
outfig="len_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(len ~ status, data = cnvtable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        xaxt = "n",  xlab = "",
        #ylim=c(2000,52000),
        main = 'Length of CNA (Mb)')
# x axis with ticks but without labels
axis(1, labels = FALSE)
labels=c("Dulak_gain","Dulak_LOH","Dulak_loss","Our_gain","Our_LOH","Our_loss")
text(x =  seq_along(labels), y = par("usr")[3] - 1.5, srt = 45, adj = 1,
     labels = labels, xpd = TRUE, cex=1.2)
beeswarm(len ~ status, data = cnvtable,
         col = color1, pch = 16,add = TRUE, cex=2)
dev.off()

#work on cnv numbers
cnvnumfile1="/fh/scratch/delete30/dai_j/freec/dulakcnvnumcount.txt"
grp1cnvnum=countcnvlength(cnvnumfile1,grp1tumors)
cnvnumfile2="/fh/scratch/delete30/dai_j/henan/freec/henancnvnumcount.txt"
grp2cnvnum=countcnvlength(cnvnumfile2,grp2tumors)
tmp=rbind(grp1cnvnum,grp2cnvnum)[,2:4]*10^6

colnames(tmp)=c("numgain","numloss","numloh")
restable=cbind(restable,tmp)
res_numgain=computep(restable[1:length(grp1tumors),"numgain"],
                     restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"numgain"])
outfig="numgain_dulakhenan1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(numgain ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Number of Amplifications')
beeswarm(numgain ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

res_numloss=computep(restable[1:length(grp1tumors),"numloss"],
                     restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"numloss"])
outfig="numloss_dulakhenan1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(numloss ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Number of Deletions')
beeswarm(numloss ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

res_numloh=computep(restable[1:length(grp1tumors),"numloh"],
                    restable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp2tumors)),"numloh"])
outfig="numloh_dulakhenan1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(numloh ~ dataset, data = restable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Number of LOH')
beeswarm(numloh ~ dataset, data = restable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()
write.table(restable,file="compare_globals_dulakhenan.txt",col.names = T,row.names = F,sep="\t",quote=F)

#add num to cnvtable
cnvtable=cbind(cnvtable,c(restable[,"numgain"],restable[,"numloss"],restable[,"numloh"]))
colnames(cnvtable)[4]="num"
outfig="num_dulakhenan.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(num ~ status, data = cnvtable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        xaxt = "n",  xlab = "",
        #ylim=c(2000,52000),
        main = 'Number of CNA')
# x axis with ticks but without labels
axis(1, labels = FALSE)
labels=c("Dulak_gain","Dulak_LOH","Dulak_loss","Our_gain","Our_LOH","Our_loss")
text(x =  seq_along(labels), y = par("usr")[3] - 1.5, srt = 45, adj = 1,
     labels = labels, xpd = TRUE, cex=1.2)
beeswarm(num ~ status, data = cnvtable,
         col = color1, pch = 16,add = TRUE, cex=2)
dev.off()



#only work on dulak's data:
#compare the mutation number with dulak's study
#compare a2c rate
duplicaterestable=data.frame(matrix(NA,nrow=length(grp1tumors)*2,ncol=4))
colnames(duplicaterestable)=c("sample","dataset","mutation","a2c","a2g","c2t","c2a","a2t","c2g")
duplicaterestable[,1]=rep(grp1tumors,2)
duplicaterestable[,2]=c(rep("Dulak",15),rep("Our",15))
duplicaterestable[,3]=c(srrdulakstable2[,4],grp1mutations[,2])
duplicaterestable[,4]=c(srrdulakstable6[,3],grp1a2c[,2])
#a2g
duplicaterestable[,5]=c(srrdulakstable6[,4],grp1a2g[,2])
#c2t
duplicaterestable[,6]=c(srrdulakstable6[,5],grp1c2t[,2])
#c2a
duplicaterestable[,7]=c(srrdulakstable6[,6],grp1c2a[,2])
#a2t
duplicaterestable[,8]=c(srrdulakstable6[,7],grp1a2t[,2])
#c2g
duplicaterestable[,9]=c(srrdulakstable6[,8],grp1c2g[,2])

#add aa2ac rate
#From Dulak supplimentary table 11
tmp=c(42.7,78.2,9.96,14.41,18.67,28.9,19.16,11.73,23.5,12.11,69.5,17.45,3.62,4.3,11.35)
duplicaterestable[,10]=c(tmp,grp1aa2acrate[,2])
colnames(duplicaterestable)[10]="aa2acrate"
write.table(duplicaterestable,file="/fh/scratch/delete30/dai_j/henan/duplicaterestable.txt",col.names = T,row.names = F,sep="\t",quote=F)

res_dulak_mutation=computep(duplicaterestable[1:length(grp1tumors),"mutation"],
                            duplicaterestable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp1tumors)),"mutation"])
outfig="dulakdata_mutations_2study1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(mutation ~ dataset, data = duplicaterestable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Number of mutations')
beeswarm(mutation ~ dataset, data = duplicaterestable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

res_dulak_a2c=computep(duplicaterestable[1:length(grp1tumors),"a2c"],
                            duplicaterestable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp1tumors)),"a2c"])
outfig="dulakdata_a2c_2study1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2c ~ dataset, data = duplicaterestable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'proportion of A>C mutation')
beeswarm(a2c ~ dataset, data = duplicaterestable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

res_dulak_aa2acrate=computep(duplicaterestable[1:length(grp1tumors),"aa2acrate"],
                       duplicaterestable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp1tumors)),"aa2acrate"])
outfig="dulakdata_aa2acrate_2study1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(aa2acrate ~ dataset, data = duplicaterestable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'AA>AC frequency (per Mb)')
beeswarm(aa2acrate ~ dataset, data = duplicaterestable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#a2g
res_dulak_a2g=computep(duplicaterestable[1:length(grp1tumors),"a2g"],
                       duplicaterestable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp1tumors)),"a2g"])
outfig="dulakdata_a2g_2study1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2g ~ dataset, data = duplicaterestable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'proportion of A>G mutation')
beeswarm(a2g ~ dataset, data = duplicaterestable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#c2t
res_dulak_c2t=computep(duplicaterestable[1:length(grp1tumors),"c2t"],
                       duplicaterestable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp1tumors)),"c2t"])
outfig="dulakdata_c2t_2study1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2t ~ dataset, data = duplicaterestable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'proportion of C>T mutation')
beeswarm(c2t ~ dataset, data = duplicaterestable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#c2a
res_dulak_c2a=computep(duplicaterestable[1:length(grp1tumors),"c2a"],
                       duplicaterestable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp1tumors)),"c2a"])
outfig="dulakdata_c2a_2study1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2a ~ dataset, data = duplicaterestable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'proportion of C>A mutation')
beeswarm(c2a ~ dataset, data = duplicaterestable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#a2t
res_dulak_a2t=computep(duplicaterestable[1:length(grp1tumors),"a2t"],
                       duplicaterestable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp1tumors)),"a2t"])
outfig="dulakdata_a2t_2study1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2t ~ dataset, data = duplicaterestable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'proportion of A>T mutation')
beeswarm(a2t ~ dataset, data = duplicaterestable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#c2g
res_dulak_c2g=computep(duplicaterestable[1:length(grp1tumors),"c2g"],
                       duplicaterestable[(length(grp1tumors)+1):(length(grp1tumors)+length(grp1tumors)),"c2g"])
outfig="dulakdata_c2g_2study1.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2g ~ dataset, data = duplicaterestable, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'proportion of C>G mutation')
beeswarm(c2g ~ dataset, data = duplicaterestable, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()


#work on AJHG paper which has 14 WGS, comparing them with henan data
oldres=read.table(file="compare_globals_dulakhenan.txt",sep="\t",header=T)
grp3tumors=c("1N01-VS-1T01","1N02-VS-1T02","1N03-VS-1T03","1N04-VS-1T04","1N05-VS-1T05","3N01-VS-3T01","3N02-VS-3T02","3N03-VS-3T03",
"3N04-VS-3T04","3N05-VS-3T05","3N06-VS-3T06","3N07-VS-3T07","3N08-VS-3T08","3N09-VS-3T09")
#from table MMC3
nummuations=c(16687,4382,13264,7559,13563,14119,15286,9335,7504,11055,5470,7204,12122,17319)
grp3muations=data.frame(matrix(NA,ncol=2,nrow=length(grp3tumors)))
colnames(grp3muations)=c("tumor","num_mutations")
grp3muations[,1]=grp3tumors
grp3muations[,2]=nummuations
#The AJHG data only have mutation in coding region, so we need to compute mutations of henan data in coding region
grp2tumors=c('1A','3A','5A','7A','9A','11A')
mutectdir2="/fh/scratch/delete30/dai_j/henan/mutect"
grp2aa2ac_coding=countaa2ac_coding(mutectdir2,grp2tumors)
grp2aa2acrate_coding=computeaa2acrate_coding(mutectdir2,grp2tumors)

#form dataframe contains 3 datasets
oldres$dataset=as.character(oldres$dataset)
threedata=oldres[,1:3]
threedata[which(threedata$dataset=="Dulak"),"dataset"]="US EA"
threedata[which(threedata$dataset=="Our"),"dataset"]="China EA"
adddata=data.frame(matrix(NA,ncol=ncol(threedata),nrow=length(grp3tumors)))
colnames(adddata)=colnames(threedata)
adddata[,1]="China ESCC"
adddata[,2]=grp3tumors
adddata[,3]=grp3muations[,2]
threedata=rbind(threedata,adddata)
color3=c("green","blue","red")
outfig="number_mutations_dulak_henan_escc.png"
png(outfig, width = 8, height = 4, units = 'in', res=300)
boxplot(mutation ~ dataset, data = threedata, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        ylim=c(2000,68000),
        main = 'Number of mutations')
beeswarm(mutation ~ dataset, data = threedata, 
         col = color3, pch = 16,add = TRUE, cex=2)
dev.off()

counta2c_varscan=function(data,opt="a2c")
{
  tmp1=unlist(strsplit(opt,"2"))
  ref=toupper(tmp1[1])
  alt=toupper(tmp1[2])
  res=nrow(data[data$ref==ref & data$var==alt,])/nrow(data)
}

#count mutation rate using varscan2---------------------------------------------
countmutations_varscan=function(varscandir,tumors)
{
  res=data.frame(matrix(NA,ncol=2,nrow=length(tumors)))
  colnames(res)=c("tumor","num_mutations")
  for (numpair in 1:length(tumors))
  {
    tumor=tumors[numpair]
    maffile=paste0(varscandir,'/',tumor,".somatic_keep.txt")
    tmp=read.table(maffile,header = T)
    res[numpair,1]=tumor
    res[numpair,2]=nrow(tmp)
  }
  return(res)
}
varscandir="/fh/scratch/delete30/dai_j/henan/varscan"
grp3normals=paste0((c(3,11,13,15,17,23,25,33,35,37,41)+1),"A")
grp3tumors=paste0(c(3,11,13,15,17,23,25,33,35,37,41),"A") #add 29 later
#the overall table, col1:dataset,col2:samplename
restable1=data.frame(matrix(NA,nrow=length(grp1tumors)+length(grp3tumors),ncol=2))
colnames(restable1)=c("dataset","sample")
restable1[1:length(grp1tumors),1]="Dulak"
restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp3tumors)),1]="Our"
restable1[,2]=c(grp1tumors,grp3tumors)
color=c("green","blue")
mutectdir1="/fh/scratch/delete30/dai_j/mutect1"
grp1mutations=countmutations(mutectdir1,grp1tumors)
mutectdir2="/fh/scratch/delete30/dai_j/henan/mutect"
grp3mutations=countmutations_varscan(varscandir,grp3tumors)
#add mutations
restable1=cbind(restable1,c(grp1mutations[,2],grp3mutations[,2]))
colnames(restable1)[3]="mutation"

res_mutation=computep(restable1[1:length(grp1tumors),"mutation"],
                      restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp3tumors)),"mutation"])
outfig="number_mutations_dulakhenan3.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(mutation ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        ylim=c(0,62000),
        main = 'Number of mutations')
beeswarm(mutation ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#AA>AC mutations
grp1aa2ac=countaa2ac(mutectdir1,grp1tumors)
grp3aa2ac=countaa2ac(mutectdir2,paste0(grp3tumors,".varscan_somatic"))
#add mutations
restable1=cbind(restable1,c(grp1aa2ac[,2],grp3aa2ac[,2]))
colnames(restable1)[4]="aa2ac"
res_aa2ac=computep(restable1[1:length(grp1tumors),"aa2ac"],
                   restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp3tumors)),"aa2ac"])
outfig="aa2ac_dulakhenan3.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(aa2ac ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Percentage of A>C at AA sites (%)')
beeswarm(aa2ac ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#AA>AC mutation rates
grp1aa2acrate=computeaa2acrate(mutectdir1,grp1tumors)
grp3aa2acrate=computeaa2acrate(mutectdir2,paste0(grp3tumors,".varscan_somatic"))
#add mutations
restable1=cbind(restable1,c(grp1aa2acrate[,2],grp3aa2acrate[,2]))
colnames(restable1)[5]="aa2acrate"
res_aa2acrate=computep(restable1[1:length(grp1tumors),"aa2acrate"],
                       restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp3tumors)),"aa2acrate"])
outfig="aa2acrate_dulakhenan3.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(aa2acrate ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=1.5,
        #ylim=c(2000,52000),
        main = 'Frequency of A>C at AA sites (per MB)')
beeswarm(aa2acrate ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()


#A>C mutations
grp1a2c=counta2c(mutectdir1,grp1tumors)
grp3a2c=counta2c(mutectdir2,paste0(grp3tumors,".varscan_somatic"))
#add mutations
restable1=cbind(restable1,c(grp1a2c[,2],grp3a2c[,2]))
colnames(restable1)[6]="a2c"
res_a2c=computep(restable1[1:length(grp1tumors),"a2c"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp3tumors)),"a2c"])
outfig="a2c_dulakhenan3.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2c ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of A>C mutations')
beeswarm(a2c ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#A>G mutations
grp1a2g=counta2c(mutectdir1,grp1tumors,opt="a2g")
grp3a2g=counta2c(mutectdir2,paste0(grp3tumors,".varscan_somatic"),opt="a2g")
#add mutations
restable1=cbind(restable1,c(grp1a2g[,2],grp3a2g[,2]))
colnames(restable1)[7]="a2g"
res_a2g=computep(restable1[1:length(grp1tumors),"a2g"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp3tumors)),"a2g"])

outfig="a2g_dulakhenan3.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2g ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of A>G mutations')
beeswarm(a2g ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#C>T mutations
grp1c2t=counta2c(mutectdir1,grp1tumors,opt="c2t")
grp3c2t=counta2c(mutectdir2,paste0(grp3tumors,".varscan_somatic"),opt="c2t")
#add mutations
restable1=cbind(restable1,c(grp1c2t[,2],grp3c2t[,2]))
colnames(restable1)[8]="c2t"
res_c2t=computep(restable1[1:length(grp1tumors),"c2t"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp3tumors)),"c2t"])
outfig="c2t_dulakhenan3.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2t ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of C>T mutations')
beeswarm(c2t ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#C>A mutations
grp1c2a=counta2c(mutectdir1,grp1tumors,opt="c2a")
grp3c2a=counta2c(mutectdir2,paste0(grp3tumors,".varscan_somatic"),opt="c2a")
#add mutations
restable1=cbind(restable1,c(grp1c2a[,2],grp3c2a[,2]))
colnames(restable1)[9]="c2a"
res_c2a=computep(restable1[1:length(grp1tumors),"c2a"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp3tumors)),"c2a"])
outfig="c2a_dulakhenan3.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2a ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of C>A mutations')
beeswarm(c2a ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#A>T mutations
grp1a2t=counta2c(mutectdir1,grp1tumors,opt="a2t")
grp3a2t=counta2c(mutectdir2,paste0(grp3tumors,".varscan_somatic"),opt="a2t")
#add mutations
restable1=cbind(restable1,c(grp1a2t[,2],grp3a2t[,2]))
colnames(restable1)[10]="a2t"
res_a2t=computep(restable1[1:length(grp1tumors),"a2t"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp3tumors)),"a2t"])
outfig="a2t_dulakhenan3.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(a2t ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of A>T mutations')
beeswarm(a2t ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()

#C>G mutations
grp1c2g=counta2c(mutectdir1,grp1tumors,opt="c2g")
grp3c2g=counta2c(mutectdir2,paste0(grp3tumors,".varscan_somatic"),opt="c2g")
#add mutations
restable1=cbind(restable1,c(grp1c2g[,2],grp3c2g[,2]))
colnames(restable1)[11]="c2g"
res_c2g=computep(restable1[1:length(grp1tumors),"c2g"],
                 restable1[(length(grp1tumors)+1):(length(grp1tumors)+length(grp3tumors)),"c2g"])
outfig="c2g_dulakhenan3.png"
png(outfig, width = 6, height = 4, units = 'in', res=300)
boxplot(c2g ~ dataset, data = restable1, 
        outline = FALSE,     ## avoid double-plotting outliers, if any
        cex.axis=2,
        cex.main=2,
        #ylim=c(2000,52000),
        main = 'Proportion of C>G mutations')
beeswarm(c2g ~ dataset, data = restable1, 
         col = color, pch = 16,add = TRUE, cex=2)
dev.off()