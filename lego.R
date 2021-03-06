## This script creates a "legoplot" similar to those produced by the Broad Institute
## The plot shows the relative abundance of each of the 6 possible mutations in the
## 16 sequence contexts

## Load packages
library(rgl)

#### START OF FUNCTIONS

## Functions modified from the "demo(hist3d)" examples in the rgl package:
# library(rgl)
# demo(hist3d)
## Note; would it have killed the original author to comment their code?

## Draws a single "column" or "stack".
## X and Y coordinates determine the area of the column
## The Z coordinate determines the height of the column
## We include "lit=FALSE" arguments to remove the nasty shiny surfaces caused by lighting
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
# Example:
# stackplot.3d(c(0,1),c(0,1),3,alpha=0.6)

## Calls stackplot.3d repeatedly to create a barplot
## z is the heights of the columns and must be an appropriately named vector
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
  filename=paste0(filename,".pdf")
  rgl.postscript(filename, fmt="pdf")
}

compute_trinucleotide_freq=function(counts,trinucleotidetable,numsample)
{
  res=counts
  for (i in 1:length(counts))
  {
    myname=names(counts)[i]
    tmp=unlist(strsplit(myname,"_"))
    tmp1=unlist(strsplit(tmp,".",fixed=T))
    tmp2=unlist(strsplit(tmp1[5],"x"))
    tri1=paste0(tmp2[1],tmp1[1],tmp2[2])
    tri2=paste0(tmp2[1],tmp1[3],tmp2[2])
    count1=trinucleotidetable[trinucleotidetable[,1]==tri1,2]
    count2=trinucleotidetable[trinucleotidetable[,1]==tri2,2]
    count=count1+count2
    res[i]=counts[i]/count/numsample*10^6
  }
  return(res)
}
# Example:
# context3d(counts)

#### END OF FUNCTIONS

## Read in example data and cast to an appropriate vector
#rawdata=read.table("snvspermegabase.txt",header=TRUE)
trinucleotidetable=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/trinucleotide_count.txt",header=T,sep="\t")

rawdata=read.table("1A.exonic.lego.txt",header=TRUE)
rawdata=read.table("all.exonic.lego.txt",header=TRUE)
rawdata=read.table("all.lego.txt",header=TRUE)
rawdata=read.table("/fh/scratch/delete30/dai_j/henan/varscan1/henan_varscan_lego.txt",header=T)
#henan data:
rawdata=read.table("/fh/scratch/delete30/dai_j/henan/varscan2/henan_varscan_lego.txt",header=T)
numsample=10
#golden4:
rawdata=read.table("/fh/scratch/delete30/dai_j/henan/varscan2/henan_golden4_varscan_lego.txt",header=T)
numsample=4

rawdata=read.table("/fh/scratch/delete30/dai_j/henan/mutect3/henan_mutect3_lego.txt",header=T)
numsample=10

#escc data:
rawdata=read.table("/fh/scratch/delete30/dai_j/escc/varscan2/escc_varscan_lego.txt",header=T)
numsample=17

#dulak data:
rawdata=read.table("/fh/scratch/delete30/dai_j/varscan2/dulak_varscan_lego.txt",header=T)
rawdata=read.table("/fh/scratch/delete30/dai_j/mutect1/dulak.lego.txt",header=T)
numsample=16

rawdata=read.table("/fh/scratch/delete30/dai_j/mutect1/dulak_mutect1_lego.txt",header=T)

counts0=as.numeric(rawdata)
names(counts0)=colnames(rawdata)

counts=compute_trinucleotide_freq(counts0,trinucleotidetable,numsample)
## Example plots

context3d(z=counts,filename="dulak.all.exonic.lego" )
context3d(z=counts,filename="dulak.all.lego")
context3d(z=counts,filename="/fh/scratch/delete30/dai_j/henan/varscan1/henan_varscan_lego" )
#henan data
context3d(z=counts,filename="/fh/scratch/delete30/dai_j/henan/varscan2/henan_varscan_lego")
context3d(z=counts,filename="/fh/scratch/delete30/dai_j/henan/mutect3/henan_mutect3_lego")

#golden4
context3d(z=counts,filename="/fh/scratch/delete30/dai_j/henan/varscan2/henan_golden4_varscan_lego")

#escc data
context3d(z=counts,filename="/fh/scratch/delete30/dai_j/escc/varscan2/escc_varscan_lego")

#dulak data
context3d(z=counts,filename="/fh/scratch/delete30/dai_j/varscan2/dulak_varscan_lego")
context3d(z=counts,filename="/fh/scratch/delete30/dai_j/mutect1/dulak_lego")

context3d(z=counts,filename="/fh/scratch/delete30/dai_j/mutect1/dulak_mutect1_lego")


context3d(counts,alpha=0.4)

## Save your images to files if you wish
rgl.snapshot(filename="example.png")
rgl.postscript("graph.svg", fmt="svg")