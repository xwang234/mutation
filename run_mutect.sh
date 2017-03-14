#! /usr/bin/env bash

data=dulak
echo $data
sleep 10

if [[ $data == "dulak" ]]
then
  wgsnormals=(SRR1002719 SRR999433 SRR999687 SRR1000378 SRR1002786 SRR999559 SRR1001730 SRR10018461 SRR1002703 SRR999599 SRR1002792 SRR1001839 SRR9994281 SRR1002710 SRR9995631 SRR1021476)
  wgstumors=(SRR1001842 SRR1002713 SRR999423 SRR1001466 SRR1002670 SRR1001823 SRR999489 SRR1002343 SRR1002722 SRR1002656 SRR1002929 SRR999438 SRR1001915 SRR999594 SRR1001868 SRR1001635)
  bamdir=/fh/scratch/delete30/dai_j/dulak
  outputdir=/fh/scratch/delete30/dai_j/mutect3
  #mutect3: qmin=10,mapmin=10,q1
  for ((i=0;i<${#wgstumors[@]};i++))
  do
    normalbam=$bamdir/${wgsnormals[$i]}.dedup.realigned.recal.bam
    tumorbam=$bamdir/${wgstumors[$i]}.dedup.realigned.recal.bam
    tumorname=${wgstumors[$i]}
    echo sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/mutect.sh $normalbam $tumorbam $tumorname $outputdir
    sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/mutect.sh $normalbam $tumorbam $tumorname $outputdir
    sleep 1
  done

fi

if [[ $data == "henan" ]]
then
  wgsnormals=(4A 12A 14A 16A 18A 26A 30A 34A 38A 42A)
  wgstumors=(3A 11A 13A 15A 17A 25A 29A 33A 37A 41A)
  #bamdir=/fh/scratch/delete30/dai_j/henan
  bamdir=/fh/scratch/delete30/dai_j/henan/newbam
  #outputdir=/fh/scratch/delete30/dai_j/henan/mutect3
  outputdir=/fh/scratch/delete30/dai_j/henan/mutectnewq1
  #mutect3: qmin=10,mapmin=10,q1
  qmin=10 
  mapmin=10

  #for ((i=0;i<${#wgstumors[@]};i++))
  for ((i=0;i<1;i++))
  do
    normalbam=$bamdir/${wgsnormals[$i]}.merged.deduprealigned.bam
    tumorbam=$bamdir/${wgstumors[$i]}.merged.deduprealigned.bam
    tumorname=${wgstumors[$i]}
    echo sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/mutect.sh $normalbam $tumorbam $tumorname $outputdir $qmin $mapmin
    sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/mutect.sh $normalbam $tumorbam $tumorname $outputdir $qmin $mapmin
    sleep 1
  done
fi
