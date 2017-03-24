#!/usr/bin/env bash

data="henan"
echo $data
sleep 10
if [[ $data == "henan" ]]
then
  #chr=16
  mutectdir=/fh/scratch/delete30/dai_j/henan/mutect3
  normals=(4A 12A 14A 16A 18A 26A 30A 34A 38A 42A)
  tumors=(3A 11A 13A 15A 17A 25A 29A 33A 37A 41A)
  for ((i=0;i<${#normals[@]};i++))
  do
    wigfile=$mutectdir/${tumors[$i]}.wig.txt
    bedfile=$mutectdir/${tumors[$i]}.wig.bed.txt
    echo /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/wig2bed.sh $wigfile $bedfile
    sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/wig2bed.sh $wigfile $bedfile
    sleep 1
 done
fi

if [[ $data == "dulak" ]]
then
  mutectdir=/fh/scratch/delete30/dai_j/mutect3
  normals=(SRR1002719 SRR999433 SRR999687 SRR1000378 SRR1002786 SRR999559 SRR1001730 SRR10018461 SRR1002703 SRR999599 SRR1002792 SRR1001839 SRR9994281 SRR1002710 SRR9995631 SRR1021476)
  tumors=(SRR1001842 SRR1002713 SRR999423 SRR1001466 SRR1002670 SRR1001823 SRR999489 SRR1002343 SRR1002722 SRR1002656 SRR1002929 SRR999438 SRR1001915 SRR999594 SRR1001868 SRR1001635)
  for ((i=0;i<${#normals[@]};i++))
  do
      wigfile=$mutectdir/${tumors[$i]}.wig.txt
      bedfile=$mutectdir/${tumors[$i]}.wig.bed.txt
      echo /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/wig2bed.sh $wigfile $bedfile
      sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/wig2bed.sh $wigfile $bedfile
      sleep 1
  done
fi

