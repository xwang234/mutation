#!/usr/bin/env bash
data="escc"
if [[ $data == "dulak" ]]
then
  wgsnormals=(SRR1002719 SRR999433 SRR999687 SRR1000378 SRR1002786 SRR999559 SRR1001730 SRR10018461 SRR1002703 SRR999599 SRR1002792
 SRR1001839 SRR9994281 SRR1002710 SRR9995631 SRR1021476)
  wgstumors=(SRR1001842 SRR1002713 SRR999423 SRR1001466 SRR1002670 SRR1001823 SRR999489 SRR1002343 SRR1002722 SRR1002656 SRR1002929
 SRR999438 SRR1001915 SRR999594 SRR1001868 SRR1001635)
  bamdir=/fh/scratch/delete30/dai_j/dulak
  outputdir=/fh/scratch/delete30/dai_j/mutect2
  if [ ! -d $outputdir ]; then mkdir $outputdir; fi
  declare -a normalfiles
  declare -a tumorfiles
  declare -a outputs
  for ((i=0;i<${#wgstumors[@]};i++))
  do
    normalfiles[$i]=$bamdir/${wgsnormals[$i]}.dedup.realigned.recal.bam
	tumorfiles[$i]=$bamdir/${wgstumors[$i]}.dedup.realigned.recal.bam
	outputs[$i]=$outputdir/${wgstumors[$i]}.Mutect2_out.txt
    echo "sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/mutect_GATK.sh ${normalfiles[$i]} ${tumorfiles[$i]} ${outputs[$i]}"
    sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/mutect_GATK.sh ${normalfiles[$i]} ${tumorfiles[$i]} ${outputs[$i]}
    sleep 1
  done
fi


if [[ $data == "escc" ]]
then
  id=(1 2 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18)
  bamdir=/fh/scratch/delete30/dai_j/escc
  outputdir=/fh/scratch/delete30/dai_j/escc/mutect2
  if [ ! -d $outputdir ]; then mkdir $outputdir; fi
  declare -a normalfiles
  declare -a tumorfiles
  declare -a outputs
  for ((i=0;i<${#id[@]};i++))
  do
    normalfiles[$i]=$bamdir/N${id[$i]}.merged.deduprealigned.bam
    tumorfiles[$i]=$bamdir/T${id[$i]}.merged.deduprealigned.bam
    outputs[$i]=$outputdir/T${id[$i]}.Mutect2_out.txt
    echo "sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/mutect_GATK.sh ${normalfiles[$i]} ${tumorfiles[$i]} ${outputs[$i]}"
    sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/mutect_GATK.sh ${normalfiles[$i]} ${tumorfiles[$i]} ${outputs[$i]}
    sleep 1
  done
fi  
