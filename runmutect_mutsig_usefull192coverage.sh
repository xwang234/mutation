#!/usr/bin/env bash

data="henan"
echo $data
sleep 10
if [[ $data == "henan" ]]
then
  outputpreface=mutect5_henan
  mutsigoutputfolder=/fh/scratch/delete30/dai_j/henan/mutsig
  wgstumors=(3A 11A 13A 15A 17A 25A 29A 33A 37A 41A)
  mutectfolder=/fh/scratch/delete30/dai_j/henan/mutect5
fi
if [[ $data == "dulak" ]]
then
  outputpreface=mutect4_dulak
  mutsigoutputfolder=/fh/scratch/delete30/dai_j/mutsig
  wgstumors=(SRR1001842 SRR1002713 SRR999423 SRR1001466 SRR1002670 SRR1001823 SRR999489 SRR1002343 SRR1002722 SRR1002656 SRR1002929 SRR999438 SRR1001915 SRR999594 SRR1001868 SRR1001635)
  mutectfolder=/fh/scratch/delete30/dai_j/mutect4
fi
if [[ $data == "escc" ]]
then
  outputpreface=mutect2_escc
  mutsigoutputfolder=/fh/scratch/delete30/dai_j/escc/mutsig
  wgstumors=(T1 T2 T3 T4 T5 T6 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18)
  mutectfolder=/fh/scratch/delete30/dai_j/escc/mutect4
fi
cd $mutectfolder
declare -a mutationfiles
mutationfiles=()
for ((i=0;i<${#wgstumors[@]};i++))
do
  mutationfiles[$i]=$mutectfolder/${wgstumors[$i]}.Mutect_out_keep.txt
done
/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/mutect_mutsig_usefull192coverage.sh $outputpreface $mutectfolder $mutsigoutputfolder ${wgstumors[@]} ${mutationfiles[@]} 
