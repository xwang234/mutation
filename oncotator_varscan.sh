#!/usr/bin/env bash
#SBATCH -t 1-0
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

arguments=( "$@" )


echo ${arguments[@]}
((numarg=$#))

echo $numarg
let numsample=$numarg/2
echo $numsample

wgstumors=()
for ((i=0;i<$numsample;i++))
do
  wgstumors[$i]=${arguments[$i]}
done

mutationfiles=()
for ((i=$numsample;i<=$numarg;i++))
do
  mutationfiles[$i-$numsample]=${arguments[$i]}
done

echo ${wgstumors[@]}
echo ${mutationfiles[@]}

echo ${wgstumors[0]}
echo ${mutationfiles[0]}
sleep 1

echo "oncotator ...."
for ((i=0;i<${#wgstumors[@]};i++))
do
  echo ${wgstumors[$i]}
  if [[ ! -f ${mutationfiles[$i]} ]]
  then
    /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/combine_somatic_chr.sh ${wgstumors[$i]}
  fi
  /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/oncotator.sh ${mutationfiles[$i]}
done

sleep 10

for ((i=0;i<${#wgstumors[@]};i++))
do
  /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/process_varscan_oncotatorout.sh ${mutationfiles[$i]}.oncotator ${wgstumors[$i]}
done
