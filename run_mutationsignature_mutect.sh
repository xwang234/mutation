#! /usr/bin/env bash

data="golden4"
colchr=3
colstart=4
colref=6
colalt=7
colsamplename=27

if [[ $data == "henan" ]]
then
  mutectdir=/fh/scratch/delete30/dai_j/henan/mutect3
  pref=henan_mutect3
  tumors=()
  tumors=(3A 11A 13A 15A 17A 25A 29A 33A 37A 41A)
fi

if [[ $data == "dulak" ]]
then
  mutectdir=/fh/scratch/delete30/dai_j/mutect1
  pref=dulak_mutect1
  tumors=()
  tumors=(SRR1001842 SRR1002713 SRR999423 SRR1001466 SRR1002670 SRR1001823 SRR999489 SRR1002343 SRR1002722 SRR1002656 SRR1002929 SRR999438 SRR1001915 SRR999594 SRR1001868 SRR1001635)
fi

if [[ $data == "escc" ]]
then
  mutectdir=/fh/scratch/delete30/dai_j/escc/mutect2
  pref=escc_mutect2
  tumors=()
  tumors=(T1 T2 T3 T4 T5 T6 T8 T9 T10 T11 T12 T13 T14 T15 T16 T17 T18)
fi

if [[ $data == "golden4" ]]
then
  mutectdir=/fh/scratch/delete30/dai_j/henan/mutect2
  pref=henan_golden4
  tumors=()
  tumors=(3A 11A 13A 15A)
fi

echo $mutectdir
echo $pref
echo ${tumors[@]}
sleep 10

inputmfile=/fh/fast/dai_j/CancerGenomics/Tools/mutationsignature/input/${pref}_inputformutationsignature.m
#matfile:
inputFile=/fh/fast/dai_j/CancerGenomics/Tools/mutationsignature/input/${pref}_inputformutationsignature.mat
#allOutputFile
allOutputFile=/fh/fast/dai_j/CancerGenomics/Tools/mutationsignature/output/${pref}_allres.mat

declare -a maffiles
maffiles=()
for ((i=0;i<${#tumors[@]};i++))
do
  maffiles[$i]=$mutectdir/${tumors[$i]}.variant_function
done

echo "/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/form_input_for_mutationsignature.R $mutectdir $colchr $colstart $colref $colalt $pref ${tumors[@]} ${maffiles[@]}"

Rscript /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/form_input_for_mutationsignature.R $mutectdir $colchr $colstart $colref $colalt $pref ${tumors[@]}  ${maffiles[@]}

cd /fh/fast/dai_j/CancerGenomics/Tools/mutationsignature/input/
matlab -nodesktop -nodisplay -nojvm -r "${pref}_inputformutationsignature;exit;"
cd /fh/fast/dai_j/CancerGenomics/Tools/mutationsignature
  #matlabpool only works with this version ?
module load matlab/R2014b
matlab -nodesktop -nodisplay -r "inputFile='$inputFile';allOutputFile='$allOutputFile';pref='$pref';run_mutationsignature2;exit;"

#check results
#matlab -nodesktop -nodisplay -r "load $allOutputFile;disp(reconstructionError);disp(stability);exit;"

#plot result
#number=3
#outputFileName=output/${pref}_${number}_signatures.mat
#matlab -nodesktop -nodisplay -r "allOutputFile='$allOutputFile';plotoutputmat;exit;"
cd $mutectdir

