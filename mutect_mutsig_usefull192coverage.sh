#!/usr/bin/env bash
#SBATCH -t 1-0
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

arguments=( "$@" )

outputpreface=${arguments[0]}
mutectfolder=${arguments[1]}
mutsigoutputfolder=${arguments[2]}
echo ${arguments[@]}
((numarg=$#))

echo $numarg
let numsample=($numarg-3)/2
echo $numsample
declare -a wgstumors
for ((i=3;i<=$numsample+2;i++))
do
  wgstumors[(($i-3))]=${arguments[$i]}
done

declare -a mutationfiles
for ((i=$numsample+3;i<=$numarg-1;i++))
do
  mutationfiles[$i-$numsample-3]=${arguments[$i]}
done



echo $outputpreface
echo $mutectfolder

echo ${wgstumors[@]}
echo ${mutationfiles[@]}

echo ${wgstumors[0]}
echo ${mutationfiles[0]}

sleep 15

oncotator=/fh/fast/dai_j/CancerGenomics/Tools/oncotator1/oncotator-1.8.0.0/oncotator/bin/oncotator
oncotator_db=/fh/fast/dai_j/CancerGenomics/Tools/oncotator/oncotator_v1_ds_June112014
#mutsig configure data folder:
mutsigfolder=/fh/fast/dai_j/CancerGenomics/Tools/MutSigCV_1.4
#mutsigcmdfolder=/fh/fast/dai_j/CancerGenomics/Tools/MutSigCV_1.4
#the newest version
mutsigcmdfolder=/fh/fast/dai_j/CancerGenomics/Tools/MutSigCV_1.41
#coveragefile=./full192.coverage.txt
coveragefile=$mutsigfolder/exome_full192.coverage.txt
covariatesfile=$mutsigfolder/gene.covariates.txt
dictfile=$mutsigfolder/mutation_type_dictionary_file.txt
hg19file=$mutsigfolder/chr_files
#outputfile1=$mutsigoutputfolder/${outputpreface}_a.mutsiga.txt
#outputfile2=$mutsigoutputfolder/${outputpreface}_b.mutsiga.txt
outputfile=$mutsigoutputfolder/${outputpreface}.mutsiga.txt
#reference=/fh/fast/dai_j/CancerGenomics/Tools/MutSigCV_1.4/exome_full192.coverage.txt
maffile=$mutsigoutputfolder/${outputpreface}.mutectmaf

#generate maf file
echo "oncotator..."
for ((i=0;i<${#wgstumors[@]};i++))
do
  echo $i
  echo ${wgstumors[$i]}
  if [[ ! -s ${wgstumors[$i]}.mutectmaf.ancotator ]]
  then
    $oncotator -v --db-dir $oncotator_db ${mutationfiles[$i]} ${wgstumors[$i]}.mutectmaf.ancotator hg19
  fi
done


if [ -f $maffile ]; then rm $maffile; fi

for ((i=0;i<${#wgstumors[@]};i++))
do
#   /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/process_mutect_oncotatorout.sh ${mutationfiles[$i]}.oncotator ${wgstumors[$i]}
   cat ${wgstumors[$i]}.mutectmaf.ancotator >> $maffile
done

echo "run mutsig..."
sleep 10
cd $mutsigcmdfolder
module load matlab/R2013b
#matlab -nodesktop -nodisplay -nojvm -r "MutSigCV('$maffile','$coveragefile','$covariatesfile','$outputfile1');exit;"
matlab -nodesktop -nodisplay -nojvm -r "MutSigCV('$maffile','$coveragefile','$covariatesfile','$outputfile','$dictfile','$hg19file');exit;"
cd $mutsigoutputfolder

#MutSigCV('./nobarrett.maf','./exome_full192.coverage.txt','exampledata/gene.covariates.txt','./nobarrett_mutsigout.txt')
#MutSigCV('/mnt/rhodium/scratch/delete30/dai_j/test/mutsig/allinfo.reduced.maf','/mnt/rhodium/scratch/delete30/dai_j/test/mutsig/full192.coverage.txt','exampledata/gene.covariates.txt','./allinfo_mutsigout.txt')
