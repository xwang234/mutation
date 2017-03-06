#!/usr/bin/env bash
#SBATCH -t 1-0
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

arguments=( "$@" )

outputpreface=${arguments[0]}
varscanfolder=${arguments[1]}
mutsigoutputfolder=${arguments[2]}
echo ${arguments[@]}
((numarg=$#))

echo $numarg
let numsample=($numarg-3)/3
echo $numsample
declare -a wgstumors
for ((i=3;i<=$numsample+2;i++))
do
  wgstumors[(($i-3))]=${arguments[$i]}
done

declare -a mutationfiles
for ((i=$numsample+3;i<=$numarg-$numsample-1;i++))
do
  mutationfiles[$i-$numsample-3]=${arguments[$i]}
done

declare -a indelfiles
for ((i=$numarg-$numsample;i<=$numarg;i++))
do
  indelfiles[$i-$numarg+$numsample]=${arguments[$i]}
done

echo $outputpreface
echo $varscanfolder

echo ${wgstumors[@]}
echo ${mutationfiles[@]}
echo ${indelfiles[@]}

echo ${wgstumors[0]}
echo ${mutationfiles[0]}
echo ${indelfiles[0]}

sleep 15

echo "check outputs..."
chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
for ((j=0;j<${#wgstumors[@]};j++))
do
  tumorname=${wgstumors[$j]}
  for ((i=0;i<24;i++))
  do
    chr=chr${chrs[$i]}
    if [[ ! -f ${tumorname}_$chr.snp ]];then echo "$j $i ${tumorname}_$chr.snp is missing!!!";fi
    if [[ ! -f ${tumorname}_$chr.indel ]];then echo "$j $i ${tumorname}_$chr.indel is missing!!!";fi
  done
done
echo "check done"
sleep 20
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
maffile=$mutsigoutputfolder/${outputpreface}.maf
indelfile=$mutsigoutputfolder/${outputpreface}.indel

#generate maf file
echo "oncotator..."
/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/oncotator_varscan.sh ${wgstumors[@]} ${mutationfiles[@]}
/fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/oncotator_varscan.sh ${wgstumors[@]} ${indelfiles[@]}

if [ -f $maffile ]; then rm $maffile; fi
if [ -f $indelfile ]; then rm $indelfile; fi

for ((i=0;i<${#wgstumors[@]};i++))
do
   /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/process_varscan_oncotatorout.sh ${mutationfiles[$i]}.oncotator ${wgstumors[$i]}
   cat ${mutationfiles[$i]}.oncotator.reduced >> $maffile
   /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/process_varscan_oncotatorout.sh ${indelfiles[$i]}.oncotator ${wgstumors[$i]}
   cat ${indelfiles[$i]}.oncotator.reduced >> $indelfile
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
