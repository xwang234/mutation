#!/usr/bin/env bash
#SBATCH -t 1-0
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

snpfile=${1?"filename"}
echo $snpfile

oncotator=/fh/fast/dai_j/CancerGenomics/Tools/oncotator1/oncotator-1.8.0.0/oncotator/bin/oncotator
oncotator_db=/fh/fast/dai_j/CancerGenomics/Tools/oncotator/oncotator_v1_ds_June112014

#make sure the column names are suitable for oncotator
title=()
title=($(awk '{if (NR==1) print }' $snpfile))
title[0]=chr
title[1]=position
title[2]=ref_allele
title[3]=alt_allele
newtitle=""
for ((i=0;i<${#title[@]};i++))
do
  if [[ $i -lt ${#title[@]}-1 ]]
  then
    newtitle="$newtitle${title[$i]}\t"
  else
	newtitle="$newtitle${title[$i]}\n"
  fi
done
echo -e $newtitle > ${snpfile}_tmp

awk '{if (NR>1) print }' $snpfile >> ${snpfile}_tmp
#to deal with indel files column 4 has "+" and "-"
awk 'BEGIN{FS="\t"; OFS="\t"} {sub(/\+/,"",$4);print;}' ${snpfile}_tmp > ${snpfile}_tmp1
awk 'BEGIN{OFS="\t"} {sub(/^\-[A,T,G,C]+/,"-",$4);print;}' ${snpfile}_tmp1 > ${snpfile}_tmp2
#sleep 10
$oncotator -v --db-dir $oncotator_db ${snpfile}_tmp2 $snpfile.oncotator hg19
rm ${snpfile}_tmp ${snpfile}_tmp1 ${snpfile}_tmp2


