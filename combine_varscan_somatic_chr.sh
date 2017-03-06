#!/usr/bin/env bash
#run it in the varscanoutput folder

tumors=( "$@" )
echo ${tumors[@]}

chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

varscandir=/fh/fast/dai_j/CancerGenomics/Tools/VarScan
varscancmd=$varscandir/VarScan.v2.3.9.jar

#check results 
echo "check outputs..."
for ((i=0;i<${#tumors[@]};i++))
do
  tumorname=${tumors[$i]}
  for ((chrid=0;chrid<24;chrid++))
  do
    chr=chr${chrs[$chrid]}
    #echo $chr
    if [[ ! -f ${tumorname}_$chr.snp ]];then echo "${tumorname}_$chr.snp is missing!!!";fi
  done
done
sleep 10

#first process 
#echo "filter mutations..."
#for ((i=0;i<${#tumors[@]};i++))
#do
#  tumorname=${tumors[$i]}
#  for chrid in {0..23}
#  do
#    chr=chr${chrs[$chrid]}
#    java -Xmx256m -jar $varscancmd processSomatic ${tumorname}_$chr.snp --min-tumor-freq 0.1 --max-normal-freq 0.2 --p-value 0.05
#  done
#done

echo "combine results..."
for ((i=0;i<${#tumors[@]};i++))
do
  tumorname=${tumors[$i]}
  for chrid in {0..23}
  do
    chr=chr${chrs[$chrid]}
    if [[ $chrid -eq 0 ]]
    then 
      cat ${tumorname}_$chr.snp.Somatic.hc > ${tumorname}.snp.Somatic.hc
    else
      awk '{if (NR>1) print}' ${tumorname}_$chr.snp.Somatic.hc >> ${tumorname}.snp.Somatic.hc  #remove title
    fi
  done
done

