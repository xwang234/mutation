#!/usr/bin/env bash

normals=($(awk '{print $2}' /fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/TCGA_WESBAM/bamfilelist.txt))
tumors=($(awk '{print $1}' /fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/TCGA_WESBAM/bamfilelist.txt))

bamdir=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/data/TCGA_WESBAM
outputdir=/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/tmp
for ((i=0;i<${#normals[@]};i++))
  do
      normalbam=$bamdir/${normals[$i]}
      tumorbam=$bamdir/${tumors[$i]}
      output=$outputdir/${tumors[$i]}.vcf.txt
      echo /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/GATK_HaplotypeCaller.sh $normalbam $tumorbam $output
      sbatch /fh/fast/dai_j/CancerGenomics/Tools/wang/mutation/GATK_HaplotypeCaller.sh $normalbam $tumorbam $output  
      sleep 1
done

