#! /usr/bin/env bash
#SBATCH -t 0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

#annotation using annovar
#script used to transfer mutect ouput to annovar format
transformat=/fh/fast/dai_j/CancerGenomics/Tools/wang/trans_Varscan_Annovar.pl
annovardir=/fh/fast/dai_j/CancerGenomics/Tools/annovar

#wgsnormals=(SRR1002719 SRR999433 SRR999687 SRR1000378 SRR1002786 SRR999559 SRR1001730 SRR10018461 SRR1002703 SRR999599 SRR1002792 SRR1001839 SRR9994281 SRR1002710 SRR9995631 SRR1021476)
#wgstumors=(SRR1001842 SRR1002713 SRR999423 SRR1001466 SRR1002670 SRR1001823 SRR999489 SRR1002343 SRR1002722 SRR1002656 SRR1002929 SRR999438 SRR1001915 SRR999594 SRR1001868 SRR1001635)
#wgsnormals=(2A 4A 6A 8A 10A 12A)
#wgstumors=(1A 3A 5A 7A 9A 11A)

#wgsnormals=(4A 12A 14A 16A 18A 26A 30A 34A 38A 42A)
#wgstumors=(3A 11A 13A 15A 17A 25A 29A 33A 37A 41A)
#	    0  1   2   3   4   5   6   7   8   9

mutationfile=${1?"snpfile"}
samplename=${2?"samplename"}

echo $mutationfile
echo $samplename

perl $transformat $mutationfile >${mutationfile}_annovar.txt
perl $annovardir/annotate_variation.pl -out ${samplename}_varscan -build hg19 ${mutationfile}_annovar.txt $annovardir/humandb/	
    
for line in ${samplename}_varscan.variant_function
do
  sed -i "s/$/$samplename/" $line
done
for line in ${samplename}_varscan.exonic_variant_function
do
   sed -i "s/$/$samplename/" $line
done

