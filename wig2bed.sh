#!/usr/bin/env bash
#SBATCH -t 1-5
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

wigfile=${1?"wigfile"}
bedfile=${2?"outputfile"}
echo ${SCRATCH_LOCAL}
echo $wigfile

lbedfile=$(basename $bedfile)
export PATH=$PATH:/fh/fast/dai_j/CancerGenomics/Tools/BEDOPS/bin
/fh/fast/dai_j/CancerGenomics/Tools/BEDOPS/bin/wig2bed --do-not-sort < $wigfile > ${SCRATCH_LOCAL}/$lbedfile.1
ls -l ${SCRATCH_LOCAL}/$lbedfile.1
awk '{if ($5 == 1) print $1"\t"$2"\t"$3}' ${SCRATCH_LOCAL}/$lbedfile.1 > ${SCRATCH_LOCAL}/$lbedfile.2
ls -l ${SCRATCH_LOCAL}/$lbedfile.2
/fh/fast/dai_j/CancerGenomics/Tools/bedtools2-master/bin/bedtools merge -i ${SCRATCH_LOCAL}/$lbedfile.2 > $bedfile

 
echo "done"

 
