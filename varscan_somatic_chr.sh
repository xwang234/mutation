#!/usr/bin/env bash

#SBATCH -t 0-4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

normalbam=${1?"normalbam"}
tumorbam=${2?"tumorbam"}
chr=${3?"chrid,0-23"}
normalname=$(basename $normalbam)
normalname=${normalname%%.*}
tumorname=$(basename $tumorbam)
tumorname=${tumorname%%.*}

echo $normalbam
echo $tumorbam

chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
chr=chr${chrs[$chr]}

varscandir=/fh/fast/dai_j/CancerGenomics/Tools/VarScan
varscancmd=$varscandir/VarScan.v2.3.9.jar
reference=/fh/fast/dai_j/CancerGenomics/Tools/database/reference/full/ucsc.hg19.fasta
#source /etc/profile.d/fh_path.sh
#module load java/jdk1.7.0_25

echo "samtools mpileup -B -q 1 -f $reference -r $chr $normalbam >${SCRATCH_LOCAL}/${normalname}_$chr.mpileup"
samtools mpileup -B -q 1 -f $reference -r $chr $normalbam >${SCRATCH_LOCAL}/${normalname}_$chr.mpileup
samtools mpileup -B -q 1 -f $reference -r $chr $tumorbam >${SCRATCH_LOCAL}/${tumorname}_$chr.mpileup
echo "java -Xmx512m -jar $varscancmd somatic ${SCRATCH_LOCAL}/${normalname}_$chr.mpileup ${SCRATCH_LOCAL}/${tumorname}_$chr.mpileup ${tumorname}_$chr -–min-coverage 10 --min-var-freq 0.1 -–somatic-p-value 0.05 --strand-filter 0  --min-avg-qual 0"
java -Xmx512m -jar $varscancmd somatic ${SCRATCH_LOCAL}/${normalname}_$chr.mpileup ${SCRATCH_LOCAL}/${tumorname}_$chr.mpileup ${tumorname}_$chr --min-coverage 10 --min-var-freq 0.1 --somatic-p-value 0.05 --strand-filter 0 --min-avg-qual 0
echo "java -Xmx256m -jar $varscancmd processSomatic ${tumorname}_$chr.snp --min-tumor-freq 0.1 --p-value 0.05"
java -Xmx256m -jar $varscancmd processSomatic ${tumorname}_$chr.snp --min-tumor-freq 0.1 --p-value 0.05

echo "java -Xmx256m -jar $varscancmd processSomatic ${tumorname}_$chr.indel --min-tumor-freq 0.1 --p-value 0.05"
java -Xmx256m -jar $varscancmd processSomatic ${tumorname}_$chr.indel --min-tumor-freq 0.1 --p-value 0.05





