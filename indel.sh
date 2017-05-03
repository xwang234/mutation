#!/usr/bin/env bash
#SBATCH -t 3-5
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

normalbam=${1?"normalbam"}
tumorbam=${2?"tumorbam"}
tumorname=${3?"tumorname"} #preface of output
outputdir=${4?"outputdir"}

reference=/fh/fast/dai_j/CancerGenomics/Tools/database/reference/full/ucsc.hg19.fasta

gatk_dir=/fh/fast/dai_j/CancerGenomics/Tools/GATK
java_opts="-Xms64m -Xmx4g -Djava.io.tmpdir=`pwd`/tmp"

java $java_opts -jar ${gatk_dir}/IndelGenotyper.36.3336-GenomeAnalysisTK.jar -T IndelGenotyperV2 --input_file:normal $normalbam --input_file:tumor $tumorbam -R $reference --somatic --window_size 300 -et NO_ET -o $outputdir/$tumorname.indel

