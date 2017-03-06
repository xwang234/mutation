#! /usr/bin/env bash
#SBATCH -t 3-5
#SBATCH --mail-type=FAIL
#SBATCH --mem=32G
#SBATCH --mail-user=xwang234@fhcrc.org
#use mutect2 in GATK
normalfile=${1?"normalbamfile"}
tumorfile=${2?"tumorbamfile"}
output=${3?"output"}
#change from mutect_new.sh
source /etc/profile.d/fh_path.sh
module load java/jdk1.8.0_31 #GATK3.6
cosmic=/fh/fast/dai_j/CancerGenomics/Tools/database/vcf/hg19_cosmic_v54_120711.sorted.vcf
dbsnp=/fh/fast/dai_j/CancerGenomics/Tools/database/vcf/dbsnp_138.hg19.compact.sorted.vcf
gatk_dir=/fh/fast/dai_j/CancerGenomics/Tools/GATK/GATK3.6
exon_interval_file=/fh/fast/dai_j/CancerGenomics/Tools/database/exomecapture/S0293689_Padded.intervals
#mutect=/fh/fast/dai_j/CancerGenomics/Tools/Mutect/muTect-1.1.4.jar
reference=/fh/fast/dai_j/CancerGenomics/Tools/database/reference/full/ucsc.hg19.fasta
java_opts="-Djava.io.tmpdir=`pwd`/tmp -Xmx8g"
echo $normalfile
echo $tumorfile
echo $output
java $java_opts -jar $gatk_dir/GenomeAnalysisTK.jar -T MuTect2 --reference_sequence $reference --cosmic $cosmic --dbsnp $dbsnp --input_file:normal $normalfile --input_file:tumor $tumorfile --out $output
