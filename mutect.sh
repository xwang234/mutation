#! /usr/bin/env bash
#SBATCH -t 5-5
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

#source /etc/profile.d/fh_path.sh
#Mutect not work with jdk1.7
module load java/jdk1.6.0_45

cosmic=/fh/fast/dai_j/CancerGenomics/Tools/database/vcf/hg19_cosmic_v54_120711.sorted.vcf
dbsnp=/fh/fast/dai_j/CancerGenomics/Tools/database/vcf/dbsnp_138.hg19.compact.sorted.vcf
gatk_dir=/fh/fast/dai_j/CancerGenomics/Tools/GATK
exon_interval_file=/fh/fast/dai_j/CancerGenomics/Tools/database/exomecapture/S0293689_Padded.intervals
mutect=/fh/fast/dai_j/CancerGenomics/Tools/Mutect/muTect-1.1.4.jar
reference=/fh/fast/dai_j/CancerGenomics/Tools/database/reference/full/ucsc.hg19.fasta
transformat=/fh/fast/dai_j/CancerGenomics/Tools/wang/trans_Mutect_Annovar.pl
annovardir=/fh/fast/dai_j/CancerGenomics/Tools/annovar

java_opts="-Djava.io.tmpdir=`pwd`/tmp -Xmx8g"


normalbam=${1?"normalbam"}
tumorbam=${2?"tumorbam"}
tumorname=${3?"tumorname"} #preface of output
outputdir=${4?"outputdir"}
qmin=${5:-"1"} #min base quality value filter set as 1 10 20
mapmin=${6:-"1"} #min mapping quality score filter

cd $outputdir
echo $normalbam
echo $tumorbam
echo $tumorname
echo ${SLURM_CPUS_ON_NODE}
echo ${SCRATCH_LOCAL}

output=$outputdir/${tumorname}.Mutect_out.txt
outputvcf=$outputdir/${tumorname}.Mutect.vcf
outputwig=$outputdir/${tumorname}.wig.txt

filtermap0=0
if [[ $filtermap0 -eq 1 ]]
then
  samtools view -hb -F 4 -q 1 $normalbam >${SCRATCH_LOCAL}/N_${tumorname}.bam
  ls -l ${SCRATCH_LOCAL}/N_${tumorname}.bam
  samtools index ${SCRATCH_LOCAL}/N_${tumorname}.bam
  samtools view -hb -F 4 -q 1 $tumorbam >${SCRATCH_LOCAL}/${tumorname}.bam
  ls -l ${SCRATCH_LOCAL}/${tumorname}.bam
  samtools index ${SCRATCH_LOCAL}/${tumorname}.bam
  java $java_opts -jar $mutect --analysis_type MuTect --reference_sequence $reference --cosmic $cosmic --dbsnp $dbsnp --min_qscore $qmin --required_maximum_alt_allele_mapping_quality_score $mapmin --input_file:normal ${SCRATCH_LOCAL}/N_${tumorname}.bam --input_file:tumor ${SCRATCH_LOCAL}/${tumorname}.bam --out $output --vcf $outputvcf --coverage_file $outputwig --enable_extended_output #-nct ${SLURM_CPUS_ON_NODE}
else
  java $java_opts -jar $mutect --analysis_type MuTect --reference_sequence $reference --cosmic $cosmic --dbsnp $dbsnp --min_qscore $qmin --required_maximum_alt_allele_mapping_quality_score $mapmin --input_file:normal $normalbam --input_file:tumor $tumorbam --out $output --vcf $outputvcf --coverage_file $outputwig --enable_extended_output #-nct ${SLURM_CPUS_ON_NODE}
fi

grep -v REJECT $output | awk '{if (NR>1) print }' - > ${tumorname}.Mutect_out_keep.txt
perl $transformat ${tumorname}.Mutect_out_keep.txt >${tumorname}.Mutect_annovar.txt
perl $annovardir/annotate_variation.pl -out ${tumorname} -build hg19 ${tumorname}.Mutect_annovar.txt $annovardir/humandb/	

#java $java_opts -jar $mutect --analysis_type MuTect --reference_sequence $reference --cosmic $cosmic --dbsnp $dbsnp --min_qscore 1 --input_file:normal $normalfile --input_file:tumor $tumorfile --out $output --vcf $outputvcf --coverage_file $outputwig -ip 50 --intervals $exon_interval_file
#java $java_opts -jar $mutect --analysis_type MuTect --reference_sequence $reference --min_qscore 20 --input_file:normal $normalfile --input_file:tumor $tumorfile --out $output --vcf $outputvcf --coverage_file $outputwig --fraction_contamination 0.25 -ip 50 --intervals $exon_interval_file

