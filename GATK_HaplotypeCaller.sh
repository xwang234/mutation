#! /usr/bin/env bash
#SBATCH -t 1-5
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

gatk_dir=/fh/fast/dai_j/CancerGenomics/Tools/GATK/GATK3.7
java_opts="-Xms64m -Xmx12g -Djava.io.tmpdir=`pwd`/tmp"
dbsnp=/fh/fast/dai_j/CancerGenomics/Tools/database/vcf/resource_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf
reference=/fh/fast/dai_j/CancerGenomics/Tools/database/reference/hg38/hg38.fa
exon_interval_file=/fh/fast/dai_j/CancerGenomics/Tools/database/exomecapture/S0293689_Padded.intervals

normalbam=${1?"normalbam"}
tumorbam=${2?"tumorbam"}
output=${3?"output"}





java -jar $gatk_dir/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference -I $normalbam -I $tumorbam \
--dbsnp $dbsnp \
#-stand_call_conf 30 \
-L $exon_interval_file \
-o $output
