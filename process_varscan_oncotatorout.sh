#! /usr/bin/env bash

#only keep the necessary info of oncotator output to reduce the file size. not used.
maffile=${1?"the name of maffile"}
wgstumor=${2?"samplename"}
#oncotator output
echo $wgstumor
echo "process oncotator on $maffile"

if [[ ! -f $maffile && -f $maffile.gz ]];then gunzip $maffile.gz ; fi

output=$maffile.reduced

# $261:normal_reads1,$264:normal_reads2,$267:normal_var_freq,$270:somatic_p_value,$273:tumor_reads1,$276:tumor_reads2,$279:tumor_var_freq
#awk 'BEGIN { FS="\t"; OFS="\t" }; { if (NR==4 || (NR>3 && $1 !="Unknown" && ($261+$264)>=10 && ($273+$276)>=10 && $270 < 0.2) ) {print $1,$5,$6,$7,$9,$10,$11,$12,$13,$16,$66,$253,$261,$262,$263,$264,$265,$266,$267,$269,$270,$271,$273,$274,$275,$276,$277,$278,$279,$280}}' $maffile >$output
awk 'BEGIN { FS="\t"; OFS="\t" }; { if (NR>3) {print $1,$5,$6,$7,$9,$10,$11,$12,$13,$16,$66,$253,$261,$262,$263,$264,$265,$266,$267,$269,$270,$271,$273,$274,$275,$276,$277,$278,$279,$280}}' $maffile >$output
for line in $output
do
  sed -i "s/__UNKNOWN__/$wgstumor/" $line
done

