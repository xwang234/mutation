#!/usr/bin/env bash

#sigfile=output/henan_varscan1_3_signatures.mat
#outfile=henan_varscan1_3_signatures.txt
#sigfile=output/henan_varscan_wgs_3_signatures.mat
#outfile=henan_varscan_3_signatures.txt
#sigfile=output/dulak_varscan_2_signatures.mat
#outfile=dulak_varscan_2_signatures.txt
#sigfile=output/henan_varscan1_2_signatures.mat
#outfile=henan_varscan1_2_signatures.txt
sigfile=output/henan_varscan_wgs_2_signatures.mat
outfile=henan_varscan_2_signatures.txt
sigfile=output/henan_varscan_wgs_1_signatures.mat
outfile=henan_varscan_1_signatures.txt
sigfile=output/henan_varscan1_1_signatures.mat
outfile=henan_varscan1_1_signatures.txt
sigfile=output/dulak_varscan_1_signatures.mat
outfile=dulak_varscan_1_signatures.txt
sigfile=output/dulak_varscan2_2_signatures.mat
outfile=dulak_varscan2_2_signatures.txt
sigfile=output/henan_varscan2_2_signatures.mat
outfile=henan_varscan2_2_signatures.txt
sigfile=output/escc_varscan2_2_signatures.mat
outfile=escc_varscan2_2_signatures.txt
sigfile=output/dulak_varscan2_1_signatures.mat
outfile=dulak_varscan2_1_signatures.txt
sigfile=output/henan_varscan2_1_signatures.mat
outfile=henan_varscan2_1_signatures.txt
sigfile=output/escc_varscan2_1_signatures.mat
outfile=escc_varscan2_1_signatures.txt
sigfile=output/henan_golden4_1_signatures.mat
outfile=henan_golden4_1_signatures.txt
matlab -nodesktop -nodisplay -r "savesignature('$sigfile','$outfile');exit;"
