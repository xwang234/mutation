#! /usr/bin/env bash

varscandir=/fh/fast/dai_j/CancerGenomics/Tools/VarScan

#wgstumors=(3A 11A 13A 15A 17A 23A 25A 29A 33A 35A 37A 41A 113A 115A 133A 135A)
#	    0  1   2   3   4   5   6   7   8   9   10  11  12   13   14   15
snpfile=${1?"varscansnpfile"}
java -Xmx256m -jar $varscandir/VarScan.v2.3.9.jar processSomatic $snpfile --min-tumor-freq 0.1 --max-normal-freq 0.1 --p-value 0.05

