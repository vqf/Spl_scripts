#!/bin/bash
c=$PWD"/*.bam"
s=$(dirname $0)"/"
for l in $c; do ( perl $s"RNAseq_raw.pl" $l & ); done 
