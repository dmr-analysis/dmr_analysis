#!/bin/bash

#make gene annontation files by using hmst-seq-analyzer
#https://hmst-seq.github.io/hmst/
OUT_FOLDER='../../final_demo_data/fl_12samples/out_data/DMR_CpG_context/'
/usr/bin/time -p hmst_seq_analyzer gene_annotation -F ${OUT_FOLDER} -i no -l 10 \
	-xL 50000000 -X 5000 -Y 1000 -M 5000 -N 1000000 -hu yes -n no \
	-r ../../final_demo_data/genomes/hg19/hg19.refFlat.txt \
	-g ../../final_demo_data/genomes/hg19/hg19.chrom.sizes.clear.sorted
echo gene_annotation-DON

