#! /bin/bash

source ~/.zshrc

WORK_DIR=`pwd`
SCRIPT_DIR=`dirname $0`

cat $WORK_DIR/out_gap.tsv.* > $WORK_DIR/out_gap.tsv
python $SCRIPT_DIR/aggregate_gap.py $WORK_DIR/out_gap.tsv $WORK_DIR/result_gap.txt
rm params_gap.tsv.* evaluate_gap.o* evaluate_gap.e* out_gap.tsv.* *.maf *.fastq *.ref
