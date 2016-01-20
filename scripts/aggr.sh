#! /bin/bash

source ~/.zshrc

WORK_DIR=`pwd`
SCRIPT_DIR=`dirname $0`

cat $WORK_DIR/out.tsv.* > $WORK_DIR/out.tsv
python $SCRIPT_DIR/aggregate_id.py $WORK_DIR/out.tsv $WORK_DIR/result.txt
rm params.tsv.* evaluate.o* evaluate.e* out.tsv.*
