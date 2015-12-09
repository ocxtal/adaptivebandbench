#! /bin/bash

source ~/.zshrc

WORK_DIR=`pwd`
SCRIPT_DIR=`dirname $0`

WORKER_NUM=2000
python $SCRIPT_DIR/generate_params.py > $WORK_DIR/params.tsv
TASK_NUM=`cat $WORK_DIR/params.tsv | wc -l`
split -d -a4 -l `expr $TASK_NUM / $WORKER_NUM + 1` $WORK_DIR/params.tsv $WORK_DIR/params.tsv.
JOB_NUM=`ls -1 $WORK_DIR/params.tsv.* | wc -l`
# python $SCRIPT_DIR/evaluate.py $WORK_DIR/params.tsv $WORK_DIR/out.tsv.0
qsub -t 1-$JOB_NUM:1 -v SCRIPT_DIR=$SCRIPT_DIR -v WORK_DIR=$WORK_DIR $SCRIPT_DIR/evaluate.sh
# SGE_TASK_ID=3 SCRIPT_DIR=$SCRIPT_DIR WORK_DIR=$WORK_DIR $SCRIPT_DIR/evaluate.sh
