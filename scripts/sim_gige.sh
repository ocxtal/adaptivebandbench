#! /bin/bash

source ~/.zshrc

WORK_DIR=`pwd`
SCRIPT_DIR=`dirname $0`

WORKER_NUM=2000
python $SCRIPT_DIR/generate_params_gige.py > $WORK_DIR/params_gige.tsv
TASK_NUM=`cat $WORK_DIR/params_gige.tsv | wc -l`
split -d -a4 -l `expr $TASK_NUM / $WORKER_NUM + 1` $WORK_DIR/params_gige.tsv $WORK_DIR/params_gige.tsv.
JOB_NUM=`ls -1 $WORK_DIR/params_gige.tsv.* | wc -l`
qsub -t 1-$JOB_NUM:1 -v SCRIPT_DIR=$SCRIPT_DIR -v WORK_DIR=$WORK_DIR $SCRIPT_DIR/evaluate_gige.sh
