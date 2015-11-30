#! /bin/bash


WORKER_NUM=100
python generate_params.py > params.tsv
TASK_NUM=`cat params.tsv | wc -l`
split -d -a4 -l `expr $TASK_NUM / $WORKER_NUM` params.tsv params.tsv.
python evaluate.py params.tsv out.tsv.0

cat out.tsv.* > out.tsv
