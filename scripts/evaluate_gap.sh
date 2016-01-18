#! /bin/bash
#$ -S /bin/bash
#$ -l s_vmem=2G
#$ -N evaluate_gap
#$ -cwd

source ~/.zshrc

TASK_NUM=`expr $SGE_TASK_ID - 1`
echo $TASK_NUM
SUFFIX=`printf %04d $TASK_NUM`
python $SCRIPT_DIR/evaluate_gap.py $WORK_DIR/params_gap.tsv.$SUFFIX $WORK_DIR/out_gap.tsv.$SUFFIX
