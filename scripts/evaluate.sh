#! /bin/bash
#$ -S /bin/bash
#$ -l s_vmem=2G
#$ -N evaluate
#$ -cwd

source ~/.zshrc

TASK_NUM=`expr $SGE_TASK_ID - 1`
echo $TASK_NUM
SUFFIX=`printf %04d $TASK_NUM`
python $SCRIPT_DIR/evaluate.py $WORK_DIR/params.tsv.$SUFFIX $WORK_DIR/out.tsv.$SUFFIX
