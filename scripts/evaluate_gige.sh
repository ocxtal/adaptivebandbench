#! /bin/bash
#$ -S /bin/bash
#$ -l s_vmem=2G
#$ -N evaluate_gige
#$ -cwd

source ~/.zshrc

TASK_NUM=`expr $SGE_TASK_ID - 1`
echo $TASK_NUM
SUFFIX=`printf %04d $TASK_NUM`
python $SCRIPT_DIR/evaluate_gige.py $WORK_DIR/params_gige.tsv.$SUFFIX $WORK_DIR/out_gige.tsv.$SUFFIX
