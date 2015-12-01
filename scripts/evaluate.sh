#! /bin/bash
#$ -S /bin/bash
#$ -l s_vmem=2G -l mem_req=2G
#$ -t 0-99:1
#$ -N evaluate
#$ -cwd

SUFFIX=`print %04d $SGE_TASK_ID`
python evaluate.py params.tsv.$SUFFIX out.tsv.$SUFFIX