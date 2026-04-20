#!/bin/bash
#SBATCH --job-name=bbmerge
#SBATCH --time=72:00:00
#SBATCH --output=%j_bbmerge.out



bbmerge.sh in1=$1 in2=${1%1.fq.gz}2.fq.gz out=merge_${1%_1.fq.gz}.fq.gz outu=${1%_1.fq.gz}.unmerged1.fq outu2=${1%_1.fq.gz}.unmerged2.fq.gz ihist=${1%_1.fq.gz}.bbmerge_hist.txt minoverlap=20 mismatches=5
