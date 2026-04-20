#!/bin/bash


for i in *1.fq.gz; do cutadapt -a AGATCGGAAGAGCACACGTCT -A GATCGTCGGACTGTAGAACTCTG --minimum-length 20:20 -o trim_$i -p trim_${i%_1.fq.gz}_2.fq.gz $i ${i%_1.fq.gz}_2.fq.gz; done

