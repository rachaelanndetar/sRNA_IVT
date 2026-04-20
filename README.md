# sRNA analysis

1. Trim adaptor sequences from fastq files in a paired fashion using scripts. 
```
  bash cutadapt_sRNA.sh
```
2. Merge forward and reverse reads using script:
```
  for i in *1.fq.gz; do sbatch ../bbmerge_pairedend.sh $i; done
```
3. Map to references with Bowtie2, try both local and paired end using scripts. Calculate statistics
```
  for i in merge*fq.gz; do temp=${i#merge_trim_}; fa=${temp%%.fq.gz}.fa; sbatch ../bowtie2_local.slm ../reference_sequences/$fa $i; done
  for i in merge*fq.gz; do temp=${i#merge_trim_}; fa=${temp%%.fq.gz}.fa; sbatch ../bowtie2_local.slm ../reference_sequences/$fa $i; done
  for i in *sam; do samtools stat $i > ${i%.sam}_stat.txt; done
```
4. Filter, sort, convert to bam file and index. 
```
  sbatch ../.slm
```
5. Use custom python script to plot 5' and 3' ends, depth at each position in transcript, count fragments for end2end alignments. Wrapped all commands in a SLURM script:
```
  sbatch endparse.slm
```
Running python script looks like this:
```
  endparse.py   -s "fastq/samfiles_end2end/merge_trim_Eco_fmet_filtered.sam"     -d "fastq/samfiles_end2end/merge_trim_Eco_fmet_filtered_sorted_depth.txt"     --ref_len 240     --hh 58     --hdv 135
```
6. Download sorted *bam* files to local, and look at both the end2end and local fragments. Also download plots and tables.
