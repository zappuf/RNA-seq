#!/bin/bash

###############################################################################
# path to genome dir
# study_dir=/media/frank/Photography/Big_Data_Overflow/Shin_Lab_Data/RNA-seq
study_dir=~/RNA-seq
# path to fastq dir
fastq=/fastq
# path to fastqc
fastqc=/opt/rna-seq/FastQC/fastqc

# path to user installed
user_installed=/opt/rna-seq/

# path to trimmomatic jar
trim=${user_installed}/Trimmomatic-0.39/trimmomatic-0.39.jar

# The minimum length you will accept reads
minlen=32

# adapters=/home/frank/downloads/trimmomatic/Trimmomatic-0.39/adapters
adapters=${user_installed}/Trimmomatic-0.39/adapters
adapter_file=TruSeq2-PE.fa
adapter=`echo $adapter_file | cut -d "." -f 1`

# path to results dir
fastqc_dir_pre=${study_dir}/fastqc_pretrim_results
fastqc_dir_post=${study_dir}/fastqc_posttrim_results_${minlen}
trimmed_dir=${study_dir}/trimmed_${minlen}
salmon_results=${study_dir}/salmon_mapped_${minlen}
####### Change threads to something reasonable for your machine!!! #######
threads=32
# path to salmon
salmon=salmon
###############################################################################
