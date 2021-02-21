#!/bin/bash

###############################################################################

####### Change threads to something reasonable for your machine #######
threads=32
#######################################################################

# path to genome dir
genome=/home/prime/Documents/Shin_lab/mm
study_dir=/media/frank/Photography/Big_Data_Overflow/Shin_Lab_Data/RNA-seq
# study_dir=/mnt/i/Big_Data_Overflow/Shin_Lab_Data/RNA-seq
# path to fastq dir
fastq=${study_dir}/fastq
phred=phred33

# path to fastqc
# fastqc=/home/prime/Downloads/FastQC/fastqc
fastqc=/home/frank/FastQC/fastqc

# path to trimmomatic jar
trim=/home/prime/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar
# trim=/home/frank/Trimmomatic-0.39/trimmomatic-0.39.jar
minlen=32
out_format=fastq
# adapters=/home/frank/downloads/trimmomatic/Trimmomatic-0.39/adapters
adapters=/home/prime/Downloads/Trimmomatic-0.39/adapters
# adapters=/home/frank/Trimmomatic-0.39/adapters
adapter_file=TruSeq2-PE.fa
adapter=`echo $adapter_file | cut -d "." -f 1`

# path to results dirs
fastqc_dir_pre=${study_dir}/fastqc_pretrim_results
fastqc_dir_post=${study_dir}/fastqc_posttrim_results_${minlen}
trimmed_dir=${study_dir}/trimmed_${minlen}
star_results=${study_dir}/star_mapped_${minlen}
salmon_results=${study_dir}/salmon_mapped_${minlen}

# path to install of STAR
star=/home/prime/Downloads/STAR-2.7.6a/bin/Linux_x86_64_static/STAR

# path to salmon
salmon=/home/prime/miniconda3/envs/salmon/bin/salmon
salmon_index=/home/prime/Downloads/salmon_index_mm10_cdna/salmon_index
salmon_index2=/home/prime/Downloads/salmon_sa_index_mm10
###############################################################################

echo "Started at `date`"
mkdir -p ${trimmed_dir}
touch trimmed_checksums.txt

while read line; do
	fwd_read=`echo ${line} | cut -d "," -f 1`;
	rev_read=`echo ${line} | cut -d "," -f 2`;
	base_dir=`echo ${line} | cut -d "_" -f 1`_`echo ${line} | cut -d "_" -f 2`;
	echo -e "fwd: ${fwd_read}\nrev: ${rev_read}\nbase_dir: ${base_dir}\n";
	# echo -e "fwd_read: ${fwd_read}\nrev_read ${rev_read}\nbase_dir ${base_dir}\nfastqc_pre dir: ${fastqc_dir_pre}/${base_dir}\nfastqc_post dir: ${fastqc_dir_post}/${base_dir}\n"

	# mkdir -p ${fastqc_dir_pre}/${base_dir}
	# ${fastqc} ${fastq}/${fwd_read} ${fastq}/${rev_read} --outdir=${fastqc_dir_pre}/${base_dir}

	# echo -e "\n-- Trimming files ${fwd_read} and ${rev_read} with ${adapter} --\n";
	# java -jar ${trim} PE -threads $threads -${phred} ${fastq}/${fwd_read} ${fastq}/${rev_read} -baseout ${trimmed_dir}/${base_dir} ILLUMINACLIP:${adapters}/${adapter_file}:2:30:10 LEADING:3 TRAILING:3 MINLEN:${minlen}
	# java -jar ${trim} PE -threads $threads -${phred} -trimlog ${trimmed_dir}/${base_dir}_log.txt ${fastq}/${fwd_read} ${fastq}/${rev_read} -baseout ${trimmed_dir}/${base_dir}_${adapter}.${out_format} ILLUMINACLIP:${adapters}/${adapter_file}:2:30:10 LEADING:3 TRAILING:3 MINLEN:${minlen}

	trimmed_fwd=${trimmed_dir}/${base_dir}_1P_${adapter}_trimmed.${out_format}
	trimmed_rev=${trimmed_dir}/${base_dir}_2P_${adapter}_trimmed.${out_format}
	# echo $trimmed_fwd
	# echo $trimmed_rev

	# md5sum $trimmed_fwd >> trimmed_checksums.txt
	# md5sum $trimmed_rev >> trimmed_checksums.txt

	# mkdir -p ${fastqc_dir_post}/${base_dir}
	# ${fastqc} $trimmed_fwd $trimmed_rev --outdir=${fastqc_dir_post}/${base_dir}

	# ${star} --runThreadN $threads --genomeDir ${genome} --readFilesIn ${trimmed_fwd} ${trimmed_rev} --outFileNamePrefix ${star_results}/${base_dir}/${base_dir}_
	# ${salmon} quant -i ${salmon_index} -l A -1 ${trimmed_fwd} -2 ${trimmed_rev} --validateMappings -o ${salmon_results}/${base_dir}
	# ${salmon} quant -i ${salmon_index2} -l A -1 ${trimmed_fwd} -2 ${trimmed_rev} --validateMappings -o ${salmon_results}_sa_mm10/${base_dir}

	# break

done < fastq_file_pairs.txt

echo "Ended at `date`"
echo "Done!";
