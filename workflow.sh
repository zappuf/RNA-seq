#!/bin/bash

##### Read in config file
source ./config.sh
source ./util.sh
###############################################################################

echo "Started at `date`"
mkdir -p ${trimmed_dir}

while read line; do
	fwd_read=`echo ${line} | cut -d "," -f 1`;
	rev_read=`echo ${line} | cut -d "," -f 2`;
	base_dir=`echo ${line} | cut -d "_" -f 1`_`echo ${line} | cut -d "_" -f 2`;
	echo -e "fwd: ${fwd_read}\nrev: ${rev_read}\nbase_dir: ${base_dir}\n";
	# echo -e "fwd_read: ${fwd_read}\nrev_read ${rev_read}\nbase_dir ${base_dir}\nfastqc_pre dir: ${fastqc_dir_pre}/${base_dir}\nfastqc_post dir: ${fastqc_dir_post}/${base_dir}\n"

	mkdir -p ${fastqc_dir_pre}/${base_dir}
	${fastqc} ${fastq}/${fwd_read} ${fastq}/${rev_read} --outdir=${fastqc_dir_pre}/${base_dir}

	echo -e "\n-- Trimming files ${fwd_read} and ${rev_read} with ${adapter} --\n";
	java -jar ${trim} PE -threads $threads ${fastq}/${fwd_read} ${fastq}/${rev_read} -baseout ${trimmed_dir}/${base_dir} ILLUMINACLIP:${adapters}/${adapter_file}:2:30:10 LEADING:3 TRAILING:3 MINLEN:${minlen}
	mv ${trimmed_dir}/${base_dir}_1P ${trimmed_dir}/${base_dir}_1P_${adapter}_trimmed.fastq
	trimmed_fwd=${trimmed_dir}/${base_dir}_1P_${adapter}_trimmed.fastq
	mv ${trimmed_dir}/${base_dir}_2P ${trimmed_dir}/${base_dir}_2P_${adapter}_trimmed.fastq
	trimmed_rev=${trimmed_dir}/${base_dir}_2P_${adapter}_trimmed.fastq

	mkdir -p ${fastqc_dir_post}/${base_dir}
	${fastqc} $trimmed_fwd $trimmed_rev --outdir=${fastqc_dir_post}/${base_dir}

	${salmon} quant -i ${salmon_index} -l A -1 ${trimmed_fwd} -2 ${trimmed_rev} --validateMappings -o ${salmon_results}/${base_dir}

	# break

done < fastq_gz_file_pairs.txt

echo "Ended at `date`"
echo "Done!";
