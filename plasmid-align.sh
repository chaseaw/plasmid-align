#!/bin/bash

# initiate log file
log_file="align-log.txt"

# activate conda environment
conda activate plasmid-align

# iterates through each folder specified
for folder in "$@"
do

	# looks for required fasta and fastq files
	fafiles=`echo ./$folder/*.fa`
	R1files=`echo ./$folder/*_R1_001.fastq`
	R2files=`echo ./$folder/*_R2_001.fastq`

	# counts how many of each file type is provided
	faVAL=`echo $fafiles | wc -w`
	R1VAL=`echo $R1files | wc -w`
	R2VAL=`echo $R2files | wc -w`

	# halts analysis if more than 1 fasta was provided
	if [ "$faVAL" -gt "1" ]; then
		echo "More than 1 fasta file is present in $folder" >> ${log_file}
		break 1

	# allows for de novo alignment only if no reference fasta is provided
	elif [ "$faVAL" -le "0" ]; then
		echo "de-novo $folder=$( date +%s )" >> ${log_file}
		name1="de-novo"
		
	# confirms fasta selection		
	else
		name1=$fafiles
		
	fi

	# checks that the same number of fastq files were provided, as they have to be matched paired end reads
	if [ "$R1VAL" -ne "$R2VAL" ]; then
		echo "Matched paired-end reads were not provided for $folder" >> ${log_file}
		break 1
	fi

	# stops analysis of folder that has no valid fastqs
	if [ "$R1VAL" -le "0" ]; then
		echo "Illumina-formatted fastqs (_R1_001.fastq) were not found in $folder" >> ${log_file}
		break 1
	fi

	# combines fastqs if multiple runs are provided
	if [ "$R1VAL" -gt "1" ]; then
		cat $R1files > ./$folder/tempR1.fastq
		cat $R2files > ./$folder/tempR2.fastq
		name2=`echo "./$folder/tempR1.fastq"`
		name3=`echo "./$folder/tempR2.fastq"`

	# moves forward if only 1 fastq each was provided
	else
		name2=$R1files
		name3=$R2files
	fi

	# output directory is specified as de_novo
	outname=`echo ./$folder/de_novo`

	# runs de novo alignment using megahit (default is -t 4 cores/threads, -m half (0.5) system memory, and target contig length --min-contig-len > 1000)	
	megahit -t 4 -m 0.5 --min-contig-len 1000 -1 "${name2}" -2 "${name3}" -o "${outname}"
	echo "$folder de novo assembly completed=$( date +%s )" >> ${log_file}

	# if a reference fasta was provided, uses bwa-mem2 to align reads to reference
	if [ ${name1: -3} == ".fa" ]; then

		# indexes the fasta
		bwa-mem2 index "${name1}"
		echo "$folder fasta indexed=$( date +%s )" >> ${log_file}

		# aligns reads and uses samtools to convert to .bam file
		bwa-mem2 mem -t 2 "${name1}" "${name2}" "${name3}" | samtools sort -o "$name1.bam" 
		echo "$folder reads aligned=$( date +%s )" >> ${log_file}

		# creates the .bam index file (.bai) necessary to load into IGV
		samtools index "$name1.bam"
		echo "$folder bam indexed=$( date +%s )" >> ${log_file}
	
		# creates output directory and moves final files there
		workdir=`echo ./$folder`
		cd "$workdir"
		mkdir ref_alignment
		mv *.fa.* ref_alignment
		cd ../
	fi

	# remove temporary files
	if [ "$name2" == `echo "./$folder/tempR1.fastq"` ]; then
		rm "$name2"
		rm "$name3"
	fi

done

# returns to base environment
conda deactivate
