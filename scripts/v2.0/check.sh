#!/bin/bash

rm -rf undone.list checked.output done.list

cat test_patch1_pathogen_validation_TRUE.txt | while read line
do
	sample_id=$(echo $line | awk '{print $1}')
	seqtype=$(echo $line | awk '{print $2}')

	if [ $seqtype == "cfDNA" ];then
        sample_id="$sample_id-CF"
    else
		if [ $seqtype == "RNA" ];then
			sample_id="$sample_id-R"
		fi
    fi

	if [ -s /opt/ossfs/mNGS_pipeline_output/${sample_id}/results/*.abundance.merged.count ];then
		echo $sample_id "done" >>checked.output

		echo $sample_id >>done.list
	else
		raw_fastq_nonexist=$(cat raw.fastq.nonexist.list | grep $sample_id)
		if [ -n "$raw_fastq_nonexist" ];then
			echo $sample_id "raw_fastq_non-exist" >>checked.output
		else
			echo $sample_id "fail" "fastp not done" >>checked.output

			if [ -d /opt/ossfs/mNGS_pipeline_output/${sample_id} ];then
				mv /opt/ossfs/mNGS_pipeline_output/${sample_id} /opt/ossfs/mNGS_pipeline_output/undone
			fi
			cat test_patch1_pathogen_validation_TRUE.txt | grep $sample_id >>undone.list
		fi
	fi
done
