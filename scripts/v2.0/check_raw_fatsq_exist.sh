#!/bin/bash

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

    sample_fastq_dir=$(find /opt/ossfs/mNGS-data/RawData -name "$sample_id.fq.gz")
    if [ ! $sample_fastq_dir ];then
        echo $sample_id "raw_fastq_non-exist" >>checked.raw.fastq.output
    else
        echo $sample_id "ok" >>checked.raw.fastq.output
    fi
done