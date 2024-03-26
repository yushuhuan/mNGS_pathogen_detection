#!/bin/bash

## 设定参数
sample_id="NCBJ231001"
seqtype="DNA"

# 如果OSS中没有以该样本名命名的文件夹，则创建;有则退出
cd /opt/ossfs/mNGS_pipeline_output/NC_samples
if [ ! -d /opt/ossfs/mNGS_pipeline_output/$sample_id ];then
    mkdir -p "$sample_id"
    cd $sample_id
else
    echo "NC sample:"$sample_id"analysied done!"
    exit 0
fi

# 获取原始数据路径
sample_fastq_dir=$(find /opt/ossfs/mNGS-data/RawData -name $sample_id".fq.gz")
if [ ! $sample_fastq_dir ];then
    echo "Error: can't found raw fastq files for $sample_id! Please make sure it has been uploaded to server!"
    exit 1
else
    echo "Raw fastq exist."
fi

## 原始数据预处理，去除接头及低质量序列
/home/yushuhuan/anaconda3/bin/fastp -i $sample_fastq_dir --thread 8\
            -o $sample_id.fastp.filtered.fastq.gz \
            -h $sample_id.qc_summary.html \
            -j $sample_id.qc_summary.json \
            -n 3 \
            -q 10 \
            -u 30 \
            --length_required 5 \
            --complexity_threshold 10 \
            --adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA

# 去除PCR重复
minirmd -i $sample_id.fastp.filtered.fastq.gz -o $sample_id.rmDup.fastq -d 1

## 分类
/opt/softwares/Kraken2/kraken2 --db /opt/mNGS/kraken2_index/nt \
                --threads 8 \
                --output $sample_id.classified.nt.output \
                -report $sample_id.classified.nt.report \
                --memory-mapping \
                --use-names \
                $sample_id.rmDup.fastq

## 提取Species水平的比对情况
/opt/softwares/Bracken/bracken -d /opt/mNGS/kraken2_index/nt \
        -i $sample_id.classified.nt.report \
        -o $sample_id.species.abundance.nt.count \
        -w $sample_id.species.abundance.nt.report
cat $sample_id.species.abundance.nt.count | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$6,$7}' >$sample_id.species.reads.info

echo "Done!"