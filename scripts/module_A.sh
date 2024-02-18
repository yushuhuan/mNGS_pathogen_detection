#!/bin/bash

# 激活虚拟环境
eval "$(conda shell.bash hook)"
conda activate mNGS_pipeline

# 从终端获取参数
sample_id=$1

cd /opt/ossfs/mNGS_pipeline_output/$sample_id

# 根据sample id获取原始数据所在目录，若找不到其原始数据则退出
sample_fastq_dir=$(find /opt/ossfs/mNGS-data/RawData -name $sample_id".fq.gz")
if [ ! -s $sample_fastq_dir ];then
        echo "Error: can't found raw fastq files for $sample_id! Please make sure it has been uploaded to server!"
        exit 1
fi

# 样本质量控制，剔除低质量序列和接头
fastp -i $sample_fastq_dir --thread 16 \
        -o $sample_id.fastp.filtered.fq \
        -h $sample_id.qc_summary.html \
        -j $sample_id.qc_summary.json \
        -n 5 \
        -q 15 \
        -u 40 \
        --length_required 10 \
        --length_limit 60 \
        --complexity_threshold 20 \
        --adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA

# 去除人源序列
bowtie2 -x /opt/ossfs/pathogen_reference_databases/bowtie2_index/homo_sapiens/GRCh38_noalt_as \
        -U $sample_id.fastp.filtered.fq \
        --un $sample_id.removedHuman.fq \
        --sam-no-qname-trunc \
        -S $sample_id.align.result.sam \
        --threads 16

# 比对至参考库
# archaea: RefSeq complete archaeal genomes/proteins
# bacteria: RefSeq complete bacterial genomes/proteins（下载中）
# plasmid: RefSeq plasmid nucleotide/protein sequences（下载失败）
# viral: RefSeq complete viral genomes/proteins（下载中）
# human: GRCh38 human genome/proteins
# fungi: RefSeq complete fungal genomes/proteins
# protozoa: RefSeq complete protozoan genomes/proteins
# nt: NCBI non-redundant nucleotide database(未下载，单独成库)
# UniVec_Core: A subset of UniVec chosen to minimize false positive hits to the vector database

# 比对至standard参考库
kraken2 --db /opt/mNGS/kraken2_index/standard \
        --threads 16 \
        --unclassified-out unclassified_sequences.fa \
        --output $sample_id.classified.standard.output \
        -report $sample_id.classified.standard.report \
        --report-zero-counts \
        --memory-mapping \
        --use-names \
        $sample_id.removedHuman.fq
# 比对至pluspf参考库
kraken2 --db /opt/mNGS/kraken2_index/pluspf \
        --threads 16 \
        --unclassified-out unclassified_sequences.fa \
        --output $sample_id.classified.pluspf.output \
        -report $sample_id.classified.pluspf.report \
        --report-zero-counts \
        --memory-mapping \
        --use-names \
        $sample_id.removedHuman.fq
# 比对至eupathdb参考库
kraken2 --db /opt/mNGS/kraken2_index/eupathdb \
        --threads 16 \
        --unclassified-out unclassified_sequences.fa \
        --output $sample_id.classified.eupathdb.output \
        -report $sample_id.classified.eupathdb.report \
        --report-zero-counts \
        --memory-mapping \
        --use-names \
        $sample_id.removedHuman.fq

# 物种丰度估计(默认在species水平上进行估计)
bracken -d /opt/mNGS/kraken2_index/standard \
        -i $sample_id.classified.standard.report \
        -o $sample_id.abundance.standard.count \
        -w $sample_id.abundance.standard.report

bracken -d /opt/mNGS/kraken2_index/pluspf \
        -i $sample_id.classified.pluspf.report \
        -o $sample_id.abundance.pluspf.count \
        -w $sample_id.abundance.pluspf.report

bracken -d /opt/mNGS/kraken2_index/eupathdb \
        -i $sample_id.classified.eupathdb.report \
        -o $sample_id.abundance.eupathdb.count \
        -w $sample_id.abundance.eupathdb.report

# 结束
echo "Mudule A done!"