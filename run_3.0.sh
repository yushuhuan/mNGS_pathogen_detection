#!/bin/bash

### 传入定参数
sample_id=$1
NC_id=$2

printf "\033[40:32m正在处理样本id为:$sample_id, 阴性对照样本(NC)id为$NC_id的任务……\033[0m\n\n"

### 首先处理NC样本
printf "\033[40:32m正在分析阴性对照样本$NC_id……\033[0m\n"
if [ -s /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.species.reads.info ];then
    printf "\033[40:32m$NC_id分析已完成\033[0m\n\n"
else
    printf "\033[40:32m开始分析……\033[0m\n"

    cd /opt/ossfs/mNGS_pipeline_output/NC_samples
    if [ ! -d $NC_id ];then
        mkdir -p "$NC_id"
    fi
    cd $NC_id

    printf "\033[40:32m正在搜索云服务器并寻找原始fastq文件……\033[0m\n"
    NC_fq_file=$(find /opt/ossfs/mNGS-data/RawData -name $NC_id".fq.gz")
    if [ ! $NC_fq_file ];then
        printf "\033[41:30mError：未找到$NC_id.fq.gz,请确认其是否存在，若不存在，请上传至云服务器！\033[0m\n"
        exit 1
    fi

    printf "\033[40:32m开始质控……\033[0m\n"
    fastp -i $NC_fq_file --thread 8 \
            -o /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.fastp.filtered.fastq.gz \
            -h /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.qc_summary.html \
            -j /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.qc_summary.json \
            -n 5 \
            -q 15 \
            -u 30 \
            --length_required 5 \
            --complexity_threshold 10 \
            --adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA
    printf "\033[40:32m质控完成\033[0m\n"

    # 去除PCR重复
    printf "\033[40:32m去除PCR重复……\033[0m\n"
    minirmd -i /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.fastp.filtered.fastq.gz \
            -o /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.rmDup.fastq \
            -d 1
    printf "\033[40:32m去除PCR重复完成\033[0m\n"

    # 分类
    printf "\033[40:32m开始比对……\033[0m\n"
    kraken2 --db /opt/mNGS/kraken2_index/nt --threads 8 \
                    --output /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.classified.nt.output \
                    -report /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.classified.nt.report \
                    --memory-mapping \
                    --use-names \
                    $NC_id.rmDup.fastq

    ## 提取Species水平的比对情况
    Rscript /opt/mNGS/ZhiDe-mNGS-analysis-V3/get_genius_and_kingdom_info.r \
            /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.classified.nt.report \
            /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id \
            $NC_id

    cat /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.pathogen_with_genius_and_domain.txt | awk -v FS="\t" -v OFS="\t" '{print $5,$2,$6}' >/opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.species.reads.info

    printf "\033[40:32m比对完成\033[0m\n\n"
fi

### 开始处理样本
printf "\033[40:32m开始处理样本$sample_id……\033[0m\n\n"

# 寻找服务器上原始fastq文件 ###
printf "\033[40:32m正在搜索云服务器并寻找原始fastq文件……\033[0m\n\n"
sample_fq_file=$(find /opt/ossfs/mNGS-data/RawData -name $sample_id".fq.gz")
if [ ! $sample_fq_file ];then
    printf "\033[41:30mError：未找到$sample_id.fq.gz,请确认其是否存在，若不存在，请上传至云服务器！\033[0m\n"
    exit 1
fi

cd /opt/ossfs/mNGS_pipeline_output/V3
if [ ! -d $sample_id ];then
    mkdir -p "$sample_id"
fi
cd $sample_id

# 计算原始序列数
raw_reads_number=$(pigz -dc $sample_fq_file | awk 'NR % 4 == 2 {c++} END {print c}')
printf "\033[43:30m原始序列reads数为$raw_reads_number\033[0m\n\n"
# 如果原始序列数 <15M，则报异常
if [ $raw_reads_number -lt 15000000 ];then
    printf "\033[43:30m请注意，原始序列reads数小于15M不合格！\033[0m\n\n"    
fi


# 原始数据预处理，去除接头及低质量序列
printf "\033[40:32m开始质控……\033[0m\n"
fastp -i $sample_fq_file --thread 8 \
            -o /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.fastp.filtered.fastq.gz \
            -h /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.qc_summary.html \
            -j /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.qc_summary.json \
            -n 5 \
            -q 15 \
            -u 30 \
            --length_required 5 \
            --complexity_threshold 10 \
            --adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA
printf "\033[40:32m质控完成\033[0m\n\n"

# 计算预处理后的序列数
fastp_filtered_reads_number=$(pigz -dc /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.fastp.filtered.fastq.gz | awk 'NR%4==2{c++} END{print c}')
printf "\033[43:30m预处理后序列数为:$fastp_filtered_reads_number\033[0m\n"
# 获取Q30比例：
percent_Q30=$(cat /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.qc_summary.html | grep "Q30" | awk 'NR==1{print}' | grep -o "(.*%" | awk -v FS="(" '{print $2}' | awk -v FS="%" '{print $1}')

printf "\033[43:30mQ30比例为:$percent_Q30\033[0m\n\n"

## 去除人源宿主序列，参考基因组hg38（bwa mem的敏感性比bowtie2高一些，能去除更多的人源序列，提高微生物富集水平）
printf "\033[40:32m开始去除人源序列……\033[0m\n"
bwa mem -t 8 /opt/ossfs/pathogen_reference_databases/bwa_index/homo_sapiens/GRCh38_latest_genomic.fna /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.fastp.filtered.fastq.gz >/opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.mapping.sam

samtools flagstat -@ 8 /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.mapping.sam >/opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.mapping.stat

# 提取人源序列比例
percent_human=$(cat /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.mapping.stat | awk 'NR==5{print}' | grep -o "(.*%" | awk -v FS="(" '{print $2}' | awk -v FS="%" '{print $1}')

printf "\033[43:30m人源序列比例为$percent_human\033[0m\n\n"

# 使用 bc 进行浮点数比较
if (( $(echo "$percent_human < 90" | bc -l) )); then
    printf "\033[43:30m请注意，人源比例小于90%不合格！\033[0m\n\n"
fi

# 获取性别信息
chrY_info=$(cat /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.mapping.sam | grep "^@SQ" | grep "NC_000024.10")
if [ ! -n "$chrY_info" ];then
    printf "\033[43:30m样本性别为：女\033[0m\n"
else
    printf "\033[43:30m样本性别为：男\033[0m\n"
fi

# 过滤掉比对上的reads，并在排序后转为fastq
samtools view -@ 8 -b -f 4 /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.mapping.sam >/opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.bam
samtools sort -@ 8 /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.bam >/opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.sorted.bam
samtools flagstat -@ 8 /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.bam >/opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.stat
# 计算去除人源序列后的序列数
microbiome_reads_number=$(cat $sample_id.unmapped.stat | awk -v FS=" " 'NR==1{print $1}')
printf "\033[43:30m去除人源序列后的序列数为：$microbiome_reads_number\033[0m\n\n"
printf "\033[40:32m人源序列去除完成\033[0m\n\n"

# picard去除PCR重复(其实如果是RNA的话一般不去重复，因为RNA数据起始量高，且存在某些基因表达高某些基因表达低的情况；但考虑到mNGS实验过程中起始量低，PCR扩增数较高，因此去除重复能更好的校正之后的微生物相对丰度)
printf "\033[40:32m开始去重复……\033[0m\n"
java -jar /opt/softwares/Picard/picard.jar MarkDuplicates \
    -I /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.sorted.bam \
    -O /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.rmDup.bam \
    -M /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.rmDup.metrics.txt \
    --REMOVE_DUPLICATES TRUE
printf "\033[40:32m去重复完成\033[0m\n\n"

samtools flagstat -@ 8 /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.rmDup.bam >/opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.rmDup.stat 
# 计算去除PCR重复后的序列数
PCR_filtered_reads_number=$(cat $sample_id.unmapped.rmDup.stat | awk -v FS=" " 'NR==1{print $1}')
printf "\033[43:30m去重复后序列数为：$PCR_filtered_reads_number\033[0m\n\n"

samtools sort -@ 8 /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.rmDup.bam -o /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.rmDup.sorted.bam
samtools fastq -@ 8 -n /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.rmDup.sorted.bam >/opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.rmDup.sorted.fastq

#################################################################################################################

printf "\033[40:32m开始分类……\033[0m\n"
kraken2 --db /opt/mNGS/kraken2_index/nt \
                --threads 8 \
                --output /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.classified.nt.output \
                -report /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.classified.nt.report \
                --memory-mapping \
                --use-names \
                /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.unmapped.rmDup.sorted.fastq
printf "\033[40:32m分类完成\033[0m\n"
# kraken2分类序列数
classified_reads_kraken2=$(cat $sample_id.classified.nt.report | awk 'NR==2{print $2}')
printf "\033[43:30mkraken2 分类序列数：$classified_reads_kraken2\033[0m\n\n"

## 对分类结果进行过滤，标注和中英文标注
# 1  提取全部species分类病原的genius,domain,reads_num,rank等信息，并分别计算基于domain和species的相对丰度
Rscript /opt/mNGS/ZhiDe-mNGS-analysis-V3/get_genius_and_kingdom_info.r \
        /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.classified.nt.report \
        /opt/ossfs/mNGS_pipeline_output/V3/$sample_id \
        $sample_id             
# 2.对每个水平下的所有病原进行标注，然后过滤
# 2.1 联合NC的分类结果
# 2.2 检查眼科高关注病原和常见背景微生物，与NC比对，报出关注度：高，低或背景
Rscript /opt/mNGS/ZhiDe-mNGS-analysis-V3/estimate_attention_index.r \
        /opt/ossfs/mNGS_pipeline_output/V3/$sample_id/$sample_id.pathogen_with_genius_and_domain.txt \
        /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id/$NC_id.species.reads.info \
        /opt/ossfs/mNGS_pipeline_output/V3/$sample_id \
        $sample_id

printf "\033[40:32m任务完成！\033[0m\n"