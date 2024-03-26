#!/bin/bash

## 设定参数
sample_id="BZD214240-CF"
seqtype="cfDNA"
NC_sample_id="NCBJ231001"

# if [ $seqtype == "cfDNA" ];then
#     sample_id="$sample_id-CF"
# else
#     if [ $seqtype == "RNA" ];then
#         sample_id="$sample_id-R"
#     fi
# fi

# 如果OSS中没有以该样本名命名的文件夹，则创建
cd /opt/ossfs/mNGS_pipeline_output/test_V2.1
if [ ! -d /opt/ossfs/mNGS_pipeline_output/$sample_id ];then
    mkdir -p "$sample_id"
    cd $sample_id
else
    cd $sample_id
fi

# 获取原始数据路径
sample_fastq_dir=$(find /opt/ossfs/mNGS-data/RawData -name $sample_id".fq.gz")
if [ ! $sample_fastq_dir ];then
    echo "Error: can't found raw fastq files for $sample_id! Please make sure it has been uploaded to server!"
    exit 1
else
    echo "Raw fastq exist."

    # 计算原始序列数
    raw_reads_number=$(pigz -dc $sample_fastq_dir | awk 'NR%4==2{c++} END{print c}')
    echo "原始序列数为：" $raw_reads_number
    # 如果原始序列数 <15M，则报异常
    if [ $raw_reads_number -lt 15000000 ];then
        echo "Reads数异常：" $raw_reads_number
    else
        echo "Reads数正常：" $raw_reads_number
    fi
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
# 计算预处理后的序列数
fastp_filtered_reads_number=$(pigz -dc $sample_id.fastp.filtered.fastq.gz | awk 'NR%4==2{c++} END{print c}')
echo "预处理后序列数为：" $fastp_filtered_reads_number
# 获取Q30比例：
percent_Q30=$(cat $sample_id.qc_summary.html | grep "Q30" | awk 'NR==1{print}' | grep -o "(.*%" | awk -v FS="(" '{print $2}' | awk -v FS="%" '{print $1}')
if [ "0"$percent_Q30 > 85 ];then
    echo "Q30 比例合格，为：" $percent_Q30"%"
else
    echo "Q30 比例不合格为：" $percent_Q30"%"
fi

## 去除人源宿主序列，参考基因组hg38（bwa mem的敏感性比bowtie2高一些，能去除更多的人源序列，提高微生物富集水平）
bwa mem -t 8 /opt/ossfs/pathogen_reference_databases/bwa_index/homo_sapiens/GRCh38_latest_genomic.fna $sample_id.fastp.filtered.fastq.gz >$sample_id.mapping.sam
samtools flagstat $sample_id.mapping.sam >$sample_id.mapping.stat
samtools view -bSh $sample_id.mapping.sam > $sample_id.mapping.bam
samtools sort -@ 8 $sample_id.mapping.bam -o $sample_id.mapping.sorted.bam
samtools index $sample_id.mapping.sorted.bam

# 提取人源序列比例
percent_human=$(cat $sample_id.mapping.stat | awk 'NR==5{print}' | grep -o "(.*%" | awk -v FS="(" '{print $2}' | awk -v FS="%" '{print $1}')
if [ "0"$percent_human > 90 ];then
    echo "人源比例合格，为：" $percent_human"%"
else
    echo "人源比例不合格为：" $percent_human"%"
fi

# 获取性别信息
chrY_info=$(cat $sample_id.mapping.sam | grep "^@SQ" | grep "NC_000024.10")
if [ ! -n "$chrY_info" ];then
    echo "样本性别为：女"
    sex="female"
else
    echo "样本性别为：男"
    sex="male"
fi

# 如果是DNA样本，则保留非宿主序列；如果是RNA样本，则保留宿主序列
if [ $seqtype == "RNA" ];then
    # 保留人源序列，后面做定量
    samtools view -b -F 4 $sample_id.mapping.sorted.bam >$sample_id.human.bam
    samtools sort -n -@ 8 $sample_id.human.bam -O BAM >$sample_id.human.sorted.bam
    # samtools fastq -@ 8 -n $sample_id.human.sorted.bam $sample_id.human.fastq.gz
else
    # 过滤掉比对上的reads，并在排序后转为fastq
    samtools view -b -f 4 $sample_id.mapping.sorted.bam >$sample_id.unmapped.bam
    samtools flagstat $sample_id.unmapped.bam >$sample_id.unmapped.stat
    # 计算去除人源序列后的序列数
    microbiome_reads_number=$(cat $sample_id.unmapped.stat | awk -v FS=" " 'NR==1{print $1}')
    echo "去除人源序列后的序列数为：" $microbiome_reads_number

    # picard去除PCR重复(其实如果是RNA的话一般不去重复，因为RNA数据起始量高，且存在某些基因表达高某些基因表达低的情况；但考虑到mNGS实验过程中起始量低，PCR扩增数较高，因此去除重复能更好的校正之后的微生物相对丰度)
    java -jar /opt/softwares/Picard/picard.jar MarkDuplicates -I $sample_id.unmapped.bam \
        -O $sample_id.unmapped.rmDup.bam \
        -M $sample_id.unmapped.rmDup.metrics.txt \
        --REMOVE_DUPLICATES TRUE
    samtools flagstat $sample_id.unmapped.rmDup.bam >$sample_id.unmapped.rmDup.stat 
    # 计算去除PCR重复后的序列数
    PCR_filtered_reads_number=$(cat $sample_id.unmapped.rmDup.stat | awk -v FS=" " 'NR==1{print $1}')
    echo "去重复后序列数为：" $PCR_filtered_reads_number

    samtools sort -@ 8 $sample_id.unmapped.rmDup.bam -o $sample_id.unmapped.rmDup.sorted.bam
    samtools fastq -@ 8 -n $sample_id.unmapped.rmDup.sorted.bam >$sample_id.unmapped.rmDup.sorted.fastq
fi
#################################################################################################################

mkdir microbiome_results
cd microbiome_results

## 处理微生物序列

## 比对到微生物数据库
# cd /opt/ossfs/mNGS_pipeline_output/$sample_id/microbiome_results
# mkdir centrifuge
# cd centrifuge
# option 1:centrifuge, 参考库：NCBI nt库


# option 2: krakenUniq，参考库：archaea, bacteria, viral, human, UniVec_Core, Eukaryotic pathogen genomes (EuPathDB54) with contaminants removed
# cd /opt/ossfs/mNGS_pipeline_output/test_V2.1/$sample_id/microbiome_results
# mkdir krakenUniq
# cd krakenUniq
# /opt/softwares/KrakenUniq/krakenuniq --report-file $sample_id.classified.nt.report \
#                 --db /opt/mNGS/krakenUniq_index/microbialdb \
#                 --threads 8 \
#                 --output $sample_id.classified.nt.output \
#                 /opt/ossfs/mNGS_pipeline_output/test_V2.1/$sample_id/$sample_id.unmapped.rmDup.sorted.fastq
# classified_reads_krakenUniq=$(cat $sample_id.classified.nt.report | awk 'NR==6{print $2}')
# echo "krakenUniq 分类序列数：" $classified_reads_krakenUniq

# option 3: kraken2, 参考库：inclusive of GenBank, RefSeq, TPA and PDB
cd /opt/ossfs/mNGS_pipeline_output/test_V2.1/$sample_id/microbiome_results
mkdir kraken2
cd kraken2 
/opt/softwares/Kraken2/kraken2 --db /opt/mNGS/kraken2_index/nt \
                --threads 8 \
                --output $sample_id.classified.nt.output \
                -report $sample_id.classified.nt.report \
                --memory-mapping \
                --use-names \
                /opt/ossfs/mNGS_pipeline_output/test_V2.1/$sample_id/$sample_id.unmapped.rmDup.sorted.fastq
# kraken2分类序列数
classified_reads_kraken2=$(cat $sample_id.classified.nt.report | head -1 | awk '{print $2}')
echo "kraken2 分类序列数：" $classified_reads_kraken2

cd /opt/ossfs/mNGS_pipeline_output/$sample_id/microbiome_results/kraken2
## 对分类结果进行过滤，标注和中英文标注
# 1  提取全部species分类病原的genius,domain,reads_num,rank等信息，并分别计算基于domain和species的相对丰度
Rscript /opt/mNGS/run2/get_genius_and_kingdom_info.r \
        "/opt/ossfs/mNGS_pipeline_output/test_V2.1/$sample_id/microbiome_results/kraken2/$sample_id.classified.nt.report" \
        "/opt/ossfs/mNGS_pipeline_output/test_V2.1/$sample_id/microbiome_results/kraken2"              
# 2.对每个水平下的所有病原进行标注，然后过滤
# 2.1 联合NC的分类结果
# 2.2 检查眼科高关注病原和常见背景微生物，与NC比对，报出关注度：高，低或背景
Rscript /opt/mNGS/run2/estimate_attention_index.r \
        "/opt/ossfs/mNGS_pipeline_output/test_V2.1/$sample_id/microbiome_results/kraken2/pathogen_with_genius_and_domain.txt" \
        "/opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_sample_id/NCBJ231001.species.reads.info" \
        "/opt/ossfs/mNGS_pipeline_output/test_V2.1/$sample_id/microbiome_results/kraken2"


####################################################################################################

if [ $seqtype == "RNA" ];then

    mkdir human_expression_profile
    cd human_expression_profile

    # 宿主基因表达定量
    RSEM 

    # 免疫因子表达信息提取（与正常样本相比）

    # 宿主免疫信息特征标注
    # 如果真菌相关表达特征提高，则报“疑似真菌感染”

fi