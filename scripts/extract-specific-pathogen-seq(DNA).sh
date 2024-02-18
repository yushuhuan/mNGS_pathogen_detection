#!/bin/bash

cd /opt/ossfs/mNGS_pipeline_output/$sample_id

# 用户上传病原原始参考序列、病原 taxid、病原科学名称
patho_taxid="10335"
patho_raw_name="Human alphaherpesvirus 3"
group="viral" # Choose from: ['all', 'archaea', 'bacteria', 'fungi', 'invertebrate', 'metagenomes', 'plant', 'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral']


# 将病原原始名中的空格替换为下划线
patho_name=$(echo ${patho_raw_name// /_})

mkdir $patho_name
cd $patho_name

# 比对
# 是否有该病原bowtie2的索引？没有则创建
if [ ! -s /opt/ossfs/pathogen_reference_databases/bowtie2_index/$patho_name ];then
    cd /opt/ossfs/pathogen_reference_databases/bowtie2_index
    mkdir $patho_name
    cd $patho_name
    ncbi-genome-download --formats fasta --taxids $patho_taxid -o ./ --flat-output $group
    gunzip *.gz

    reference_in=$(ls ./ | grep "fna")

    # 建立bowtie2索引
    bowtie2-build -f --threads 16  -q $reference_in $patho_name
fi

cd /opt/ossfs/mNGS_pipeline_output/$sample_id/$patho_name

bowtie2 -x /opt/ossfs/pathogen_reference_databases/bowtie2_index/$patho_name/$patho_name \
        -U ../$sample_id.removedHuman.fq \
        --al  $sample_id.align.$patho_name.fa \
        --sam-no-qname-trunc \
        -S $sample_id.$patho_name.align.result.sam \
        --threads 16 \
        --quiet
# 统计比对结果
samtools flagstat $sample_id.$patho_name.align.result.sam >$sample_id.align.stat.txt
# 过滤sam文件，去除未比对上的比对结果(第二列为4表示未比对上)
cat $sample_id.$patho_name.align.result.sam | awk '{if($2 != 4) print}' >$sample_id.$patho_name.sam
# 去除mapping quality <10的比对结果(第五列为比对质量要求大于10)
cat $sample_id.$patho_name.sam | awk '{if($5 > 10) print}' >$sample_id.$patho_name.mq10.sam
# sam文件转为bam
samtools view -@ 16 -b -T /opt/ossfs/pathogen_reference_databases/bowtie2_index/$patho_name/*.fna \
            $sample_id.$patho_name.mq10.sam >$sample_id.$patho_name.mq10.bam
# 按照原Fasta在文件中的顺序排序
samtools sort -@ 16 $sample_id.$patho_name.mq10.bam >$sample_id.$patho_name.mq10.sorted.bam
# 对bam文件建立索引方便快速处理
samtools index -@ 16 $sample_id.$patho_name.mq10.sorted.bam
# 去除重复序列，包括PCR扩增和一些重复序列
java -jar /opt/softwares/Picard/picard.jar MarkDuplicates -I $sample_id.$patho_name.mq10.sorted.bam \
    -O $sample_id.$patho_name.mq10.rmDup.bam \
    -M $sample_id.$patho_name.mq10.rmDup.metrics.txt \
    --REMOVE_DUPLICATES TRUE
# 再次排序
samtools sort -@ 16 $sample_id.$patho_name.mq10.rmDup.bam >$sample_id.$patho_name.mq10.rmDup.sorted.bam

# 转为fasta格式输出
samtools fasta -@ 16 $sample_id.$patho_name.mq10.rmDup.sorted.bam >$sample_id.$patho_name.mq10.rmDup.sorted.fa

# 该病原基因组覆盖度报出
# 转为bed文件
# bedtools bamtobed -i $sample_id.$patho_name.mq10.rmDup.sorted.fa >$sample_id.$patho_name.mq10.rmDup.sorted.bed
mkdir depth_coverage_stat
java -jar -Xmx6g /opt/softwares/BAMStats-1.25/BAMStats-1.25.jar \
        -d -l -m -q -s -v html \
        -i $sample_id.$patho_name.mq10.rmDup.sorted.bam \
        -o ./depth_coverage_stat

# Coverage.html
# EditDistances.html
# MappedCoverage.html
# MappingQuality.html
# ReadLengths.html
# StartPositions.html

# 分析完成
echo "extract specific pathogen sequences done!"