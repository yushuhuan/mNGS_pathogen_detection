#!/bin/bash

set -e  # 确保任何命令失败时脚本退出

# 使用 getopt 解析命令行选项
ARGS=$(getopt -o s:l:n:d:o:t:k:a:b:h --long sampleID:,library:,NCid:,NCdir:,outputDir:,threads:,kraken2DB:,NASdir:,bwaIndex:,help -- "$@")

# 显示帮助信息的函数
show_help() {
  echo "使用方法: bash $0 -s|--sampleID [sampleID] -l|--library [single|paired] -n|--NCid [NCid] -d|--NCdir [NCdir] -o|--outputDir [outputDir] -t|--threads [threads] -k|--kraken2DB [kraken2DB_dir]"
  echo "选项:"
  echo "  -s, --sampleID  样本id"
  echo "  -l, --library  建库方式，single或者paired"
  echo "  -n, --NCid  对照样本id，如果没有NC，则设为none即可"
  echo "  -d, --NCdir  对照样本输出文件存储目录"
  echo "  -o, --outputDir  输出文件存储目录"
  echo "  -t, --threads  线程数"
  echo "  -k, --kraken2DB  kraken2所用数据库对应目录"
  echo "  -a, --NASdir  NAS目录"
  echo "  -b, --bwaIndex  bwa index hg38参考基因组对应目录+prefix,eg."
  echo "  -h, --help    显示帮助信息"
  exit 0
}

# 如果没有传递任何参数，则显示帮助信息
if [ $# -eq 0 ]; then
  show_help
  exit 0
fi

# 如果参数解析成功，传递给 eval 处理
eval set -- "$ARGS"

# 循环处理每个选项
while true; do
  case "$1" in
    -s | --sampleID)
      sampleID=$2
      shift 2
      ;;
    -l | --library)
      library=$2
      shift 2
      ;; 
    -n | --NCid)
      NCid=$2
      shift 2
      ;;
    -d | --NCdir)
      NCdir=$2
      shift 2
      ;;
    -o | --outputDir)
      outputDir=$2
      shift 2
      ;;
    -t | --threads)
      threads=$2
      shift 2
      ;;
    -k | --kraken2DB)
      kraken2DB=$2
      shift 2
      ;;
    -a | --NASdir)
      NASdir=$2
      shift 2
      ;;
    -b | --bwaIndex)
      bwaIndex=$2
      shift 2
      ;;
    --)
      shift
      break
      ;;
    -h | --help)
      show_help
      ;;
    *)
      echo "无效的选项 $1"
      show_help
      ;;
  esac
done

# 参数验证：确保所有必需的参数都传入
if [ -z "$sampleID" ] || [ -z "$library" ] || [ -z "$NCid" ] || [ -z "$NCdir" ] || [ -z "$outputDir" ] || [ -z "$threads" ] || [ -z "$kraken2DB" ] || [ -z "$bwaIndex" ] || [ -z "$NASdir" ]; then
  echo "错误: 必须提供所有必需的参数。"
  show_help
  exit 1
fi

# 打印已解析的参数
echo "解析的参数："
echo "样本ID: $sampleID"
echo "建库方式：$library"
echo "对照样本ID: $NCid"
echo "NAS目录：$NASdir"
echo "对照样本目录: $NCdir"
echo "输出目录: $outputDir"
echo "线程数: $threads"
echo "Kraken2 数据库: $kraken2DB"
echo "BWA 索引: $bwaIndex"

echo -e "\n"

start_time=$(date "+%Y-%m-%d %H:%M:%S")
cd $outputDir
if [ ! -d $sampleID ];then
    mkdir -p "$sampleID"
fi
cd $sampleID

printf "正在处理样本id为:$sampleID, 阴性对照样本(NC)id为:$NCid的任务……\n\n"

### 首先处理NC样本
printf "$sampleID: Step1-对照样本分析\n" >$outputDir/$sampleID/$sampleID.step.log

printf "正在分析阴性对照样本$NCid……\n"
if [ $NCid == "none" ];then
    printf "NCid不存在，将直接开始处理样本\n\n"
else
  if [ -s $NCdir/$NCid/$NCid.species.reads.info ];then
      printf "$NCid分析已完成!\n\n"
  else
      bash NC_analysis.sh -s $NCid -o $NCdir -n $NASdir -t $threads -k $kraken2DB
  fi
fi

### 开始处理样本
printf "开始处理样本$sampleID……\n"

# 寻找服务器上原始fastq文件 ###
printf "$sampleID: Step2-样本原始fastq文件搜寻\n" >>$outputDir/$sampleID/$sampleID.step.log

### 单端与双端的代码略有不同
if [ $library == "single" ];then
    printf "正在搜索NAS并寻找原始fastq文件……\n"
    sample_fq_file=$(find $NASdir -name $sampleID".fq.gz")
    if [ ! $sample_fq_file ];then
            printf "Error：未找到$sampleID.fq.gz,请确认其是否存在，若不存在，请上传至NAS或确认样本id！\n"
            exit 1
    else
            priintf "样本fastq文件目录为：$sample_fq_file\n\n"
    fi

    cd $outputDir/$sampleID

    # 计算原始序列数
    raw_reads_number=$(pigz -dc $sample_fq_file | awk 'NR % 4 == 2 {c++} END {print c}')
    printf "原始序列reads数为$raw_reads_number\n" >>$outputDir/$sampleID/$sampleID.step.log
    # 如果原始序列数 <15M，则报异常
    if [ $raw_reads_number -lt 15000000 ];then
        printf "请注意，原始序列reads数小于15M不合格！\n\n" >>$outputDir/$sampleID/$sampleID.step.log    
    fi

    # 原始数据预处理，去除接头及低质量序列
    printf "$sampleID: Step3-样本质控\n" >>$outputDir/$sampleID/$sampleID.step.log
    printf "开始质控……\n"
    timeout 1800 fastp -i $sample_fq_file --thread $threads \
                -o $outputDir/$sampleID/$sampleID.fastp.filtered.fastq.gz \
                -h $outputDir/$sampleID/$sampleID.qc_summary.html \
                -j $outputDir/$sampleID/$sampleID.qc_summary.json \
                -n 5 \
                -q 15 \
                -u 30 \
                --length_required 5 \
                --complexity_threshold 10 \
                --adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
                2>/dev/null
    # 保存 timeout 的退出状态
    timeout_status=$?
    # 检查 timeout 的退出状态
    if [ $timeout_status -eq 124 ]; then
        echo "FASTP command exceeded 0.5 hours and was terminated."
        exit 1
    fi
    # 保存 fastp 的退出状态
    fastp_status=$?
    # 检查 fastp 的退出状态
    if [ $fastp_status -ne 0 ]; then
        echo "FASTP command failed with exit status $fastp_status!"
        exit 1
    fi
    # 如果 fastp 执行成功
    echo "FASTP completed successfully."
    printf "质控完成\n\n"

    # 计算预处理后的序列数
    fastp_filtered_reads_number=$(pigz -dc $outputDir/$sampleID/$sampleID.fastp.filtered.fastq.gz | awk 'NR%4==2{c++} END{print c}')
    printf "预处理后序列数为:$fastp_filtered_reads_number\n" >>$outputDir/$sampleID/$sampleID.step.log
    # 获取Q30比例：
    percent_Q30=$(cat $outputDir/$sampleID/$sampleID.qc_summary.html | grep "Q30" | awk 'NR==1{print}' | grep -o "(.*%" | awk -v FS="(" '{print $2}' | awk -v FS="%" '{print $1}')

    printf "Q30比例为:$percent_Q30\n" >>$outputDir/$sampleID/$sampleID.step.log

    ## 去除人源宿主序列，参考基因组hg38（bwa mem的敏感性比bowtie2高一些，能去除更多的人源序列，提高微生物富集水平）
    printf "$sampleID: Step4-去除人源序列\n" >>$outputDir/$sampleID/$sampleID.step.log
    printf "开始去除人源序列……\n"
    timeout 3600 bwa mem -t $threads $bwaIndex $outputDir/$sampleID/$sampleID.fastp.filtered.fastq.gz >$outputDir/$sampleID/$sampleID.mapping.sam
    # 保存 timeout 的退出状态
    timeout_status=$?
    # 检查 timeout 的退出状态
    if [ $timeout_status -eq 124 ]; then
        echo "BWA command exceeded 0.5 hours and was terminated."
        exit 1
    fi
    # 保存 bwa 的退出状态
    bwa_status=$?
    # 检查 bwa 的退出状态
    if [ $bwa_status -ne 0 ]; then
        echo "BWA command failed with exit status $fastp_status!"
        exit 1
    fi
    # 如果 bwa 执行成功
    echo "BWA completed successfully."
else
    printf "正在搜索NAS并寻找原始fastq文件……\n"
    sample_fq_file_1=$(find $NASdir -name $sampleID"_raw_1.fq.gz")
    if [ ! $sample_fq_file_1 ];then
            printf "Error：未找到${sampleID}_raw_1.fq.gz,请确认其是否存在，若不存在，请上传至NAS或确认样本id！\n"
            exit 1
    else
            echo "样本fastq_1文件目录为：$sample_fq_file_1"
    fi
    sample_fq_file_2=$(find $NASdir -name $sampleID"_raw_2.fq.gz")
    if [ ! $sample_fq_file_2 ];then
            printf "Error：未找到${sampleID}_raw_2.fq.gz,请确认其是否存在，若不存在，请上传至NAS或确认样本id！\n"
            exit 1
    else
            printf "样本fastq_2文件目录为：$sample_fq_file_2\n\n"
    fi

    cd $outputDir/$sampleID

    # 原始数据预处理，去除接头及低质量序列
    printf "$sampleID: Step3-样本质控\n" >>$outputDir/$sampleID/$sampleID.step.log
    printf "开始质控……\n"
    timeout 1800 fastp -i $sample_fq_file_1 -I $sample_fq_file_2 --thread $threads \
                -o $outputDir/$sampleID/${sampleID}_1.fastp.filtered.fastq.gz \
                -O $outputDir/$sampleID/${sampleID}_2.fastp.filtered.fastq.gz \
                -h $outputDir/$sampleID/$sampleID.qc_summary.html \
                -j $outputDir/$sampleID/$sampleID.qc_summary.json \
                -n 5 \
                -q 15 \
                -u 30 \
                --length_required 5 \
                --complexity_threshold 10 \
                --adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
                2>/dev/null
    # 保存 timeout 的退出状态
    timeout_status=$?
    # 检查 timeout 的退出状态
    if [ $timeout_status -eq 124 ]; then
        echo "FASTP command exceeded 0.5 hours and was terminated."
        exit 1
    fi
    # 保存 fastp 的退出状态
    fastp_status=$?
    # 检查 fastp 的退出状态
    if [ $fastp_status -ne 0 ]; then
        echo "FASTP command failed with exit status $fastp_status!"
        exit 1
    fi
    # 如果 fastp 执行成功
    echo "FASTP completed successfully."
    printf "质控完成\n\n"

    ## 去除人源宿主序列，参考基因组hg38（bwa mem的敏感性比bowtie2高一些，能去除更多的人源序列，提高微生物富集水平）
    printf "$sampleID: Step4-去除人源序列\n" >>$outputDir/$sampleID/$sampleID.step.log
    printf "开始去除人源序列……\n"
    timeout 3600 bwa mem -t $threads $bwaIndex $outputDir/$sampleID/${sampleID}_1.fastp.filtered.fastq.gz $outputDir/$sampleID/${sampleID}_2.fastp.filtered.fastq.gz >$outputDir/$sampleID/$sampleID.mapping.sam
    # 保存 timeout 的退出状态
    timeout_status=$?
    # 检查 timeout 的退出状态
    if [ $timeout_status -eq 124 ]; then
        echo "BWA command exceeded 0.5 hours and was terminated."
        exit 1
    fi
    # 保存 bwa 的退出状态
    bwa_status=$?
    # 检查 bwa 的退出状态
    if [ $bwa_status -ne 0 ]; then
        echo "BWA command failed with exit status $fastp_status!"
        exit 1
    fi
    # 如果 bwa 执行成功
    echo "BWA completed successfully."
fi

samtools flagstat -@ 8 $outputDir/$sampleID/$sampleID.mapping.sam >$outputDir/$sampleID/$sampleID.mapping.stat

# 提取人源序列比例
percent_human=$(cat $outputDir/$sampleID/$sampleID.mapping.stat | awk 'NR==5{print}' | grep -o "(.*%" | awk -v FS="(" '{print $2}' | awk -v FS="%" '{print $1}')
printf "人源序列比例为$percent_human\n" >>$outputDir/$sampleID/$sampleID.step.log
# 使用 bc 进行浮点数比较
if (( $(echo "$percent_human < 90" | bc -l) )); then
    printf "请注意，人源比例小于90%不合格！\n\n" >>$outputDir/$sampleID/$sampleID.step.log
    echo "请注意，人源比例小于90%不合格！"
fi

# 获取性别信息
chrY_info=$(cat $outputDir/$sampleID/$sampleID.mapping.sam | grep "^@SQ" | grep "NC_000024.10")
if [ ! -n "$chrY_info" ];then
    printf "样本性别为：女\n" >>$outputDir/$sampleID/$sampleID.step.log
    echo "样本性别为：女"
else
    printf "样本性别为：男\n" >>$outputDir/$sampleID/$sampleID.step.log
    echo "样本性别为：男"
fi

# 过滤掉比对上的reads，并在排序后转为fastq
timeout 1800 samtools view -@ $threads -b -f 4 $outputDir/$sampleID/$sampleID.mapping.sam >$outputDir/$sampleID/$sampleID.unmapped.bam
# 保存 timeout 的退出状态
timeout_status=$?
# 检查 timeout 的退出状态
if [ $timeout_status -eq 124 ]; then
    echo "SAMTOOLS command exceeded 0.5 hours and was terminated."
    exit 1
fi
# 保存 bwa 的退出状态
samtools_status=$?
# 检查 samtools 的退出状态
if [ $samtools_status -ne 0 ]; then
    echo "SAMTOOLS command failed with exit status $fastp_status!"
    exit 1
fi
# 如果 samtools 执行成功
echo "SAMTOOLS completed successfully."

samtools sort -@ $threads $outputDir/$sampleID/$sampleID.unmapped.bam >$outputDir/$sampleID/$sampleID.unmapped.sorted.bam
samtools flagstat -@ $threads $outputDir/$sampleID/$sampleID.unmapped.bam >$outputDir/$sampleID/$sampleID.unmapped.stat
samtools fastq -@ 8 -n $outputDir/$sampleID/$sampleID.unmapped.sorted.bam >$outputDir/$sampleID/$sampleID.unmapped.sorted.fastq

# 计算去除人源序列后的序列数
microbiome_reads_number=$(cat $outputDir/$sampleID/$sampleID.unmapped.stat | awk -v FS=" " 'NR==1{print $1}')
printf "去除人源序列后的序列数为：$microbiome_reads_number\n" >>$outputDir/$sampleID/$sampleID.step.log
printf "人源序列去除完成\n\n"

#################################################################################################################

printf "$sampleID: Step6-病原分类\n" >>$outputDir/$sampleID/$sampleID.step.log
printf "开始分类……\n"
timeout 1800 kraken2 --db $kraken2DB \
                --threads $threads \
                --output $outputDir/$sampleID/$sampleID.classified.nt.output \
                -report $outputDir/$sampleID/$sampleID.classified.nt.report \
                --memory-mapping \
                --use-names \
                $outputDir/$sampleID/$sampleID.unmapped.sorted.fastq
# 保存 timeout 的退出状态
timeout_status=$?
# 检查 timeout 的退出状态
if [ $timeout_status -eq 124 ]; then
    echo "KRAKEN2 command exceeded 0.5 hours and was terminated."
    exit 1
fi
# 保存 bwa 的退出状态
kraken2_status=$?
# 检查 kraken2 的退出状态
if [ $kraken2_status -ne 0 ]; then
    echo "KRAKEN2 command failed with exit status $fastp_status!"
    exit 1
fi
# 如果 kraken2 执行成功
echo "KRAKEN2 completed successfully."
printf "分类完成\n"
# kraken2分类序列数
classified_reads_kraken2=$(cat $outputDir/$sampleID/$sampleID.classified.nt.report | awk 'NR==2{print $2}')
printf "分类序列数：$classified_reads_kraken2\n" 

## 对分类结果进行过滤，标注和中英文标注
printf "$sampleID: Step7-个性化注释\n\n" >>$outputDir/$sampleID/$sampleID.step.log
# 1  提取全部species分类病原的genius,domain,reads_num,rank等信息，并分别计算基于domain和species的相对丰度
Rscript /opt/mNGS/ZhiDe_mngs_docker/get_genius_and_kingdom_info.r \
        $outputDir/$sampleID/$sampleID.classified.nt.report \
        $outputDir/$sampleID \
        $sampleID             

if [ $NCid == "none" ];then
    Rscript /opt/mNGS/ZhiDe_mngs_docker/estimate_attention_index_without_NC.r \
            $outputDir/$sampleID/$sampleID.pathogen_with_genius_and_domain.txt \
            $outputDir/$sampleID \
            $sampleID
else
    # 2.对每个水平下的所有病原进行标注，然后过滤
    # 2.1 联合NC的分类结果
    # 2.2 检查眼科高关注病原和常见背景微生物，与NC比对，报出关注度：高，低或背景
    Rscript /opt/mNGS/ZhiDe_mngs_docker/estimate_attention_index_with_NC.r \
            $outputDir/$sampleID/$sampleID.pathogen_with_genius_and_domain.txt \
            $outputDir/$NCid/$NCid.species.reads.info \
            $outputDir/$sampleID \
            $sampleID
fi

  end_time=$(date "+%Y-%m-%d %H:%M:%S")

  printf "$sampleID 分析开始时间为$start_time\n" >>$outputDir/$sampleID/$sampleID.step.log
  printf "$sampleID 分析结束时间为$end_time\n\n" >>$outputDir/$sampleID/$sampleID.step.log

  printf "任务完成！" >>$outputDir/$sampleID/$sampleID.step.log