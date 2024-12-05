#!/bin/bash

set -e  # 确保任何命令失败时脚本退出

# 使用 getopt 解析命令行选项
ARGS=$(getopt -o s:o:n:t:k:h --long NCid:,outputDir:,NASdir:,threads:,kraken2DB:,help -- "$@")

# 显示帮助信息的函数
show_help() {
  echo "使用方法: bash $0 -s|--NCid [NC_id] -o|--outputDir [output_dir] -n|--NASdir [NAS_dir] -t|--threads [threads] -k|--kraken2DB [kraken2DB_dir]"
  echo "选项:"
  echo "  -s, --NCid  对照样本id"
  echo "  -o, --outputDir  结果输出目录"
  echo "  -n, --NASdir  NAS目录，存储下机数据的目录"
  echo "  -t, --threads  线程数"
  echo "  -k, --kraken2DB  kraken2 参考数据库目录"
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
    -s | --NCid)
      NCid=$2
      shift 2
      ;;
    -o | --outputDir)
      outputDir=$2
      shift 2
      ;;
    -n | --NASdir)
      NASdir=$2
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
if [ -z "$NCid" ] || [ -z "$outputDir" ] || [ -z "$NASdir" ] || [ -z "$threads" ] || [ -z "$kraken2DB" ]; then
  echo "错误: 必须提供所有必需的参数。"
  show_help
  exit 1
fi

### 开始处理
cd $outputDir
if [ ! -d $NCid ];then
	mkdir -p "$NCid"
fi
cd $NCid

printf "正在搜索NAS并寻找原始fastq文件……\n"
NC_fq_file=$(find $NASdir -name $NCid".fq.gz")
if [ ! $NC_fq_file ];then
	printf "\033[41:30mError：未找到$NCid.fq.gz,请确认其是否存在，若不存在，请上传至云服务器！\n"
	exit 1
fi

printf "开始质控……\n"
fastp -i $NC_fq_file --thread $threads \
    -o $outputDir/$NCid/$NCid.fastp.filtered.fastq.gz \
    -h $outputDir/$NCid/$NCid.qc_summary.html \
    -j $outputDir/$NCid/$NCid.qc_summary.json \
    -n 5 \
    -q 15 \
    -u 30 \
    --length_required 5 \
    --complexity_threshold 10 \
    --adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA
printf "质控完成\n"

# 去除PCR重复
printf "去除重复……\n"
minirmd -i $outputDir/$NCid/$NCid.fastp.filtered.fastq.gz \
    -o $outputDir/$NCid/$NCid.rmDup.fastq \
    -d 1
printf "去除重复完成\n"

# 分类
printf "开始比对……\n"
kraken2 --db $kraken2DB --threads $threads \
	    --output $outputDir/$NCid/$NCid.classified.nt.output \
	    -report $outputDir/$NCid/$NCid.classified.nt.report \
	    --memory-mapping \
	    --use-names \
	    $NCid.rmDup.fastq

## 提取Species水平的比对情况
Rscript /app/get_genius_and_kingdom_info.r \
    $outputDir/$NCid/$NCid.classified.nt.report \
    $outputDir/$NCid \
    $NCid

cat $outputDir/$NCid/$NCid.pathogen_with_genius_and_domain.txt | awk -v FS="\t" -v OFS="\t" '{print $5,$2,$6}' >$outputDir/$NCid/$NCid.species.reads.info

printf "比对完成!"
