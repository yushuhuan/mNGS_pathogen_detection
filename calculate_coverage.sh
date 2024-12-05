#!/bin/bash

set -e  # 确保任何命令失败时脚本退出

# 使用 getopt 解析命令行选项
ARGS=$(getopt -o s:d:p:t:r:h --long sampleID:,sampleDir:,pathogen:,threads:,reference:,help -- "$@")

# 显示帮助信息的函数
show_help() {
  echo "使用方法: bash $0 -s|--sampleID [sample_id] -d|--sampleDir [sample_mngs_output_dir] -p|--pathogen [pathogen_name] -t|--threads [threads] -r|--reference [ZhiDe_ref_pathogen_DB]"
  echo "选项:"
  echo "  -s, --sampleID  样本id"
  echo "  -d, --sampleDir  样本mngs结果目录"
  echo "  -p, --pathogen  病原名-科学名，例如Cutibacterium acnes"
  echo "  -t, --threads  线程数"
  echo "  -r,  --reference 智德病原参考库目录"
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
    -d | --sampleDir)
      sampleDir=$2
      shift 2
      ;;
    -p | --pathogen)
      pathogen=$2
      shift 2
      ;;
    -t | --threads)
      threads=$2
      shift 2
      ;;
    -r | --reference)
      reference=$2
      shift 2
      ;;
    -h | --help)
      show_help
      ;;
    --)
      shift
      break
      ;;
    *)
      echo "无效的选项 $1"
      show_help
      ;;
  esac
done

# 参数验证：确保所有必需的参数都传入
if [ -z "$sampleID" ] || [ -z "$sampleDir" ] || [ -z "$pathogen" ] || [ -z "$reference" ] || [ -z "$threads" ]; then
  echo "错误: 必须提供所有必需的参数。"
  show_help
  exit 1
fi

echo "解析的参数："
echo "样本id为：$sampleID"
echo "样本mngs结果目录为：$sampleDir"
echo "病原名为：$pathogen"
echo "线程数为：$threads"
echo "智德病原参考库目录为：$reference"
echo -e "\n"

### 计算所有比对到智德病原数据库的病原的覆盖度
cd $sampleDir
mkdir coverage_output
cd coverage_output

pathogen_ref=${pathogen// /_}

mkdir $sampleDir/coverage_output/$pathogen_ref

echo "$sampleID-Step1:提取病原对应的reads序列"
if [ -s $sampleDir/$sampleID.classified.nt.output ];then
  cat $sampleDir/$sampleID.classified.nt.output | grep "$pathogen" | awk -v FS="\t" '{print $2}' | seqtk subseq $sampleDir/$sampleID.unmapped.sorted.fastq - > $sampleDir/coverage_output/$pathogen_ref/$sampleID.$pathogen_ref.fastq
else
  echo "样本$sampleID的病原分类结果不存在，请确认该病原的$sampleDir/$sampleID.classified.nt.output文件存在！"
  exit 1
fi

echo "$sampleID-Step2:比对获取病原序列"
bwa mem -t $threads $reference/$pathogen_ref/$pathogen_ref $sampleDir/coverage_output/$pathogen_ref/$sampleID.$pathogen_ref.fastq >$sampleDir/coverage_output/$pathogen_ref/$sampleID.$pathogen_ref.sam
samtools view -bS $sampleDir/coverage_output/$pathogen_ref/$sampleID.$pathogen_ref.sam > $sampleDir/coverage_output/$pathogen_ref/$sampleID.$pathogen_ref.bam
samtools sort $sampleDir/coverage_output/$pathogen_ref/$sampleID.$pathogen_ref.bam -o $sampleDir/coverage_output/$pathogen_ref/$sampleID.$pathogen_ref.sorted.bam

echo "$sampleID-Step3:计算覆盖度"
/opt/softwares/qualimap_v2.3/qualimap bamqc -bam $sampleDir/coverage_output/$pathogen_ref/$sampleID.$pathogen_ref.sorted.bam -outdir $sampleDir/coverage_output/$pathogen_ref

echo "$sampleID-Step4:删除中间文件，保留病原的原始序列"
rm $sampleDir/coverage_output/$pathogen_ref/$sampleID.$pathogen_ref.sam
rm $sampleDir/coverage_output/$pathogen_ref/$sampleID.$pathogen_ref.bam
rm $sampleDir/coverage_output/$pathogen_ref/$sampleID.$pathogen_ref.sorted.bam

printf "覆盖度计算完成!\n"