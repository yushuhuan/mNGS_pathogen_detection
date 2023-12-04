#!/bin/sh

### 帮助文档
helpdoc(){
    cat <<EOF
Description:
    This script is used to execute one step DNA based mNGS analysis.
Usage:
    bash -i [sample is list] -t [Seq_type:cfDNA,DNA,RNA]
Option:
    -i this is sample id
    -t seq_type, including cfDNA, DNA and RNA
    -s the state of filter, strict or loose
EOF
}
###

while getopts ":i:t:s:v" opt
do
    case "$opt" in
        i) sample=$OPTARG ;;
        t) seq_type=$OPTARG;;
        s) filter_state=$OPTARG;;
        v) SHOW_DETAIL=true ;;
        ?) echo "Unkown parameters!" ; exit 2 ;;
    esac
done
# 若无指定任何参数则输出帮助文档
if [ $# = 0 ]
then
    helpdoc
    exit 1
fi

echo "Processing sample : $sample . start!"

### body ###

input_file=$(ls /opt/ossfs/mNGS-data/$seq_type | grep $sample)
if [ ! $input_file ]
then
    echo "Error: please make sure the raw fastq file exist!"
    exit 1
else
    echo $input_file
fi

fastp -i /opt/ossfs/mNGS-data/$seq_type/$input_file --threads 16 /
        -o /opt/ossfs/DNA-based-output/fastp_filtered_fq/$sample.fastp.filtered.fq /
        -h /opt/ossfs/DNA-based-output/sample_quality_reports/${sample}_qc_summary.html /
        -j /opt/ossfs/DNA-based-output/sample_quality_reports/${sample}_qc_summary.json /
        -n 5 / # default:5, quality control, limit N numbers, if >10 N bases in a read, then dicard it
        -q 15 / # default:15, quality control, a base's phred quality >=Q15 is qualified
        -u 40 / #default:40, quality control, 40% bases allowed to be unqualified
        --length_required 10 / # default:15, length filter, reads shorter than 10bp will be discard
        --length_limit 60 / # default:0, length filter, reads longer than 60bp will be dicard
        --complexity_threshold 20 / # default:30, length filter, 20% complexity is required
        --adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA
        # --adapter_fasta /opt/mNGS/files/ZhiDe_adapter_AMT02-480.fa # adapter trimming, trim the adapters in fastq file one by one
        
blastn -db human_ZhiDe /
        -query /opt/ossfs/DNA-based-output/fastp_filtered_fq/$sample.fastp.filtered.fq /
        -out /opt/ossfs/DNA-based-output/blast/$sample.blast.fq /
        --num_threads 16 /
        -html TRUE /
        -outfmt 6 /
        





bash /opt/mNGS/scripts/taxomonic_profiling.sh -i /opt/ossfs/DNA-based-output/kneaddata_remove_host_condam/${sample}_kneaddata.fastq -s $sample

echo "Analysis done! Please execute downstream analysis!"
