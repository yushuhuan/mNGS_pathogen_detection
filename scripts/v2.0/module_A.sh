#!/bin/bash

# 激活虚拟环境
#eval "$(conda shell.bash hook)"
#conda activate mNGS_pipeline

# 从终端获取参数
sample_id=$1

#cd /opt/ossfs/mNGS_pipeline_output/$sample_id

# 根据sample id获取原始数据所在目录，若找不到其原始数据则退出
sample_fastq_dir=$(find /opt/ossfs/mNGS-data/RawData -name $sample_id".fq.gz")
if [ ! -s $sample_fastq_dir ];then
        echo "Error: can't found raw fastq files for $sample_id! Please make sure it has been uploaded to server!"
        exit 1
else
        # 样本质量控制，剔除低质量序列和接头
        /home/yushuhuan/anaconda3/bin/fastp -i $sample_fastq_dir --thread 16 \
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
        /usr/bin/bowtie2 -x /opt/ossfs/pathogen_reference_databases/bowtie2_index/homo_sapiens/GRCh38_noalt_as \
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
        /opt/softwares/Kraken2/kraken2 --db /opt/mNGS/kraken2_index/standard \
                --threads 16 \
                --unclassified-out unclassified_sequences.fa \
                --output $sample_id.classified.standard.output \
                -report $sample_id.classified.standard.report \
                --report-zero-counts \
                --memory-mapping \
                --use-names \
                $sample_id.removedHuman.fq
        # 比对至pluspf参考库
        /opt/softwares/Kraken2/kraken2 --db /opt/mNGS/kraken2_index/pluspf \
                --threads 16 \
                --unclassified-out unclassified_sequences.fa \
                --output $sample_id.classified.pluspf.output \
                -report $sample_id.classified.pluspf.report \
                --report-zero-counts \
                --memory-mapping \
                --use-names \
                $sample_id.removedHuman.fq
        # 比对至eupathdb参考库
        /opt/softwares/Kraken2/kraken2 --db /opt/mNGS/kraken2_index/eupathdb \
                --threads 16 \
                --unclassified-out unclassified_sequences.fa \
                --output $sample_id.classified.eupathdb.output \
                -report $sample_id.classified.eupathdb.report \
                --report-zero-counts \
                --memory-mapping \
                --use-names \
                $sample_id.removedHuman.fq

        # 物种丰度估计(默认在species水平上进行估计)
        /opt/softwares/Bracken/bracken -d /opt/mNGS/kraken2_index/standard \
                -i $sample_id.classified.standard.report \
                -o $sample_id.abundance.standard.count \
                -w $sample_id.abundance.standard.report

        /opt/softwares/Bracken/bracken -d /opt/mNGS/kraken2_index/pluspf \
                -i $sample_id.classified.pluspf.report \
                -o $sample_id.abundance.pluspf.count \
                -w $sample_id.abundance.pluspf.report

        /opt/softwares/Bracken/bracken -d /opt/mNGS/kraken2_index/eupathdb \
                -i $sample_id.classified.eupathdb.report \
                -o $sample_id.abundance.eupathdb.count \
                -w $sample_id.abundance.eupathdb.report

        # 原始分类结果合并
        if [ -s $sample_id.abundance.eupathdb.count ];then
                /home/yushuhuan/anaconda3/bin/python /opt/softwares/Bracken/analysis_scripts/combine_bracken_outputs.py \
                        --files $sample_id.abundance.standard.count $sample_id.abundance.pluspf.count $sample_id.abundance.eupathdb.count \
                        --names standard,plusdf,eupathdb \
                        -o $sample_id.raw.abundance.merged.count
        else
                /home/yushuhuan/anaconda3/bin/python /opt/softwares/Bracken/analysis_scripts/combine_bracken_outputs.py \
                        --files $sample_id.abundance.standard.count $sample_id.abundance.pluspf.count \
                        --names standard,plusdf,eupathdb \
                        -o $sample_id.raw.abundance.merged.count
        fi

        # 分类结果过滤
        dbs=(standard pluspf eupathdb)
        for db in ${dbs[@]}
        do
                # 过滤掉人的分类结果
                /home/yushuhuan/anaconda3/bin/python /opt/softwares/KrakenTools/filter_bracken.out.py -i $sample_id.abundance.$db.count \
                        -o $sample_id.abundance.$db.exclude.human.count \
                        --exclude 9606
                # 提取眼科相关微生物丰度结果(第一版本中包含254种眼科常见微生物，修改--include参数，可改变最终报出的微生物丰度列表)
                /home/yushuhuan/anaconda3/bin/python /opt/softwares/KrakenTools/filter_bracken.out.py -i $sample_id.abundance.$db.exclude.human.count \
                        -o $sample_id.abundance.eyePathos.$db.count \
                        --include 129951 130310 130308 162425 5518 5059 5507 746128 38323 573 1352 10376 10298 10359 32603 \
                        10335 5476 562 139 160 5811 470 480 485 727 40324 287 29459 550 615 1282 1313 1351 28037 1280 1290 \
                        1283 29388 37329 1307 1292 1336 1773 1769 1334 10345 487 584 76775 5207 446 10243 11234 11041 28131 \
                        32604 32603 644 1305 478 1302 1764 46124 1343 46125 10372 1304 2702 83554 113107 571 1311 28035 29385 \
                        40215 194440 545 37734 1328 729 1363 582 550 732 34062 135487 37326 292 1338 1747 1710 1423 1275 207340 \
                        37923 1408 490 650 1396 29430 648 40214 1778 1502 2697049 24 38313 1353 546 1314 1303 84135 488 27973 \
                        40216 1396 38290 1379 1019 69218 173 293387 449 45634 117187 735 29382 48296 98668 28132 1463165 \
                        43770 110539 470931 654 39687 64104 1393 1389713 158836 1401 83552 1260 59814 54005 827 299767 33034 \
                        51101 33011 38304 1979527 28090 38286 38303 1276 43768 1725 38284 1727 401472 43990 61592 33028 161899 \
                        388357 1718 38289 43765 38301 43769 146827 161879 160386 169292 1705 29380 1270 258224 37637 134034 \
                        61015 28028 108486 556499 1697 187491 71999 441501 43771 29378 29379 42817 441500 156976 2559073 214473 \
                        225326 72000 71237 35755 337051 273371 333924 2055947 1891726 687363 308354 4909 687355 687361 106648 \
                        45972 68887 1134687 33964 170573 687349 246432 687350 687357 687358 687368 108980 1049583 1774 70255 337043 \
                        687352 108981 1530123 487316 813 548 83655 1288 76860 1338 257758 277944 1891727 687351 687359 687366 687346 \
                        1891762 1891763 2876790 3052189 2065049 1332244 641809 488241 335341 3050294
                # 是否有阴性对照样本？
                if [ $negative_control_exist == "TRUE" ];then   
                        # 是否已有NC样本的分类结果?没有则获取NC样本分类结果
                        if [ ! -s /opt/ossfs/mNGS_pipeline_output/NC_samples/$NC_id ];then
                                
                                cd /opt/ossfs/mNGS_pipeline_output/NC_samples
                                mkdir $NC_id
                                cd $NC_id

                                mkdir logs

                                bash /opt/mNGS/scripts/module_A.sh $NC_id 1>./logs/moduleA_$NC_id.log 2>&1
                                # 合并多个bracken 丰度表
                                /home/yushuhuan/anaconda3/bin/python /opt/softwares/Bracken/analysis_scripts/combine_bracken_outputs.py \
                                        --files $NC_id.abundance.eyePathos.standard.count $NC_id.abundance.eyePathos.pluspf.count $NC_id.abundance.eyePathos.eupathdb.count \
                                        --names standard,plusdf,eupathdb \
                                        -o $NC_id.abundance.merged.count
                        fi

                        # 提取NC样本 taxid 信息
                        taxid_in_NC=($(cat $NC_id.abundance.merged.count | awk -v FS="\t" 'NR>1{print $2}'))
                        taxid_in_NC_string=$(IFS=" "; echo "${taxid_in_NC[*]}")
                        # 移除同时在该样本及NC样本中出现的微生物
                        # --exclude参数移除对照样本中出现的背景微生物
                        cd /opt/ossfs/mNGS_pipeline_output/$sample_id
                        /home/yushuhuan/anaconda3/bin/python /opt/softwares/KrakenTools/filter_bracken.out.py -i $sample_id.abundance.eyePathos.$db.count \
                                -o $sample_id.abundance.eyePathos.$db.final.count \
                                --exclude $taxid_in_NC_string
                else
                        mv $sample_id.abundance.eyePathos.$db.count $sample_id.abundance.eyePathos.$db.final.count
                fi      
        done

        # 合并多个bracken 丰度表
        # cd /opt/ossfs/mNGS_pipeline_output/$sample_id

        if [ -s $sample_id.abundance.eyePathos.eupathdb.final.count ];then
                /home/yushuhuan/anaconda3/bin/python /opt/softwares/Bracken/analysis_scripts/combine_bracken_outputs.py \
                --files $sample_id.abundance.eyePathos.standard.final.count $sample_id.abundance.eyePathos.pluspf.final.count $sample_id.abundance.eyePathos.eupathdb.final.count \
                --names standard,plusdf,eupathdb \
                -o $sample_id.abundance.merged.count
        else
                /home/yushuhuan/anaconda3/bin/python /opt/softwares/Bracken/analysis_scripts/combine_bracken_outputs.py \
                        --files $sample_id.abundance.eyePathos.standard.final.count $sample_id.abundance.eyePathos.pluspf.final.count \
                        --names standard,plusdf,eupathdb \
                        -o $sample_id.abundance.merged.count
        fi

        # 合并多个bracken reports
        if [ -s $sample_id.abundance.eupathdb.report ];then
                /home/yushuhuan/anaconda3/bin/python /opt/softwares/KrakenTools/combine_kreports.py \
                        -r $sample_id.abundance.standard.report $sample_id.abundance.pluspf.report $sample_id.abundance.eupathdb.report \
                        -o $sample_id.raw.abundance.merged.report \
                        --only-combined --no-headers
        else
                /home/yushuhuan/anaconda3/bin/python /opt/softwares/KrakenTools/combine_kreports.py \
                                -r $sample_id.abundance.standard.report $sample_id.abundance.pluspf.report \
                                -o $sample_id.raw.abundance.merged.report \
                                --only-combined --no-headers
        fi

        # 眼科相关病原微生物丰度报出及可视化
        # bracken报告转为krona格式
        perl /opt/softwares/Krona/KronaTools/scripts/ImportTaxonomy.pl -m 3 -t 5 -n "cellular organisms" \
                -o $sample_id.raw.abundance.merged.html \
                $sample_id.raw.abundance.merged.report

        # 结果文件有2：
        # 1、$sample_id.abundance.merged.count：包含254种眼科常见病原在样本中的丰度信息等
        # 2、$sample_id.raw.abundance.merged.html：原始全部的比对结果，html格式
        # 3、$sample_id.qc_summary.html: 质控结果，可查看样本质量情况

        # 结果文件输出
        cp  $sample_id.abundance.merged.count ./results
        cp $sample_id.raw.abundance.merged.html ./results
        cp $sample_id.qc_summary.html ./results

        # 结束
        echo "Mudule A done!"
        exit 0
fi