#!/bin/bash

rm -rf /opt/mNGS/run/jobs/*

cat test_patch1_pathogen_validation_TRUE.txt | while read line
do
	sample_id=$(echo $line | awk '{print $1}')
	seqtype=$(echo $line | awk '{print $2}')

	if [ $seqtype == "cfDNA" ];then
        sample_id="$sample_id-CF"
    else
		if [ $seqtype == "RNA" ];then
			sample_id="$sample_id-R"
		fi
    fi

	echo "#!/bin/bash" >./jobs/$sample_id.onestep.run.sh

	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh

	echo "#PBS -q fz" >>./jobs/$sample_id.onestep.run.sh
	echo "#PBS -e /opt/mNGS/run/logs/$sample_id.onestep.run.err" >>./jobs/$sample_id.onestep.run.sh
	echo "#PBS -o /opt/mNGS/run/logs/$sample_id.onestep.run.log" >>./jobs/$sample_id.onestep.run.sh
	echo "#PBS -l nodes=1:ppn=8" >>./jobs/$sample_id.onestep.run.sh

	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh

	echo "start_time=\$(date "+%Y"-"%m"-"%d-%H":"%M":"%S")" >>./jobs/$sample_id.onestep.run.sh

	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh

	echo "negative_control_exist=\"FALSE\"" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh

	echo "cd /opt/ossfs/mNGS_pipeline_output" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh

	echo "mkdir $sample_id" >>./jobs/$sample_id.onestep.run.sh
	echo "cd $sample_id" >>./jobs/$sample_id.onestep.run.sh
	echo "mkdir results" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh

	echo "# 根据sample id获取原始数据所在目录，若找不到其原始数据则退出" >>./jobs/$sample_id.onestep.run.sh
	echo "sample_fastq_dir=\$(find /opt/ossfs/mNGS-data/RawData -name \"$sample_id.fq.gz\")" >>./jobs/$sample_id.onestep.run.sh
	echo "if [ ! \$sample_fastq_dir ];then" >>./jobs/$sample_id.onestep.run.sh
	echo "		echo \"Error : can't found raw fastq files for $sample_id! Please make sure it has been uploaded to server!\"" >>./jobs/$sample_id.onestep.run.sh
	echo "		exit 1" >>./jobs/$sample_id.onestep.run.sh
	echo "else" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 样本质量控制，剔除低质量序列和接头" >>./jobs/$sample_id.onestep.run.sh
	echo "		/home/yushuhuan/anaconda3/bin/fastp -i \$sample_fastq_dir \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--thread 16 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-o $sample_id.fastp.filtered.fq \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-h $sample_id.qc_summary.html \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-j $sample_id.qc_summary.json \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-n 5 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-q 15 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-u 40 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--length_required 10 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--length_limit 60 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--complexity_threshold 20 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		chmod 777 *" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 去除人源序列" >>./jobs/$sample_id.onestep.run.sh
	echo "		/usr/bin/bowtie2 -x /opt/ossfs/pathogen_reference_databases/bowtie2_index/homo_sapiens/GRCh38_noalt_as \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-U $sample_id.fastp.filtered.fq \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--un $sample_id.removedHuman.fq \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--sam-no-qname-trunc \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-S $sample_id.align.result.sam \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--threads 8" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		chmod 777 *" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 比对至standard参考库" >>./jobs/$sample_id.onestep.run.sh
	echo "		/opt/softwares/Kraken2/kraken2 --db /opt/mNGS/kraken2_index/standard \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--threads 8 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--unclassified-out unclassified_sequences.fa \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--output $sample_id.classified.standard.output \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-report $sample_id.classified.standard.report \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--report-zero-counts \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--memory-mapping \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--use-names \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				$sample_id.removedHuman.fq" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 比对至pluspf参考库" >>./jobs/$sample_id.onestep.run.sh
	echo "		/opt/softwares/Kraken2/kraken2 --db /opt/mNGS/kraken2_index/pluspf \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--threads 8 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--unclassified-out unclassified_sequences.fa \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--output $sample_id.classified.pluspf.output \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-report $sample_id.classified.pluspf.report \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--report-zero-counts \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--memory-mapping \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--use-names \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				$sample_id.removedHuman.fq" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 比对至eupathdb参考库" >>./jobs/$sample_id.onestep.run.sh
	echo "		/opt/softwares/Kraken2/kraken2 --db /opt/mNGS/kraken2_index/eupathdb \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--threads 8 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--unclassified-out unclassified_sequences.fa \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--output $sample_id.classified.eupathdb.output \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-report $sample_id.classified.eupathdb.report \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--report-zero-counts \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--memory-mapping \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--use-names \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				$sample_id.removedHuman.fq" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		chmod 777 *" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 物种丰度估计(默认在species水平上进行估计)" >>./jobs/$sample_id.onestep.run.sh
	echo "		/opt/softwares/Bracken/bracken-python3.11.3 -d /opt/mNGS/kraken2_index/standard \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-i $sample_id.classified.standard.report \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-o $sample_id.abundance.standard.count \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-w $sample_id.abundance.standard.report" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		/opt/softwares/Bracken/bracken-python3.11.3 -d /opt/mNGS/kraken2_index/pluspf \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-i $sample_id.classified.pluspf.report \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-o $sample_id.abundance.pluspf.count \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-w $sample_id.abundance.pluspf.report" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		/opt/softwares/Bracken/bracken-python3.11.3 -d /opt/mNGS/kraken2_index/eupathdb \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-i $sample_id.classified.eupathdb.report \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-o $sample_id.abundance.eupathdb.count \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-w $sample_id.abundance.eupathdb.report" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		chmod 777 *" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 原始分类结果合并" >>./jobs/$sample_id.onestep.run.sh
	echo "		files=\$(ls . | grep "count")" >>./jobs/$sample_id.onestep.run.sh
	echo "		/home/yushuhuan/anaconda3/bin/python /opt/softwares/Bracken/analysis_scripts/combine_bracken_outputs.py \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--files \$files \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				--names standard,plusdf,eupathdb \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-o $sample_id.raw.abundance.merged.count" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		chmod 777 *" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 分类结果过滤" >>./jobs/$sample_id.onestep.run.sh
	echo "		dbs=(standard pluspf eupathdb)" >>./jobs/$sample_id.onestep.run.sh
	echo "		for db in \${dbs[@]}" >>./jobs/$sample_id.onestep.run.sh
	echo "		do" >>./jobs/$sample_id.onestep.run.sh
	echo "				# 过滤掉人的分类结果" >>./jobs/$sample_id.onestep.run.sh
	echo "				/home/yushuhuan/anaconda3/bin/python /opt/softwares/KrakenTools/filter_bracken.out.py -i $sample_id.abundance.\$db.count \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						-o $sample_id.abundance.\$db.exclude.human.count \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						--exclude 9606" >>./jobs/$sample_id.onestep.run.sh
	echo "				# 提取眼科相关微生物丰度结果(第一版本中包含254种眼科常见微生物，修改--include参数，可改变最终报出的微生物丰度列表)" >>./jobs/$sample_id.onestep.run.sh
	echo "				/home/yushuhuan/anaconda3/bin/python /opt/softwares/KrakenTools/filter_bracken.out.py -i $sample_id.abundance.\$db.exclude.human.count \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						-o $sample_id.abundance.eyePathos.\$db.count \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						--include 129951 130310 130308 162425 5518 5059 5507 746128 38323 573 1352 10376 10298 10359 32603 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						10335 5476 562 139 160 5811 470 480 485 727 40324 287 29459 550 615 1282 1313 1351 28037 1280 1290 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						1283 29388 37329 1307 1292 1336 1773 1769 1334 10345 487 584 76775 5207 446 10243 11234 11041 28131 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						32604 32603 644 1305 478 1302 1764 46124 1343 46125 10372 1304 2702 83554 113107 571 1311 28035 29385 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						40215 194440 545 37734 1328 729 1363 582 550 732 34062 135487 37326 292 1338 1747 1710 1423 1275 207340 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						37923 1408 490 650 1396 29430 648 40214 1778 1502 2697049 24 38313 1353 546 1314 1303 84135 488 27973 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						40216 1396 38290 1379 1019 69218 173 293387 449 45634 117187 735 29382 48296 98668 28132 1463165 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						43770 110539 470931 654 39687 64104 1393 1389713 158836 1401 83552 1260 59814 54005 827 299767 33034 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						51101 33011 38304 1979527 28090 38286 38303 1276 43768 1725 38284 1727 401472 43990 61592 33028 161899 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						388357 1718 38289 43765 38301 43769 146827 161879 160386 169292 1705 29380 1270 258224 37637 134034 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						61015 28028 108486 556499 1697 187491 71999 441501 43771 29378 29379 42817 441500 156976 2559073 214473 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						225326 72000 71237 35755 337051 273371 333924 2055947 1891726 687363 308354 4909 687355 687361 106648 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						45972 68887 1134687 33964 170573 687349 246432 687350 687357 687358 687368 108980 1049583 1774 70255 337043 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						687352 108981 1530123 487316 813 548 83655 1288 76860 1338 257758 277944 1891727 687351 687359 687366 687346 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "						1891762 1891763 2876790 3052189 2065049 1332244 641809 488241 335341 3050294" >>./jobs/$sample_id.onestep.run.sh
	echo "				# 是否有阴性对照样本？" >>./jobs/$sample_id.onestep.run.sh
	echo "				if [ \$negative_control_exist == "TRUE" ];then"    >>./jobs/$sample_id.onestep.run.sh
	echo "						# 提取NC样本 taxid 信息" >>./jobs/$sample_id.onestep.run.sh
	echo "						taxid_in_NC=(\$(cat \$NC_id.abundance.merged.count | awk -v FS="\t" 'NR>1{print $2}'))" >>./jobs/$sample_id.onestep.run.sh
	echo "						taxid_in_NC_string=\$(IFS=" "; echo "\${taxid_in_NC[\*]}")" >>./jobs/$sample_id.onestep.run.sh
	echo "						# 移除同时在该样本及NC样本中出现的微生物" >>./jobs/$sample_id.onestep.run.sh
	echo "						# --exclude参数移除对照样本中出现的背景微生物" >>./jobs/$sample_id.onestep.run.sh
	echo "						cd /opt/ossfs/mNGS_pipeline_output/$sample_id" >>./jobs/$sample_id.onestep.run.sh
	echo "						/home/yushuhuan/anaconda3/bin/python /opt/softwares/KrakenTools/filter_bracken.out.py -i $sample_id.abundance.eyePathos.\$db.count \\" >>./jobs/$sample_id.onestep.run.sh
	echo "								-o $sample_id.abundance.eyePathos.\$db.final.count \\" >>./jobs/$sample_id.onestep.run.sh
	echo "								--exclude \$taxid_in_NC_string" >>./jobs/$sample_id.onestep.run.sh
	echo "				else" >>./jobs/$sample_id.onestep.run.sh
	echo "						mv $sample_id.abundance.eyePathos.\$db.count $sample_id.abundance.eyePathos.\$db.final.count" >>./jobs/$sample_id.onestep.run.sh
	echo "				fi" >>./jobs/$sample_id.onestep.run.sh    
	echo "			done" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		chmod 777 *" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "	 	# 合并多个bracken 丰度表" >>./jobs/$sample_id.onestep.run.sh
	echo "	 	cd /opt/ossfs/mNGS_pipeline_output/$sample_id" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		files=\$(ls . | grep "final.count")" >>./jobs/$sample_id.onestep.run.sh
	echo "		/home/yushuhuan/anaconda3/bin/python /opt/softwares/Bracken/analysis_scripts/combine_bracken_outputs.py \\" >>./jobs/$sample_id.onestep.run.sh
	echo "					--files \$files \\" >>./jobs/$sample_id.onestep.run.sh
	echo "					--names standard,plusdf,eupathdb \\" >>./jobs/$sample_id.onestep.run.sh
	echo "					-o $sample_id.abundance.merged.count" >>./jobs/$sample_id.onestep.run.sh
	echo "		chmod 777 *" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 合并多个bracken reports" >>./jobs/$sample_id.onestep.run.sh
	echo "		files=\$(ls . | grep "abundance.*.report")" >>./jobs/$sample_id.onestep.run.sh
	echo "		/home/yushuhuan/anaconda3/bin/python /opt/softwares/KrakenTools/combine_kreports.py \\" >>./jobs/$sample_id.onestep.run.sh
	echo "					-r \$files \\" >>./jobs/$sample_id.onestep.run.sh
	echo "					-o $sample_id.raw.abundance.merged.report \\" >>./jobs/$sample_id.onestep.run.sh
	echo "					--only-combined --no-headers" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		chmod 777 *" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 眼科相关病原微生物丰度报出及可视化" >>./jobs/$sample_id.onestep.run.sh
	echo "		# bracken报告转为krona格式" >>./jobs/$sample_id.onestep.run.sh
	echo "		perl /opt/softwares/Krona/KronaTools/scripts/ImportTaxonomy.pl -m 3 -t 5 \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				-o $sample_id.raw.abundance.merged.html \\" >>./jobs/$sample_id.onestep.run.sh
	echo "				$sample_id.raw.abundance.merged.report" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		chmod 777 *" >>./jobs/$sample_id.onestep.run.sh
	echo -e "\n" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 结果文件输出" >>./jobs/$sample_id.onestep.run.sh
	echo "		cp  $sample_id.abundance.merged.count ./results" >>./jobs/$sample_id.onestep.run.sh
	echo "		cp $sample_id.raw.abundance.merged.html ./results" >>./jobs/$sample_id.onestep.run.sh
	echo "		cp $sample_id.qc_summary.html ./results" >>./jobs/$sample_id.onestep.run.sh

	echo "		# 结束" >>./jobs/$sample_id.onestep.run.sh
	echo "		end_time=\$(date "+%Y"-"%m"-"%d-%H":"%M":"%S")" >>./jobs/$sample_id.onestep.run.sh
	echo "		# 分析完成" >>./jobs/$sample_id.onestep.run.sh
	echo "		echo "开始处理的时间为：" \$start_time" >>./jobs/$sample_id.onestep.run.sh
	echo "		echo "处理结束的时间为：" \$end_time" >>./jobs/$sample_id.onestep.run.sh
	echo "		exit 0" >>./jobs/$sample_id.onestep.run.sh
	echo "fi" >>./jobs/$sample_id.onestep.run.sh

	# sed 's/\$sample_id/'$sample_id'/g' one-step-analysis.sh >./jobs/$sample_id.onestep.run.sh
	# sed -i 's/DNA/'$seqtype'/g' ./jobs/$sample_id.onestep.run.sh
done