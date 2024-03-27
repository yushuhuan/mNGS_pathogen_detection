#!/bin/bash

#PBS -q fz
#PBS -e /opt/mNGS/run/logs/$sample_id.onestep.run.err
#PBS -o /opt/mNGS/run/logs/$sample_id.onestep.run.log
#PBS -l nodes=1:ppn=8


# 获取开始处理的时间
start_time=$(date "+%Y"-"%m"-"%d-%H":"%M":"%S")

negative_control_exist="FALSE"
seqtype="DNA"

#echo "Processing " $sample_id

cd /opt/ossfs/mNGS_pipeline_output

if [ -s /opt/ossfs/mNGS_pipeline_output/$sample_id/results/$sample_id.abundance.merged.count ];then
        echo "Sample has been analysised before, just skip it."
	exit 0
else
	if [ ! -d /opt/ossfs/mNGS_pipeline_output/$sample_id ];then
        	mkdir $sample_id
	fi
		
	cd $sample_id
	
	if [ ! -d /opt/ossfs/mNGS_pipeline_output/$sample_id/results ];then
	        mkdir results
        	mkdir logs
	fi

	# 如果测序样本是cfDNA，则样本名为”样本id-CF“
        if [ $seqtype == "cfDNA" ];then
                newSample_id="$sample_id-CF"
                bash /opt/mNGS/run/module_A.sh $newSample_id 1>./logs/moduleA_$sample_id.log 2>&1
        else
                if [ $seqtype == "RNA" ];then
                        newSample_id="$sample_id-R"
                        bash /opt/mNGS/run/module_A.sh $newSample_id 1>./logs/moduleA_$sample_id.log 2>&1
                else
                        if [ $seqtype == "DNA" ];then
                        # 执行module A 分析
                                bash /opt/mNGS/run/module_A.sh $sample_id 1>./logs/moduleA_$sample_id.log 2>&1
                        fi
                fi
        fi
fi

# 获取处理结束的时间
end_time=$(date "+%Y"-"%m"-"%d-%H":"%M":"%S")

# 分析完成
echo "开始处理的时间为：" $start_time
echo "处理结束的时间为：" $end_time

echo "Done!"

exit 0
