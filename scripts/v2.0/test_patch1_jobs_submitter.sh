#!/bin/bash

# n=0
cat undone.list | while read line
do
    sample_id=$(echo $line | awk '{print $1}')

    # 控制提交的任务数，每4个任务为一组，每组的最后一个样本提交后过15min再提交
    n=$((n+1))
    if [ $((n%4)) -eq 0 ];then
        sleep 40m
    else
        qsub /opt/mNGS/run/jobs/$sample_id.onestep.run.sh
    fi
done
