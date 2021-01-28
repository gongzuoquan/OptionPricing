#!/bin/sh

max_scale=6
if [ $# >1 ];then
    max_scale=$1
fi

project_path=$(basename `pwd`)
suffix=${project_path:9}
gtime=$(date "+%Y-%m-%d-%H-%M-%S")
log=log$gtime-$suffix
for file in $(ls)
do
    if [[ "$file" = log* ]];then
        mv $file backup
    fi
done

t=4096
for i in $(seq 1 $max_scale);do
    echo "计算规模：" ${t}
    echo "计算时间:"
    echo ${t} >log
    time ./execute $t 300 >log
    echo ""


    t=$[$t*4] 
done
#echo $suffix
