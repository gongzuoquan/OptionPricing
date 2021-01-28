#!/bin/sh
dir="./dir_prof"
exe="./execute"
time=$(date "+%Y-%m-%d-%H:%M:%S")

if [ ! -d $prof ];then
    mkdir $prof
fi
nvprof --export-profile $dir/timeline-$time.prof $exe
nvprof --metrics achieved_occupancy,executed_ipc -o $dir/metrics-$time.prof $exe
