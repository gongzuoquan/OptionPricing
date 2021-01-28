#!/bin/sh
# 用于批量log转化csv
dir=("GPU1" "GPU2" "CPUs" "CPU")
base="workspace"
for d in ${dir[@]}
do
    argv=$base$d/log*${d,,}
    python log2csv.py $argv
done
target=""
for d in ${dir[@]}
do
    #target=$target$base$d/log*${d,,}.csv
    target=$target$base$d/log*$d.csv
    target=$target" "
done
echo $target
zip all_csv.zip -j $target
