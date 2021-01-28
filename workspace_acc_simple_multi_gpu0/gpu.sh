t=4096
for i in {1..10};do
    ./execute_multiGPU $t 500
    t=$[$t*4]
    #echo ${t}
done

