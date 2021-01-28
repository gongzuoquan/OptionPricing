t=4096
for i in {1..6};do
    echo ${t}
    ./execute $t 300
    t=$[$t*4]
    #echo ${t}
done

