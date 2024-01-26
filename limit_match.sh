# 1 memory limit
data_path=/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin
dataset=friendster_csr
group=GraphPi
mem_limit=$1
query=$2
optimizations=$3
threads=$4
threshold=$5
density=$6
interval=$7
sudo cgset -r memory.limit_in_bytes=${mem_limit} ${group}
sudo cgset -r memory.memsw.limit_in_bytes=${mem_limit} ${group}
sudo cgset -r memory.limit_in_bytes=${mem_limit} ${group}

sudo cgexec -g memory:$group perf stat -e minor-faults,major-faults,cycles,instructions,block:block_rq_issue,block:block_rq_complete bin/pattern_match ${data_path}/$dataset $query $optimizations $threads $threshold $density $interval > ${dataset}_${mem_limit}_O${optimizations}_T${threads}_Q${query}_R${threshold}_D${density}_I${interval}.log 2>&1 &

