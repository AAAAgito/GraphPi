# 1 memory limit
data_path=/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin
dataset=friendster_csr
group=GraphPi
mem_limit=$1
query=$2
optimizations=$3
threads=$4
sudo cgset -r memory.limit_in_bytes=${mem_limit} ${group}
sudo cgset -r memory.memsw.limit_in_bytes=${mem_limit} ${group}
sudo cgset -r memory.limit_in_bytes=${mem_limit} ${group}

sudo cgexec -g memory:$group perf stat -e minor-faults,major-faults,cycles,instructions bin/pattern_match ${data_path}/$dataset $query $optimizations $threads > perf_cycles${dataset}_${mem_limit}_prefetch.log 2>&1 &

