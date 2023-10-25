export input=/home/yanglaoyuan/AsyncSubGraphStorage/WebGraph_Trans/docker-webgraph/gsh-2015.edge
export output=/home/yanglaoyuan/AsyncSubGraphStorage/dataset/large_graph/gsh-2015
export v=988490691
export mem=32000000000

# syrupy.py bin/export_large_graph $input $output $v $mem
bin/export_large_graph $input $output $v $mem