export graph=patents
export input_path=/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/$graph
export ins_r=0.1
export del_r=0.1
export split=5
bin/generate_update $input_path $ins_r $del_r $split 20