export base=/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/
export threads=40

g=(patents live_journal orkut friendster)
tri=(7515023 177820130 627584181 4173724142)

# echo "set "$setting
# setting=2
# graph=/home/yanglaoyuan/AsyncSubGraphStorage/dataset/large_graph/gsh-2015
# query=1
# bin/graph_client $graph $threads $experiment $setting $query 1


# bash export_graph.sh patents
# bash export_graph.sh live_journal
# bash export_graph.sh orkut
# bash export_graph.sh friendster

# export experiment=0
# for query in 4
# do
#     echo q$query
#     for setting in 0 1
#     do
#         echo "set "$setting
#         for ((i=0;i<4;i++))
#         do
#             export graph=$base${g[i]}
#             export tricnt=${tri[i]}
#             echo "$graph"
#             echo ""
#             bin/graph_client $graph $threads $experiment $setting $query $tricnt
#         done
#     done
# done

export experiment=1
export rate=1
export split=0
for rate in 1 2 3
do
echo "=====   $rate   ======="
for setting in 0 1
do
    echo "set "$setting
    
    for ((i=0;i<4;i++))
    do
        export graph=$base${g[i]}
        # bash export_graph.sh ${g[i]}
        # cp ${graph}_dup $graph
        echo ${g[i]}
        bin/graph_client $graph $threads $experiment $setting $rate $split
        bash export_graph.sh ${g[i]}
    done
done
done
# bash export_graph.sh patents
# bash export_graph.sh live_journal
# bash export_graph.sh orkut
# bash export_graph.sh friendster