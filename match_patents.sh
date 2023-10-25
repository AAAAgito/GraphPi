export base=/home/yanglaoyuan/AsyncSubGraphStorage/dataset/bin/
export threads=40

g=(patents live_journal orkut friendster)
tri=(7515023 177820130 627584181 4173724142)

export experiment=0
for query in {1..4}
do
    echo q$query
    for setting in 0
    do
    i=3
    echo "set "$setting
    
        export graph=$base${g[i]}
        export tricnt=${tri[i]}
        echo "$graph"
        echo ""
        syrupy.py bin/graph_client $graph $threads $experiment $setting $query $tricnt
    
    done
done