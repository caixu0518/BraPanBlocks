for i in {0..10}
    do
    for j in {1..10}
        do
        echo   treemix -i   ECD04.snps.maf0.05.fm.sub.clean.pruned.input_for_treemix.frq.txt.gz    -root Rapini -k 500 -m ${i} -bootstrap -o migration_${i}_bootstrap_${j}
    done
done


