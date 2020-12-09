
export top_dir=/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/1fft/10cutoff
export analys_dir=$top_dir/analysis
export mod_dir=$top_dir/analysis/GSMs_cl1/
export name=1fft

# cp $top_dir/density.txt $mod_dir
# cp $mod_dir/sample_A.txt $mod_dir/Scores_A.txt
# cp $mod_dir/sample_B.txt $mod_dir/Scores_B.txt

#cp analysis3GSMs_cl/sample_A.txt Scores_A.txt
#cp analysis3GSMs_cl/sample_B.txt Scores_B.txt

nohup python ~/newIMPanalysis/Instruct_Workshop/analysis/scripts/Master_Sampling_Exhaustiveness_Analysis.py --sysname $name --path $mod_dir --extension rmf --mode cpu_omp --cores 20 --align --density density.txt --gridsize 1.0 --gnuplot  --scoreA sample_A.txt --scoreB sample_B.txt> clustering.log &

