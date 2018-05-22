#!/bin/bash 
script_dir="./" # directory of run_nmf.py
data_dir="./" # directory of input data
data_file="example_data.csv" # file name of input data
out_dir="./example_results/" # output directory
py_script="${script_dir}/run_nmf.py"

expr_f="${data_dir}${data_file}"
krange="[5,10]" # K (rank) for running NMF. ex. "[5]": NMF will be performed for K=5; "[3,5]": NMF will be performed for K=3, K=4, and K=5; "[3,5,6,8]": NMF will be performed for K=3, K=5, K=6 and K=8. 

mkdir -p ${out_dir} # create output directory if doesn't exist already

start=`date +%s`
python -u ${py_script} -i ${expr_f} -K ${krange} -rep 5 -o ${out_dir} -scl "median" -miter 10000 -perm True -run_perm True -tol 1e-6 -a 2 -init "nndsvd"
# see run_nmf.py for explanation of parameters
end=`date +%s`

runtime=$((end-start))
echo "Run time: ${runtime}s"