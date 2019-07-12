#!/bin/bash 
#SBATCH -J NMF_zf2dpf
#SBATCH --array=36,37,38,39%4   ## --array=X-Y%Z, where X-Y is the range of job
#SBATCH -n 5 # Number of cores 
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 0-20:00 # Runtime in D-HH:MM 
#SBATCH -p serial_requeue # Partition to submit to 
#SBATCH --mem=20000 # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o run_log/NMF_zf2dpf_%a.out # File to which STDOUT will be written 
#SBATCH -e run_log/NMF_zf2dpf_%a.err # File to which STDERR will be written 

## Modified on May10 2018 to run combined stages in order to find seven sleeper module in bud and dome stages. No median scaling was performed. No permuted data was decomposed in paralelle
source new-modules.sh
module load Anaconda3/2.1.0-fasrc01
source activate python_env1 ## sklearn 0.17 and seaborn are installed in this environment

data_dir="../Data/all_var/"

#"zf6s"
#"zf10s"
#"zf14s"
#"zf18s"
#"zf20hpf"
#"zf24hpf"
#"zf2dpf"

#stage="zf18s"
stage="zf2dpf"
py_script="./run_nmf.py"
expr_f="${data_dir}${stage}_vargene.csv"
perm_f="${data_dir}${stage}_vargene_perm.csv"
out_dir="../Results/${stage}/"
mkdir -p ${out_dir}
#python -u ${py_script} -i ${expr_f} -K ${SLURM_ARRAY_TASK_ID} -rep 5 -o ${out_dir} -scl False -miter 10000 -perm ${perm_f} -run_perm True -tol 1e-7 -a 5 -init "nndsvd"
#python -u ${py_script} -i ${expr_f} -K ${SLURM_ARRAY_TASK_ID} -rep 5 -o ${out_dir} -scl False -miter 10000 -perm ${perm_f} -run_perm True -tol 1e-7 -a 2 -init "nndsvd"
python -u ${py_script} -i ${expr_f} -K ${SLURM_ARRAY_TASK_ID} -rep 5 -o ${out_dir} -scl "median" -miter 10000 -perm ${perm_f} -run_perm True -tol 1e-7 -a 2 -init "nndsvd"

source deactivate
