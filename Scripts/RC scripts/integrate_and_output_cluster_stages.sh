#!/bin/bash 
#SBATCH -J zf6s_Integrate
#SBATCH -n 6 # Number of cores 
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 0-05:00 # Runtime in D-HH:MM 
#SBATCH -p serial_requeue # Partition to submit to 
#SBATCH --mem=10000 # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o run_log/Integrate_6s.out # File to which STDOUT will be written 
#SBATCH -e run_log/Integrate_6s.err # File to which STDERR will be written 

module load Anaconda3/5.0.1-fasrc02

#stage="zf20hpf"
stage="zf6s"

## Integrate objects from multiple jobs
echo "Integrating results into a single .pkl file..."
script_dir="./"
#in_dir="../Results/${stage}/"
in_dir="../Results/testCentos7/"
out_dir="${in_dir}/integrated/"
py_script="${script_dir}/integrate_objs.py"

mkdir -p ${out_dir}
python -u ${py_script} -i ${in_dir} -o ${out_dir} -perm True -analyze True -prt_top True

## Write out tables of result matrices
#script_dir="./"
echo "Writing result tables as .csv..."
in_f="${out_dir}/All_results_obj"
tbl_dir="${out_dir}/tables/"
py_script="${script_dir}/output_results.py"

mkdir -p ${tbl_dir}
python -u ${py_script} -i ${in_f} -o ${tbl_dir} -K None -prt_top True

module load R/3.4.2-fasrc01

## Organize result tables into R object
echo "Collecting tables into an R object..."
R_script="${script_dir}/collect_tbls.R"
Rscript ${R_script} ${tbl_dir} 
echo "Done!"
