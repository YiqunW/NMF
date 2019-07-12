#!/bin/bash 
#SBATCH -J NMF_integrate
#SBATCH -n 6 # Number of cores 
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 0-05:00 # Runtime in D-HH:MM 
#SBATCH -p serial_requeue # Partition to submit to 
#SBATCH --mem=16000 # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o run_log/Integrate_SWTallvar.out # File to which STDOUT will be written 
#SBATCH -e run_log/Integrate_SWTallvar.err # File to which STDERR will be written 

source new-modules.sh
module load Anaconda3/2.1.0-fasrc01
source activate python_env1

## Integrate objects from multiple jobs
echo "Integrating results into a single .pkl file..."
script_dir="./"
in_dir="../Results/ZFSWT_all/allvar_Parameters4_noScale/"
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

source new-modules.sh
module load R/3.4.1-fasrc01
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER

## Organize result tables into R object
echo "Collecting tables into an R object..."
R_script="${script_dir}/collect_tbls.R"
Rscript ${R_script} ${tbl_dir} 
echo "Done!"
