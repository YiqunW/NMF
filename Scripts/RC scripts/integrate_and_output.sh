#!/bin/bash 
## Integrate objects from multiple jobs
echo "Integrating results into a single .pkl file..."
script_dir="./"
in_dir="../example_results_new/"
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

## Organize result tables into R object
echo "Collecting tables into an R object..."
R_script="${script_dir}/collect_tbls.R"
Rscript ${R_script} ${out_dir} 
echo "Done!"