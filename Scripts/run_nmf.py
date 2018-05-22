#!/usr/bin/python
import argparse
import sys
sys.path.append("../")
import nmf_fxn
import ast
import copy
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import pickle
import pandas as pd
import numpy as np
import os
import seaborn as sns

### Take arguments from command line --- 
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-in_f','-i', help='expr.txt input file path', required=True)
    parser.add_argument('-out_dir','-o', help='output result dictionaries and picture file directory', required=True)
    parser.add_argument('-scale','-scl', help='whether to scale the expression matrix. Can be False, "median" or "max"', default=False)
    parser.add_argument('-num_k','-K', help='number of cluster/module/reduced dimension', required=True)
    parser.add_argument('-rand_state','-rand', help='random state(s) for running nmf (values other than default not tested)', default=None)
    parser.add_argument('-alpha','-a', help='alpha value for running nmf', default=0.25)
    parser.add_argument('-l1','-l1', help='l1 ratio for running nmf', default=0.5)
    parser.add_argument('-max_iter','-miter', help='maximum iteratiion for running nmf', default=25000)
    parser.add_argument('-tol','-tol', help='error tolerance for running nmf', default=1e-7)
    parser.add_argument('-num_rep','-rep', help='number of repeated runs of NMF with random initial conditions for each K', default=1)
    parser.add_argument('-sub_spl','-sub', help='number of cells for running NMF on a random subset of cells (values other than default not tested)', default=None)
    parser.add_argument('-permute','-perm', help='whether to generate a randomly permuted matrix as control', default=True)
    parser.add_argument('-prt_top_genes','-prt_top', help='whether to save tables of 30 top ranking genes and their weights', default=True)
    parser.add_argument('-init','-init', help='Initialization method for NMF. Default is nndsvdar (faster than nndsvda, which is better when sparseness is not desired). Could also use nndsvd (better for sparseness), or random.', default=None)
    parser.add_argument('-run_perm','-run_perm', help='Whether to perform NMF on a permuted dataset. Accepted arguments are True, False or "only"', default=False)
    parser.add_argument('-analyze','-analyze', help='Whether to analyze error, n_iter, and stabilities of the results.', default=True)
    
        
    args = parser.parse_args()
    in_file = args.in_f
    init = args.init
    try:
        scale = ast.literal_eval(args.scale)
    except:
        scale = args.scale
    ks = ast.literal_eval(args.num_k) ## k values provided in ks argument need to be in increasing order
    try:
        rand_state = ast.literal_eval(args.rand_state)
    except:
        rand_state = args.rand_state
    try:
        alpha = ast.literal_eval(args.alpha)
    except:
        alpha = args.alpha
    try:
        l1 = ast.literal_eval(args.l1)
    except:
        l1 = args.l1
    try:
        max_iter = ast.literal_eval(args.max_iter)
    except:
        max_iter = args.max_iter
    try:
        tol = ast.literal_eval(args.tol)
    except:
        tol = args.tol
    try:
        rep = ast.literal_eval(args.num_rep)
    except:
        rep = args.num_rep
    try:
        sub = ast.literal_eval(args.sub_spl)
    except:
        sub = args.sub_spl
    try:
        perm = ast.literal_eval(args.permute)
    except:
        perm = args.permute
    out_dir = args.out_dir
    try:
        prt_top_genes = ast.literal_eval(args.prt_top_genes)
    except:
        prt_top_genes = args.prt_top_genes
    try:
        run_perm = ast.literal_eval(args.run_perm)
    except:
        run_perm = args.run_perm
    try:
        analyze = ast.literal_eval(args.analyze)
    except:
        analyze = args.analyze

## Example run:
##    python ./run_nmf.py -i expr.txt -K [3,40] -rep 30 -o ./results -ngene 35
##### write a python script that runs nmf and prints out some results #####
## read in expression matrix

run_nmf=nmf_fxn.nmf_reps()

if type(ks) == int:
    ks = [ks]
elif len(ks) == 2:
    ks = range(ks[0],ks[1]+1)
elif len(ks) > 2:
    ks = ks
else:
    print('Invalid K input. Must be an integer or a list of integers.')
print("K values determined:")
print(list(ks))

print("Setting parameters for class nmf_reps...")
run_nmf.set_param(scale=scale, Ks=ks, rand_state=rand_state, alpha=alpha, l1=l1, max_iter=max_iter, tol=tol, rep=rep, sub=sub, init=init,verbose=True, permute=perm)

data = nmf_fxn.read_data(in_file)
if data.shape[0]*data.shape[1]==0:
    data=pd.read_csv(in_file,index_col=0)
print("Expression data read. Dataset dimensions:")
print(data.shape)

print("Setting datasets for class nmf_reps...")
run_nmf.set_data(data)

print("Running nmf...")
if prt_top_genes:
    prt_top_genes=out_dir+"/Top30genes"
run_nmf.nmf_results(permuted=run_perm,prt_top_genes=prt_top_genes)
if run_perm=="only":
    print("NMF runs finished. Results are stored in run_nmf.permuted_results.")
elif run_perm:
    print("NMF runs finished. Results are stored in run_nmf.results and run_nmf.permuted_results.")
else:
    print("NMF runs finished. Results are stored in run_nmf.results.")

file_name = out_dir+"/Results_obj.pkl"
expand=0
if os.path.isfile(file_name):
    expand = 1
    while True:
        expand += 1
        new_file_name = file_name.split(".pkl")[0] + "(" + str(expand) + ")" + ".pkl"
        if os.path.isfile(new_file_name):
            continue
        else:
            file_name = new_file_name
            break

print("Pickling the result object "+file_name+"...")
nmf_fxn.save_obj(run_nmf,file_name)

if analyze:
    if expand==0:
        post=""
    else:
        post="(" + str(expand) + ")"
    if run_perm is not "only":
        print("Plotting reconstruction errors and actual iteration numbers...")
        run_nmf.plot_err(stat="err",measure="mean",save=out_dir+"/Reconstruction_errors"+post)
        run_nmf.plot_err(stat="n_iter",measure="mean",save=out_dir+"/Actual_iterations"+post)
        print("Plotting ratios between the 2 top ranking genes...")
        ratio_tbl=run_nmf.top_ratio()
        pp = PdfPages(out_dir+'/TopWeightsRatio'+post+'.pdf')
        plt.figure()
        sns.swarmplot(x="K", y="ratio(rank1/rank2)", data=ratio_tbl,size=3)
        pp.savefig()
        if run_perm:
            perm_ratio=run_nmf.top_ratio(permuted=True)
            plt.figure()
            sns.swarmplot(x="K", y="ratio(rank1/rank2)", data=perm_ratio,size=3)
            plt.title("Results from Permuted Dataset")
            pp.savefig()
        pp.close()
    print("Calculating result stability measurements...")
    run_nmf.calc_stability(stats=["inconsistency_G","inconsistency_C","cophenetic_C","cophenetic_G"],permuted=run_perm)

    if run_perm is not "only":
        print("Plotting result stability...")
        run_nmf.plot_stability(permuted=None,stats=["inconsistency_G","inconsistency_C","cophenetic_C","cophenetic_G"],save=out_dir+"/Stability"+post,ttl=False,fg_sz=[7,5])

    print("Pickling the result object with stability table(s)...")
    nmf_fxn.save_obj(run_nmf,file_name)

print("All done.")

