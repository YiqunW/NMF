#!/usr/bin/python
import argparse
import sys
sys.path.append("../")
import nmf_fxn
import ast
import copy
import pickle
import pandas as pd
import numpy as np
import os

### Take arguments from command line --- 
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-in_f','-i', help='path of the Results_obj object', required=True)
    parser.add_argument('-out_dir','-o', help='output directory for the tables', required=True)
    parser.add_argument('-num_k','-K', help='out put results for which K values', default=None)
    parser.add_argument('-prt_top_genes','-prt_top', help='whether to save tables of 30 top ranking genes and their weights', default=True)
        
        
    args = parser.parse_args()
    in_f = args.in_f
    out_dir = args.out_dir
    ks = ast.literal_eval(args.num_k)
    try:
        prt_top_genes = ast.literal_eval(args.prt_top_genes)
    except:
        prt_top_genes = args.prt_top_genes

## Example run:
##    python ./run_nmf.py -i expr.txt -K [3,40] -rep 30 -o ./results -ngene 35
##### write a python script that runs nmf and prints out some results #####
## read in expression matrix

print("Loading Results_obj...")
Results_all=nmf_fxn.load_obj(in_f)

print("Writting out parameters used...")
Results_all.params2txt(out_dir+"/parameters.txt")

if Results_all.data["scaled"] is not None:
    print("Writting out scaled dataset used in NMF...")
    Results_all.data["scaled"].to_csv(out_dir+"/scaled_data.csv")
## Maybe output permuted as well?

if ks is None:
    ks=Results_all.Params.Ks
else:
    if type(ks) == int:
        ks = [ks]
    elif len(ks) == 2:
        ks = range(ks[0],ks[1]+1)
    elif len(ks) > 2:
        ks = ks
    else:
        print('Invalid K input. Must be Noe, an integer, or a list of integers.')

print("Outputting results for the following K values:")
print(list(ks))

Results_=Results_all.results
for k in ks:
    results_k=Results_["K="+str(k)]
    print("Writing tables for K="+str(k)+"...")
    dir_k=out_dir+"/K="+str(k)
    if not os.path.exists(dir_k):
        os.makedirs(dir_k)
    for rep in results_k.keys():
        dir_rep=dir_k+"/"+rep
        if not os.path.exists(dir_rep):
            os.makedirs(dir_rep)
        results_k[rep]["G"].to_csv(dir_rep+"/G.csv")
        results_k[rep]["C"].to_csv(dir_rep+"/C.csv")
        if prt_top_genes:
            nmf_fxn.print_top_genes(results_k[rep]["G"], n_top_genes = 30,component = None, prt= False, save_tbl=dir_rep+"/top30genes.csv")

print("Done.")