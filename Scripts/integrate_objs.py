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
    parser.add_argument('-in_dir','-i', help='directory for the Results_obj objects', required=True)
    parser.add_argument('-out_dir','-o', help='output directory for All_results_obj and figures', required=True)
    parser.add_argument('-analyze','-analyze', help='Whether to analyze error, n_iter, and stabilities of the results', default=True)
    parser.add_argument('-permute','-perm', help='whether to analyze results from the permuted matrix. Considered only if analyze is set to True', default=False)
    parser.add_argument('-prt_top_genes','-prt_top', help='whether to save tables of 30 top ranking genes and their weights. Considered only if analyze is set to True', default=False)
        
        
    args = parser.parse_args()
    in_dir = args.in_dir
    
    try:
        run_perm = ast.literal_eval(args.permute)
    except:
        run_perm = args.permute
    out_dir = args.out_dir
    try:
        prt_top_genes = ast.literal_eval(args.prt_top_genes)
    except:
        prt_top_genes = args.prt_top_genes
    try:
        analyze = ast.literal_eval(args.analyze)
    except:
        analyze = args.analyze

## Example run:
##    python ./run_nmf.py -i expr.txt -K [3,40] -rep 30 -o ./results -ngene 35
##### write a python script that runs nmf and prints out some results #####
## read in expression matrix

print("Looking for Results_obj in out_dir...") ## look for all files with .pkl extensions
res_objs=[result_obj for result_obj in os.listdir(in_dir) if result_obj.endswith(".pkl")]

## load in the first object in the res_objs list
Results_all=nmf_fxn.load_obj(in_dir+"/"+res_objs[0])

if len(res_objs)>1:
    ## see if the loaded object is empty (without a "results" attribute). If so, skip this object and read in the next until encountering a non-empty object.
    i=0
    while not hasattr(Results_all,"results"):
        pritn(res_objs[i]+" doesn't contain any results.")
        i=i+1
        Results_all=nmf_fxn.load_obj(in_dir+"/"+res_objs[i])

    ## if the object doesn't contain results for all the K values specified in .Params.Ks, skip it and read in the next object in the list
    while len(Results_all.results.keys())!=len(list(Results_all.Params.Ks)):
        pritn(res_objs[i]+" doesn't contain all the results matching to its Ks parameters.")
        i=i+1
        Results_all=nmf_fxn.load_obj(in_dir+"/"+res_objs[i])
        ## check if the new object has a results attribute, if not, read in the next object in the list
        while not hasattr(Results_all,"results"):
            pritn(res_objs[i]+" doesn't contain any results.")
            i=i+1
            Results_all=nmf_fxn.load_obj(in_dir+"/"+res_objs[i])

    params=list(Results_all.Params.keys())
    datasets=list(Results_all.data.keys())

    res_objs=res_objs[i:]

    if len(res_objs)>1:
        for result_obj in res_objs[1:]:
            ## load in the next object in the res_objs list
            result_add=nmf_fxn.load_obj(in_dir+"/"+result_obj)
            ## check if the object has results attribute
            if not hasattr(result_add,"results"):
                print(result_obj+" doesn't contain any results.")
            else:
                ## check if the object has results for each K specified in its Params.Ks
                if len(result_add.results.keys())!=len(list(result_add.Params.Ks)):
                    print(result_obj+" doesn't contain all the results matching to its Ks parameters.")
                else:
                    ## check if parameters agree and add in new parameters (such as new K values in .Params.Ks)
                    for par in params:
                        if par=="Ks":
                            Ks=list(set.union(set(Results_all.Params.Ks),set(result_add.Params.Ks)))
                            Results_all.Params.Ks=sorted(Ks)
                        elif par=="rand_state":
                            try:
                                rand_states=list(Results_all.Params.rand_state)+list(result_add.Params.rand_state)
                                Results_all.Params.rand_state=rand_states
                            except TypeError:
                                pass
                        elif par=="rep":
                            if Results_all.Params.rep != result_add.Params.rep:
                                if set(Results_all.Params.Ks) == set(result_add.Params.Ks):
                                    Results_all.Params.rep=Results_all.Params.rep+result_add.Params.rep
                                else:
                                    Results_all.Params.rep=None ## should I change rep to an array or a list to match each K?
                        else:
                            if Results_all.Params[par] != result_add.Params[par]: ## for these parameters, it makes less sence to change them to a list that matches Ks.
                                print("Error: parameters don't agreen between "+res_objs[0]+" and "+result_obj+". Cannot integrate Results.")
                                exit()
                    ## check if datasets agree
                    for dataset in datasets:
                        if type(Results_all.data[dataset]) != type(result_add.data[dataset]): ## should change such that if some objects don't have a permuted dataset, ignore it and integrate the objects with real data.
                            print("Error: datasets don't agreen between "+res_objs[0]+" and "+result_obj+". Cannot integrate Results.")
                            exit()
                        else:
                            if Results_all.data[dataset] is not None:
                                if np.array(Results_all.data[dataset]-result_add.data[dataset]).any():
                                    if dataset!="permuted":
                                        print("Error: datasets don't agreen between "+res_objs[0]+" and "+result_obj+". Cannot integrate Results.")
                                        exit()
                    ## if all agree, integrate results
                    if hasattr(result_add,"permuted_results"):
                        add_permuted=True
                        if not hasattr(Results_all,"permuted_results"):
                            Results_all.permuted_results={}
                    else:
                        add_permuted=False

                    ## adding results from result_add to results in Results_all
                    for res_k in result_add.results.keys():
                        if res_k in Results_all.results.keys(): ## if there are results for the same K in the two objects, append the repeated runs in the new object to the ones in Results_all
                            rep_st=1+max([int(rep.split("p")[1]) for rep in Results_all.results[res_k].keys() if rep.startswith('rep')])
                            for rep in result_add.results[res_k].keys():
                                Results_all.results[res_k]["rep"+str(rep_st)]=copy.deepcopy(result_add.results[res_k][rep])
                                rep_st+=1
                        else:
                            Results_all.results[res_k]=copy.deepcopy(result_add.results[res_k])
                    if add_permuted:
                        for res_k in result_add.permuted_results.keys():
                            if res_k in Results_all.permuted_results.keys():
                                rep_st=1+max([int(rep.split("p")[1]) for rep in Results_all.permuted_results[res_k].keys() if rep.startswith('rep')])
                                for rep in result_add.permuted_results[res_k].keys():
                                    Results_all.permuted_results[res_k]["rep"+str(rep_st)]=copy.deepcopy(result_add.permuted_results[res_k][rep])
                                    rep_st+=1
                            else:
                                Results_all.permuted_results[res_k]=copy.deepcopy(result_add.permuted_results[res_k])

    file_name = out_dir+"/All_results_obj.pkl"
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

    print("Pickling the integrated result object"+file_name+"...")
    nmf_fxn.save_obj(Results_all,file_name)
else:
    file_name = out_dir+"/Analyzed_results_obj.pkl"
    expand=0

if analyze:
    if expand==0:
        post=""
    else:
        post="(" + str(expand) + ")"
    if run_perm:
        if not hasattr(Results_all,"permuted_results"):
            run_perm=False
    if run_perm is not "only":
        if prt_top_genes:
            print("Saving tables for 30 top ranking genes and their weights...")
            for k in Results_all.Params.Ks:
                nmf_fxn.print_top_genes(Results_all.results["K="+str(k)]["rep0"]["G"], n_top_genes = 30,component = None, prt= False, save_tbl=out_dir+"/Top30genes_K"+str(k)+"_rep0")
                if run_perm:
                    nmf_fxn.print_top_genes(Results_all.permuted_results["K="+str(k)]["rep0"]["G"], n_top_genes = 30,component = None, prt= False, save_tbl=out_dir+"/Top30genes_permuted_K"+str(k)+"_rep0")
        print("Plotting reconstruction errors and actual iteration numbers...")
        Results_all.plot_err(stat="err",measure="mean",save=out_dir+"/Reconstruction_errors"+post)
        Results_all.plot_err(stat="n_iter",measure="mean",save=out_dir+"/Actual_iterations"+post)
        print("Plotting ratios between the 2 top ranking genes...")
        ratio_tbl=Results_all.top_ratio()
        pp = PdfPages(out_dir+'/TopWeightsRatio'+post+'.pdf')
        plt.figure()
        sns.swarmplot(x="K", y="ratio(rank1/rank2)", data=ratio_tbl,size=3)
        pp.savefig()
        if run_perm:
            perm_ratio=Results_all.top_ratio(permuted=True)
            plt.figure()
            sns.swarmplot(x="K", y="ratio(rank1/rank2)", data=perm_ratio,size=3)
            plt.title("Results from Permuted Dataset")
            pp.savefig()
        pp.close()

    print("Calculating result stability measurements...")
    Results_all.calc_stability(stats=["inconsistency_G","inconsistency_C","cophenetic_C","cophenetic_G"],permuted=run_perm)

    if run_perm is not "only":
        print("Plotting result stability...")
        Results_all.plot_stability(permuted=None,stats=["inconsistency_G","inconsistency_C","cophenetic_C","cophenetic_G"],save=out_dir+"/Stability"+post,ttl=False,fg_sz=[7,5])


    print("Pickling the result object with stability table(s)...")
    nmf_fxn.save_obj(Results_all,file_name)

print("All done.")