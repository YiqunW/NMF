## import necessary libraries/packages
import numpy as np
import math
from sklearn.preprocessing import MaxAbsScaler
import matplotlib
matplotlib.use('Agg')
import pylab as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.cluster import AgglomerativeClustering
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from sklearn.cluster import DBSCAN
from sklearn.decomposition import NMF, LatentDirichletAllocation
import pandas as pd
import seaborn as sns
import copy
import pickle
from scipy.stats.stats import pearsonr
import scipy
import scipy.spatial.distance as ssd

##### Functions for Parsing Data, Extracting Results, Saving, ... #####

## define functions to read in data
def read_data(file_path):
    """ 
    Parameters:
        file_path --> the path of the input file as a string. Input file is presumed to be tab (.txt) or coma (.csv) delimited expression table, with the column names being gene names, and row names being cell names.
    Returns:
        A pandas DataFrame with cell names set as column labels and gene names set as row index.
    """
    if type(file_path) is str:
        if file_path.endswith(".csv"):
            df=pd.read_csv(file_path,index_col=0)  
        elif file_path.endswith(".txt"):
            df=pd.read_table(file_path,index_col=0)
        else:
            df=pd.read_pickle(file_path)
            if type(df) is pd.sparse.frame.SparseDataFrame:
                df=df.to_dense()
        return df
    else:
        print("Invalid input file path(s).")
        return

## define functions to load and save result datasets as pkl objects
def save_obj(obj, name):
    """
    Parameters:
        obj --> the python object to be saved.
        name --> the path of the pickle file the object will be saved to.
    """
    if name.endswith(".pkl"):
        filename=name
    else:
        filename=name+".pkl"
    with open(filename, 'wb') as f:
        try:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        except:
            pickle.dump(obj, f)

def load_obj(name):
    """
    Parameters:
        name --> the path of the pickle file to be loaded.
    Returns:
        The python object saved in the pickle file.
    """
    if name.endswith(".pkl"):
        filename=name
    else:
        filename=name+".pkl"
    re_load=False
    with open(filename, 'rb') as f:
        try:
            obj=pickle.load(f)
        except:
            re_load=True
    if re_load:
        with open(filename, 'rb') as f:
            u = pickle._Unpickler(f)
            u.encoding = 'latin1'
            obj = u.load()
            #obj=pickle.load(f,encoding='latin1')
    return obj
    
## Define a function which, after nmf run, prints out the top n genes that defines each gene module
def print_top_genes(G, n_top_genes = 30,component = None, prt= False, save_tbl=False):
    """
    Parameters:
        G --> The G matrix (number of genes by K) resulted from NMF run, in format of a pandas DataFrame, with gene names as the index.
        n_top_genes --> the number of top ranking genes to display or save.
        component --> the module indices to display or save. If None, all modules will be included.
        prt --> whether to print the results in the environment in addition to returning the resultant table.
        save_tbl --> whether to save a table of the top ranking genes as a csv file. Default is False, not saving. Give the path of the csv file if saving is wanted.
    Returns:
        A pandas DataFrame with the top ranking genes and their weights for each module.
    """
    G_=copy.deepcopy(G)
    if component is None:
        components_ = np.array(G_.T)
    else:
        components_ = np.array(G_.T)[component]
    string = ""
    gene_names=list(G_.index)
    top_gene=pd.DataFrame()
    for topic_idx, topic in enumerate(components_):
        if component is not None:
            topic_idx=component[topic_idx]
        if prt:
            string += "Gene Module #%d:" % topic_idx + '\n'
            string += " ".join([gene_names[i] for i in topic.argsort()[:-n_top_genes - 1:-1]])
            string += "\n"
            string += ", ".join([str(topic[i]) for i in topic.argsort()[:-n_top_genes -1:-1]])
            string += "\n"
        top_gene['Module '+str(topic_idx)]=[gene_names[i] for i in topic.argsort()[:-n_top_genes - 1:-1]]
        top_gene['Weights '+str(topic_idx)]=topic[topic.argsort()[:-n_top_genes - 1:-1]]
    if prt:
        print(string)
    if save_tbl:
        if save_tbl.endswith(".csv"):
            top_gene.to_csv(save_tbl)
        else:
            top_gene.to_csv(save_tbl+'.csv')
    return top_gene

##### Functions for Pre-processing Datasets and Run NMF #####

## Define a function to scale the input expression matrix such that every gene has a maximun expression of 10
def max_scale(pd_df,max_value=None, log_space=True):
    """
    Parameters:
        pd_df --> an expression matrix in form of a pandas DataFrame. Each row corresponds to a gene, each column corresponds to a cell.
        max_value --> all genes will have this maximum expression value across the dataset after scaling. If set to None, the median of max vavlues of all genes will be used.
        log_space --> whether the input expression matrix is in log space. If true, the matrix will be transformed back to linear space, scaled, and then re-transformed to log space.
    Returns:
        A pandas DataFrame of the scaled expression matrix.
    """
    if type(pd_df) is pd.core.frame.DataFrame:
        pd_df_=pd_df.copy()
        df=np.array(pd_df_)
        if max_value is None:
            max_value=np.median(np.max(np.array(pd_df),axis=1))
        if log_space:
            df=np.expm1(df)
            max_value=np.expm1(max_value)
        df_fit=MaxAbsScaler().fit_transform(df.T)
        df_fit=df_fit*max_value
        if log_space:
            df_fit=np.log1p(df_fit)
        pd_df_.iloc[:,:]=df_fit.T
        return pd_df_
    else:
        print("Wrong input data type.")
        return
        
def med_scale(pd_df, med=None,log_space=True, med_only=False):
    """
    Parameters:
        pd_df --> an expression matrix in form of a pandas DataFrame. Each row corresponds to a gene, each column corresponds to a cell.
        med --> all genes will have this non-zero median expression value across the dataset after scaling.
        med_only --> whether to return the non-zero median values for the genes only. If true, only non-zero median values will be returned and the expression matrix will not be scaled.
        log_space --> whether the input expression matrix is in log space. If true, the matrix will be transformed back to linear space, scaled, and then re-transformed to log space.
    Returns:
        A pandas DataFrame of the scaled expression matrix.
    """
    if type(pd_df) is pd.core.frame.DataFrame:
        pd_df_=pd_df.copy()
        df=np.array(pd_df_)
        if log_space:
            df=np.expm1(df)
            if med is not None:
                med=np.expm1(med)
        df_masked=np.ma.masked_where(df == 0, df)
        non0med=np.ma.median(df_masked, axis=1).filled(0)
        if med is None:
            med=np.median(non0med)
        if med_only:
            return(non0med)
        else:
            df_2=df*med / non0med[:,None]
            if log_space:
                df_2=np.log1p(df_2)
            pd_df_.iloc[:,:]=df_2
            return(pd_df_)
    else:
        print("Wrong input data type.")
        return


## Define a function to randomly permute each row of a dataset
def permute(df):
    """
    Parameters:
        df --> an expression matrix in form of either a pandas DataFrame or a numpy array.
    Returns:
        A matrix of the same input format with values in each row randomly permuted.
    """
    df_=df.copy()
    data_arr=np.array(df_)
    sx,sy =data_arr.shape
    perm_data=[row[np.random.permutation(sy)] for row in data_arr]
    if type(df_) is pd.core.frame.DataFrame:
        perm_df=df_
        perm_df.iloc[:,:]=np.array(perm_data)
    else:
        perm_df=np.array(perm_data)
    return perm_df


## define functions to run nmf and visulize results
def run_nmf(data_frame, n_groups, rand_state=None, alpha=.25, l1=.5, max_iter=25000, tol=1e-7, rep=1, sub=None, init=None,verbose=True):
    """
    This function calls the sklearn.decomposition.NMF function from scikit-learn library for each individual NMF run. 
    Parameters:
        data_frame --> the expression matrix to be decomposed, in form of a pandas DataFrame, with each row corresponding to a gene and each column corresponding to a cell.
        n_groups --> number of modules to reduce the data to. See n_components argument in sklearn.decomposition.NMF.
        rand_state --> if int (for rep=1), or list of int (for rep>1), rand_state is the seed used by the random number generator. If None, the random number generator is the RandomState instance used by np.random. See random_state argument in sklearn.decomposition.NMF.
        alpha --> a non-negative constant that multiplies the regularization terms. See alpha argument in sklearn.decomposition.NMF.
        l1 --> a number between 0 and 1 that specifies the regularization mixing between l1 and l2 penalties. See l1_ratio argument in sklearn.decomposition.NMF.
        max_iter --> maximum number of iterations before terminating the nmf run. See max_iter argument in sklearn.decomposition.NMF.
        tol --> tolerance of the stopping condition. See tol argument in sklearn.decomposition.NMF.
        rep --> number of repeated nmf runs to perform. 
        sub --> whether to randomly subsample a fraction of the input data_frame for each nmf run. If int, specified number of cells will be subsampled. If float, specified fraction of cells will be subsample. If None, no sabsampling will be performed.
        init --> method used to initialize the procedure. See init argument in sklearn.decomposition.NMF.
        verbose --> whther to print out progress as the script is running.
    Returns:
        A dictionary with keys "rep0", "rep1", .... Each contain a dictionary storing the NMF result for that repeated run.
        Each "rep" dictionary has keys:
            "G": the genes by modules matrix in form of a pandas DataFrame, retaining the same index as the input data_frame. Corresponding to matrix W in sklearn.decomposition.NMF.
            "C": the modules by cells matrix in from of a pandas DataFrame, retaining the same columns as the input data_frame. Corresponding to matrix H in sklearn.decomposition.NMF.
            "err": the Frobenius norm of the matrix difference, or beta-divergence, between the input data_frame and the reconstructed data GxC from the fitted model. Corresponding to the reconstruction_err_ attribute of the sklearn.decomposition.NMF result.
            "n_iter": actual number of iterations. Corresponding to the n_iter_ attribute of the sklearn.decomposition.NMF result.
    """
    X_=data_frame.copy()
    results={}
    if rand_state is not None: ## one random state for each replicate
        if len(rand_state) != rep:
            print("Error: number of random states doesn't match the number of repeats.")
            return
    else:
        rand_state=np.repeat(None,rep)
    
    if sub is None:
        sub_spl=False
    else:
        num_cell=len(X_.columns)
        if type(sub)== int: ## sub specifies how many cells to randomly pick for each run
            if sub >= num_cell:
                print('Error: invalid sub argument. Number of subsampled cells must be less than the number of cells in the dataset.')
                return
            sub_spl=True
            num_sub=sub
        elif type(sub)== float:
            if num_sub >= 1.0:
                print('Error: invalid sub argument. Fraction of subsampled cells must be less than 1.')
                return
            sub_spl=True
            num_sub=int(num_cell*sub)
        else: 
            print('Error: invalid sub argument type.\nAllowed input: \
                      \n\tOne integer indicating the number of cells to subsample in every replication. \
                      \n\tOne float indicating the fraction of cells to subsample in every replication.')
            return
    for i in range(0,rep):
        if not sub_spl:
            X=X_
        else:
            cell_ind=np.random.permutation(num_cell)
            cell_use=cell_ind[:num_sub]
            X_arr=np.array(X_)
            X=X_arr[:,cell_use]
        if verbose:
            print("Running rep"+str(i)+":")
        nmf = NMF(n_components=n_groups, random_state=rand_state[i], alpha=alpha, l1_ratio=l1,max_iter=max_iter,tol=tol, init=init).fit(X)
        C=nmf.components_ ## n_groups X num_cells
        err=nmf.reconstruction_err_
        n_iter=nmf.n_iter_
        if verbose:
            print('error:', err)
            print('n_iter:', n_iter)
        G=pd.DataFrame(nmf.transform(X),index=X_.index)
        if not sub_spl:
            C_df=pd.DataFrame(C,columns=X_.columns)
        else:
            C_df=pd.DataFrame(C,columns=X_.columns[cell_use])
        results["rep"+str(i)]={"G":G, "C":C_df, "err":err, "n_iter":n_iter}
    return results

##### Functions for Calculating and Ploting Some Statistics from NMF Run Results #####

## Define a function that builds a consensus matrix
def calc_consens(result_k,scl=True,M='C'):
    """
    Parameters:
        result_k --> A dictionary resulted from run_nmf function that contains multiple repeated runs for a single K (n_group).
        scl --> whether to scale each matrix before calculating the consensus matrix such that all components (module) have the same maximum value.
        M --> the matrices used for calculating the consensus matrix. Acceptable values: 'C', or 'G'. If M='C' and there was sub-sampling when running nmf, the function will find the constant cells among the repeated runs and build consensus matrix with only these cells.
    Returns:
        A consensus matrix in form of a 2-D numpy array.
    """
    #result_=copy.deepcopy(result_k)
    result_=result_k
    reps=[i for i in result_.keys() if i.startswith('rep')]
    num_rep=len(reps)
    con_sum=0
    if M=='C': ## test if the cells are constant across C's in different reps
        if len(reps)>1:
            cells1=result_[reps[0]]['C'].columns
            cells2=result_[reps[1]]['C'].columns
            if all(cells1==cells2):
                sub_spl=False
            else:
                sub_spl=True
                cells=set(cells1)
                for i in reps:
                    cellsi=result_[i]['C'].columns
                    cells=cells.intersection(set(cellsi))
                cells=list(cells)
        else:
            sub_spl=False
    for i in reps:
        print('building connectivity matrix for '+i+'...')
        if M=='C':
            if not sub_spl:
                C=np.array(result_[i]['C'])
            else:
                C=np.array(result_[i]['C'][cells])
        elif M=='G':
            G=np.array(result_[i]['G'])
            C=G.T
        else:
            print('Error: invalid input for M. Use either \'C\' or \'G\'.')
        if scl:
            C_fit=MaxAbsScaler().fit_transform(C.T)
            C=C_fit.T
        ## look for the maximum in each column (only want the index)
        group=C.argmax(axis=0) ## an array of #cells integers, each specifies the class the particular cell is in
        C1, C2 = np.meshgrid(group,group)
        con_sum+=np.array(C1==C2,dtype=int)
    con=con_sum/num_rep
    con=con.astype('float32')
    return con ## Returns a consensus matrix over all the replications for K


def consis_plt(con,fg_sz=[10,4],ylim1=None, ylim2=None, save=False,plot=False,logyscale=True,bins=50):
    """
    Parameters:
        con --> the consensus matrix returned by calc_consens.
        fg_sz --> list of length 2. Size of the output figures.
        ylim1 --> if a list of length 2, the two numbers in the list specify the lower and upper bound of y-axis range to plot for the consensus matrix elements plot.
        ylim2 --> if a list of length 2, the two numbers in the list specify the lower and upper bound of y-axis range to plot for the deviation from perfect consistency plot.
        logyscale --> whether to transform y-axis to log scale in plots.
        bins --> number of bins in histograms.
    Save can be either False or a string for the prefix of output .pdf file.
    If save is False, no figures will be saved.
    If save is specified, a pdf file with the specified name will be generated
    """
    if save and not plot:
        plot=True
    if save:
        pp = PdfPages(save+'.pdf')
    con_M=np.array(con)
    if plot: 
        fig_size=np.array(fg_sz)
    con_el=con_M.flatten()
    num_cell=con_M.shape[0]
    ones=np.where(con_el==1)[0]
    zeros=np.where(con_el==0)[0]
    con_el_c=np.delete(con_el,ones[:num_cell]) ## delete the diagonal elements
    con_dev=abs(con_el_c-0.5)
    con_dev=abs(con_dev-0.5)
    mean_dev=np.mean(con_dev)
    if plot:
        plt.rcParams["figure.figsize"] = fig_size
        plt.figure()

        plt.subplot(121)
        plt.hist(con_el,bins=bins)
        plt.title('Distribution of values \nin consensus matrix')
        plt.ylabel('Counts')
        if logyscale:
            plt.yscale('log')
        if ylim1 is not None:
            plt.ylim((ylim1[0],ylim1[1]))

        plt.subplot(122)
        plt.hist(con_dev,bins=bins)
        plt.xlim((0,0.51))
        if logyscale:
            plt.yscale('log')
        if ylim2 is not None:
            plt.ylim((ylim2[0],ylim2[1]))
        plt.xlabel('Deviatioin from 0/1')
        plt.title('Mean Deviation='+str(mean_dev))
        if save:
            pp.savefig()
            pp.close()
    return mean_dev

def calc_cophenet(con,method='average'):
    """ 
    Parameters:
        con --> a consensus matrix (such as output from calc_consens).
    Returns:
        A double. The cophenetic coefficient calculated from the input consensus matrix. 
    """
    dis=1-con
    distArray = ssd.squareform(dis)
    Z=scipy.cluster.hierarchy.linkage(distArray,method=method)
    c=scipy.cluster.hierarchy.cophenet(Z,distArray)
    return c[0]

def stability_cpr(tbl, perm_tbl=None, fg_sz=[7,5], xtick_sz=14, save=False, ttl=False, stats=["inconsistency_G","inconsistency_C","cophenetic_C","cophenetic_G"]):
    """
    Parameters:
        tbl --> a pandas DataFrame with column names as the k values and the rownames corresponding to the stats to be plotted.
        fg_sz --> a list of length 2 to specify the size of the output figure.
        xtick_sz --> size of the xtick labels.
        save --> whether to save the figure as a pdf file. If a path is supplied, the function will save a pdf file to the specified path.
        ttl --> desired title on the figure.
    """
    if save:
        if save.endswith(".pdf"):
            pp = PdfPages(save)
        else:
            pp = PdfPages(save+'.pdf')
    plt.rcParams["figure.figsize"] = fg_sz
    plt.rcParams['xtick.labelsize'] = xtick_sz
    for stat2plot in stats:
        plt.figure()
        plt1, =plt.plot(tbl.columns,tbl.loc[stat2plot],label='Expression Data',color='firebrick',linewidth=2,linestyle='--')
        plt.plot(tbl.columns,tbl.loc[stat2plot],'o',mec="firebrick",color='firebrick')
        if perm_tbl is not None:
            plt2, =plt.plot(perm_tbl.columns,perm_tbl.loc[stat2plot],label='Randomly Permuted Data',color='silver',linewidth=2,linestyle='--')
            plt.plot(perm_tbl.columns,perm_tbl.loc[stat2plot],'o',mec="silver",color='silver')
        if ttl:
            plt.title(ttl,fontsize=16)
        else:
            plt.title(stat2plot,fontsize=16)
        plt.yticks(fontsize=8)
        plt.xlabel('K',fontsize=18)
        if stat2plot.startswith("inconsistency"):
            plt.ylabel('Mean Deviation from \nPerfect Consistency (0/1)',fontsize=14)
        elif stat2plot.startswith("cophenetic"):
            plt.ylabel('Cophenetic Coefficient',fontsize=14)
        if perm_tbl is not None:
            plt.legend(handles=[plt1,plt2],frameon=False)
        if save:
            pp.savefig()
    if save:
        pp.close()
    return

class AttributeDict(dict):
    #extend dictionary object to set atributes as for Matlab objects.
    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            raise AttributeError
    def __setattr__(self, attr, value):
        self[attr] = value

class nmf_reps:
    def __init__(self):
        """
        Documentation: This is the master class for computing nmf over many different K values
        """
        pass

    def set_param(self, scale="median", Ks=range(10,20), rep=3, rand_state=None, alpha=.25, l1=.5, max_iter=25000, tol=1e-7, sub=None, init=None,verbose=True, permute=False):
        """
        Sets paramaters.
        Parameters:
        scale --> whether to scale the expression matrix. Can be False, "median" or "max".
        Ks --> the range of K (ranks) to run NMF for. Can be list, np.array, or range of integers.
        rep --> an integer to specify how many repeated runs of NMF to perform for each K value specified in Ks.
        permute --> True or False. Whether to automatically generate a permuted dataset when set_data is called.
        """
        Params = AttributeDict()
        Params.rand_state = rand_state
        Params.Ks=Ks
        Params.alpha = alpha
        Params.l1 = l1
        Params.max_iter = max_iter 
        Params.tol = tol
        Params.rep = rep
        Params.sub = sub
        Params.init = init
        Params.verbose = verbose
        Params.scale=scale
        Params.permute = permute
        self.Params=Params
        return

    def set_data(self,data,scaled=None,permuted=None,verbose=True):
        """
        Parameters:
        data --> expression matrix in format of a pandas DataFrame. rows = genes. columns = cells. 
        scaled --> the scaled (with max or median) dataset used for NMF.
        permuted --> the permuted dataset run as a control.
        """
        self.data={}
        ## test to see if data contains all zero rows or columns
        if len(np.where(np.sum(data,axis=1)==0)[0])>0 or len(np.where(np.sum(data,axis=0)==0)[0])>0:
            print("Removing "+str(len(np.where(np.sum(data,axis=1)==0)[0]))+" all-zero rows and "+str(len(np.where(np.sum(data,axis=0)==0)[0])) +" all-zero columns...")
            data=data.iloc[np.where(np.sum(data,axis=1)>0)[0],np.where(np.sum(data,axis=0)>0)[0]]
        self.data["raw"] = data
        Params=self.Params
        if scaled is None:
            if Params.scale=="median":
                if verbose:
                    print("Scaling expression data by non-zero median...")
                self.data["scaled"]=med_scale(data)
            elif Params.scale=="max":
                if verbose:
                    print("Scaling expression data by max...")
                self.data["scaled"]=max_scale(data)
            else:
                self.data["scaled"]=None
        else:
            self.data["scaled"] = scaled
        if permuted is None:
            if Params.permute:
                if self.data["scaled"] is None:
                    if verbose:
                        print("Permuting expression data...")
                    self.data["permuted"]=permute(self.data["raw"])
                else:
                    if verbose:
                        print("Permuting scaled expression data...")
                    self.data["permuted"]=permute(self.data["scaled"])
            else:
                self.data["permuted"]=None
        else:
            self.data["permuted"]=permuted
        return


    def nmf_results(self, permuted=False,prt_top_genes=False):
        """
        Run NMF based on the parameters set in self.Params and store the results in self.results.
        Parameters:
            print_top_genes --> whether to print the top 30 genes and their weights in each module. Can be True, False, or a path to save the table as a .csv file. 
        """
        Params=self.Params
        Ks=Params.Ks
        self.results={}
        run_permute=False
        if permuted == "only":
            run_permute=True
            data_perm=self.data["permuted"]
            if data_perm is None:
                print("Error, no permuted dataset is in store.")
                return
        else:
            if self.data["scaled"] is None:
                data_use=self.data["raw"]
            else:
                data_use=self.data["scaled"]
            if permuted:
                run_permute=True
                data_perm=self.data["permuted"]
                if data_perm is None:
                    print("Error, no permuted dataset is in store.")
                    return
        if run_permute:
            self.permuted_results={}
        for k in Ks:
            if Params.verbose:
                print("Running NMF for K=" +str(k)+"...")
            if permuted is not "only":
                self.results["K="+str(k)]=run_nmf(data_use, k, rand_state=Params.rand_state, alpha=Params.alpha, l1=Params.l1, max_iter=Params.max_iter, tol=Params.tol, rep=Params.rep, sub=Params.sub, init=Params.init,verbose=Params.verbose)
                if prt_top_genes:
                    if type(prt_top_genes) is str:
                        print_top_genes(self.results["K="+str(k)]["rep0"]["G"],prt=False,save_tbl=prt_top_genes+"_K"+str(k)+"_rep0")
                        if Params.rep>1:
                            print_top_genes(self.results["K="+str(k)]["rep1"]["G"],prt=False,save_tbl=prt_top_genes+"_K"+str(k)+"_rep1")
                    else:
                        print_top_genes(self.results["K="+str(k)]["rep0"]["G"],prt=True,save_tbl=False)
                        if Params.rep>1:
                            print_top_genes(self.results["K="+str(k)]["rep1"]["G"],prt=True,save_tbl=False)
            if run_permute:
                if Params.verbose:
                    print("  running permuted dataset...")
                self.permuted_results["K="+str(k)]=run_nmf(data_perm, k, rand_state=Params.rand_state, alpha=Params.alpha, l1=Params.l1, max_iter=Params.max_iter, tol=Params.tol, rep=Params.rep, sub=Params.sub, init=Params.init,verbose=Params.verbose)
                if prt_top_genes:
                    if type(prt_top_genes) is str:
                        print_top_genes(self.permuted_results["K="+str(k)]["rep0"]["G"],prt=False,save_tbl=prt_top_genes+"_permuted_K"+str(k)+"_rep0")
                    else:
                        print_top_genes(self.permuted_results["K="+str(k)]["rep0"]["G"],prt=True,save_tbl=False)
        return

    def calc_stability(self,stats=["inconsistency_G","inconsistency_C","cophenetic_C","cophenetic_G"],permuted=False):
        """
        Parameters:
            stats --> a list of strings specifying the statistics to calculate.
            permuted --> whether to calculate the statistics for the nmf run with the permuted dataset. If False, only stats on real dataset will be calculated. If True, stats for both real and permuted datasets will be calculated. If "only", only stats for permuted statsets will be calculated.
        Returns:
            Adds .stability and/or .permuted_stability attribute(s) to self.
        """
        Params=self.Params
        Ks=Params.Ks
        run_permute=False
        if permuted == "only":
            run_permute=True
            try:
                permuted_results=self.permuted_results
                stats_tbl_perm=pd.DataFrame(columns=Ks,dtype="float",index=stats)
            except AttributeError:
                print("Error, no result for permuted dataset is in store.")
                return
        else:
            stats_tbl=pd.DataFrame(columns=Ks,dtype="float",index=stats)
            if permuted:
                run_permute=True
                try:
                    permuted_results=self.permuted_results
                    stats_tbl_perm=pd.DataFrame(columns=Ks,dtype="float",index=stats)
                except AttributeError:
                    print("Error, no result for permuted dataset is in store.")
                    return
        if permuted is not "only":
            if Params.verbose:
                print("Calculating stability for results...")
            self.stability=stability_tbl(self.results,Ks,stats=stats,verbose=Params.verbose)
        if run_permute:
            if Params.verbose:
                print("Calculating stability for results from permuted dataset...")
            self.permuted_stability=stability_tbl(permuted_results,Ks,stats=stats,verbose=Params.verbose)
        return

    def plot_stability(self,permuted=None,stats=["inconsistency_G","inconsistency_C","cophenetic_C","cophenetic_G"],save=False,ttl=False,fg_sz=[7,5]):
        """
        Parameters:
            permuted --> whether to plot the stats for results from permuted dataset in the same plot as results from the real dataset. Can be False, True, or None. If None, the function will plot stat for permuted dataset if results for permuted dataset is detected.
        """
        if permuted is None:
            try:
                perm_tbl=self.permuted_stability
                permuted=True
            except AttributeError:
                permuted=False
        if permuted:
            stability_cpr(tbl=self.stability, perm_tbl=perm_tbl, stats=stats, save=save, ttl=ttl, fg_sz=fg_sz)
        else:
            stability_cpr(tbl=self.stability, stats=stats, save=save, ttl=ttl, fg_sz=fg_sz)
        return

    def extract_err(self,permuted=False,Ks=None,extr='err',measure="mean"):
        """
        Parameters:
            permuted --> whether to extract 'err' or 'n_iter' from results from the permuted dataset instead of results from the real dataset.
            Ks --> the k range to integrate into a table. If None, the Ks values in self.Params will be used.
            extr --> whether to extract 'err' or 'n_iter' from the results.
            measure --> whether to plot mean or median of each K over multiple replicates. Can be 'mean' or 'median'.
        Returns:
            A pandas DataFrame with the K values as the columns names and 'err' or 'n_iter' for each repeated run as a row. 
            'mean' or 'median' will be added to the table as a row depending on the measure argument.
        """
        if Ks is None:
            Params=self.Params
            Ks=Params.Ks
        err_tbl=pd.DataFrame(columns=sorted(Ks),dtype="float")
        if permuted:
            results_=self.permuted_results
        else:
            results_=self.results
        for k in Ks:
            results_k=results_["K="+str(k)]
            reps=[i for i in results_k.keys() if i.startswith("rep")]
            for rep in reps:
                err_tbl.loc[rep,k]=results_k[rep][extr]
        if measure=="mean":
            err_tbl.loc[measure]=np.nanmean(err_tbl,axis=0)
        elif measure=="median":
            err_tbl.loc[measure]=np.nanmedian(err_tbl,axis=0)
        return err_tbl


    def plot_err(self,stat="err",Ks=None,permuted=None,measure="mean",save=False,fg_sz = [7,5], xtick_sz=14, ttl = False):
        """
        Parameters:
            stat --> what to plot. Either "err" or "n_iter".
            perm --> whether to plot the stat for results from permuted dataset in the same plot as results from the real dataset. Can be False, True, or None. If None, the function will plot stat for permuted dataset if results for permuted dataset is detected.
            save --> whether to save the plot as a pdf. Can be False or the path of the .pdf file to be saved.
        """
        ## Get tables for "err" or "n_iter"
        if Ks is None:
            Params=self.Params
            Ks=Params.Ks
        if permuted is None:
            try:
                perm_results=self.permuted_results
                permuted=True
            except AttributeError:
                permuted=False
        stat_tbl=self.extract_err(permuted=False,Ks=Ks,extr=stat,measure=measure)
        if permuted:
            perm_tbl=self.extract_err(permuted=True,Ks=Ks,extr=stat,measure=measure)
        
        ## Generate plot
        plt.rcParams["figure.figsize"] = fg_sz
        plt.rcParams['xtick.labelsize'] = xtick_sz
        if save:
            if save.endswith(".pdf"):
                pp = PdfPages(save)
            else:
                pp = PdfPages(save+'.pdf')
        all_rep=[rep for rep in stat_tbl.index if rep.startswith("rep")]
        num_rep=len(all_rep)
        Xs=np.tile(Ks,reps=num_rep)
        Ys=np.array(stat_tbl.loc[all_rep]).flatten()
        plt.figure()
        plt1, =plt.plot(Ks,np.array(stat_tbl.loc[measure]),label='Expression Data',color='firebrick',linewidth=2,linestyle='--')
        plt.plot(Xs,Ys,'o',mec="firebrick",color='firebrick')
        if permuted:
            all_rep=[rep for rep in perm_tbl.index if rep.startswith("rep")]
            num_rep=len(all_rep)
            Xs=np.tile(Ks,reps=num_rep)
            Ys=np.array(perm_tbl.loc[all_rep]).flatten()
            plt2, =plt.plot(Ks,np.array(perm_tbl.loc[measure]),label='Randomly Permuted Data',color='silver',linewidth=2,linestyle='--')
            plt.plot(Xs,Ys,'o',mec="silver",color='silver')
        if ttl:
            plt.title(ttl,fontsize=16)
        plt.xlabel('K',fontsize=18)
        if stat == 'err':
            plt.ylabel('Reconstruction Error',fontsize=16)
        elif stat == 'n_iter':
            plt.ylabel('Actual # Iteration',fontsize=16)
        if permuted:
            plt.legend(handles=[plt1,plt2],frameon=False)
        if save:
            pp.savefig()
            pp.close()
        return

    def rebuild_hex_class(self, k, rep, data_use="scaled", fg_sz=10, min_val=0.1, save=False,title=None):
        """
        Parameters:
            data_use --> "raw", "scaled", or "permuted"
        """
        ## can take both np.array and pd.df as input
        G=self.results["K="+str(k)][rep]['G']
        C=self.results["K="+str(k)][rep]['C']
        re_const=np.dot(np.array(G),np.array(C))
        real_data=np.array(self.data[data_use])
        re_const_arr=re_const.flatten('C')
        real_data_arr=real_data.flatten('C')
        ind_use=np.intersect1d(np.where(re_const_arr>=min_val)[0],np.where(real_data_arr>=min_val)[0])
        hexplot = sns.jointplot(real_data_arr[ind_use],re_const_arr[ind_use],kind="hex", size=fg_sz)
        hexplot.set_axis_labels(xlabel='Real Data', ylabel='Reconstructed Data')
        if title is not None:
            plt.title(title)
        cax = hexplot.fig.add_axes([1, 0.25, 0.04, 0.5])# size and placement of bar [left, bottom, width, height]
        plt.colorbar(cax=cax)
        if save:
            pp = PdfPages(save+'.pdf')
            pp.savefig()
            pp.close()
        return

    def top_ratio(self,permuted=False, Ks=None,reps_use=None):
        """
        Parameters:
        """
        if Ks is None:
            Ks=self.Params.Ks
        if permuted:
            results_=copy.deepcopy(self.permuted_results)
        else:
            results_=copy.deepcopy(self.results)
        ratio_tbl=pd.DataFrame(columns=["K","ratio(rank1/rank2)","rep"])
        for k in Ks:
            if reps_use is None:
                reps=[i for i in results_["K="+str(k)].keys() if i.startswith('rep')]
            else:
                reps=["rep"+str(i) for i in reps_use]
            ratio12=np.array([])
            rep_val=[]
            for rep in reps:
                result_G=results_["K="+str(k)][rep]["G"]
                sorted_G=np.sort(result_G,axis=0)
                ratio_rep=sorted_G[-1,:]/sorted_G[-2,:]
                ratio12=np.append(ratio12,ratio_rep)
                rep_val=rep_val+[rep]*len(ratio_rep)
            K_val=[k]*len(ratio12)
            append_tbl=np.array([K_val,ratio12,rep_val])
            ratio_tbl=ratio_tbl.append(pd.DataFrame(append_tbl.T,columns=["K","ratio(rank1/rank2)","rep"]))
        ratio_tbl["K"]=pd.to_numeric(ratio_tbl["K"])
        ratio_tbl["ratio(rank1/rank2)"]=pd.to_numeric(ratio_tbl["ratio(rank1/rank2)"],errors='coerce')
        return(ratio_tbl)

    def mean_GC(self,K,reps_use=None,st=0):
        """
        Parameters:
            K --> the K value whose results are to be used in this function. self.results must have "K=K" as a key.
            reps_use --> which repeated runs to use. Can be either a list of integers indicating the rep indices, or None, in which case results from all repeated runs will be used.
            st --> results from which repeated run should be used as the reference matrices for other results to be matched to. If st=0, the G and C matrices from the first repeated run (rep0) will be used as reference.
        """
        results_=copy.deepcopy(self.results["K="+str(K)])
        rep_st="rep"+str(st)
        G0=results_[rep_st]['G']
        C0=results_[rep_st]['C']
        if reps_use is None:
            reps=[i for i in results_.keys() if i.startswith('rep')]
        else:
            reps=["rep"+str(i) for i in reps_use]
        Gf=np.array(G0)
        Cf=np.array(C0)
        tot=len(reps)
        for rep in reps:  ## calculate correlation matrices for each rep with the rep0 results
            if rep != rep_st:
                result_rep=results_[rep]
                G=np.array(result_rep['G'])
                C=np.array(result_rep['C'])
                num_k=G.shape[1]
                G_corr=np.zeros((num_k,num_k))
                C_corr=np.zeros((num_k,num_k))
                for j in range(0,num_k):
                    for k in range(0,num_k): ## calculate correlations of one group in rep with all groups in the rep0
                        corrG=pearsonr(G[:,j],np.array(G0)[:,k])
                        G_corr[j,k]=corrG[0]
                        corrC=pearsonr(C[j,:],np.array(C0)[k,:])
                        C_corr[j,k]=corrC[0]
                ## match indices (which g in rep matches which g in rep0). Find the index of max in each column
                Ggroup=G_corr.argmax(axis=0)
                Cgroup=C_corr.argmax(axis=0)
                if all(Ggroup==Cgroup):
                    if len(np.unique(Ggroup))==len(Ggroup): ## if one-to-one match
                        ## reorder the rep(i) G and C
                        G_new=G[:,Ggroup]
                        C_new=C[Ggroup,:]
                        Gf=Gf+G_new
                        Cf=Cf+C_new
                    else: ## if not one-to-one match
                        tot=tot-1
                        print('Skipping '+rep)
                        print('Multiple modules matched to a single one. Group assignments:')
                        print(Ggroup)
                else:
                    tot=tot-1
                    print('Skipping '+rep)
                    print('Group mapping by G and C do not agree.')
                    print('G group:')
                    print(Ggroup)
                    print('C group:')
                    print(Cgroup)
        ## divide the summation matrics by number of replicates
        Gf=Gf/tot
        Cf=Cf/tot
        G0.iloc[:,:]=Gf
        C0.iloc[:,:]=Cf
        mean={'G':G0,'C':C0}
        return mean

    def params2txt(self,save):
        """
        Parameters:
            save --> path for the output .txt file. If set to False, the parameters will be printed in the current environment.
        """
        if save:
            if save.endswith(".txt"):
                save_path=save
            else:
                save_path=save+".txt"
            with open(save_path, 'w') as param_txt:
                for param in self.Params.keys():
                    param_txt.write(param+": "+str(self.Params[param])+"\n")
        else:
            for param in self.Params.keys():
                print(param+": "+str(self.Params[param]))
        return

def stability_tbl(results,Ks=None,stats=["inconsistency_G","inconsistency_C","cophenetic_C","cophenetic_G"],verbose=True):
    if Ks is None:
        Ks_key=[i for i in results.keys() if i.startswith('K=')]
        Ks=[int(i.split("=")[1]) for i in Ks_key]
        Ks=sorted(Ks)
    stats_tbl=pd.DataFrame(columns=Ks,dtype="float",index=stats)
    if "inconsistency_G" in stats or "cophenetic_G" in stats:
        consis=[]
        cophe=[]
        for k in Ks:
            if verbose:
                print("Calculating stats on matrices G for K="+str(k)+"...")
            consens=calc_consens(results["K="+str(k)],M='G')
            if "inconsistency_G" in stats:
                consis.append(consis_plt(consens,save=False,plot=False))
            if "cophenetic_G" in stats:
                try:
                    cophe.append(calc_cophenet(consens))
                except ValueError:
                    print("Error: ValueError in scipy.cluster.hierarchy.cophenet(Z,distArray). This bug is fixed in Scipy version 0.19. Converting problematic value to nan...")
                    cophe.append(np.nan)
        if "inconsistency_G" in stats:
            stats_tbl.loc["inconsistency_G"]=np.array(consis)
        if "cophenetic_G" in stats:
            stats_tbl.loc["cophenetic_G"]=np.array(cophe)
    if "inconsistency_C" in stats or "cophenetic_C" in stats:
        consis=[]
        cophe=[]
        for k in Ks:
            if verbose:
                print("Calculating stats on matrices C for K="+str(k)+"...")
            consens=calc_consens(results["K="+str(k)],M='C')
            if "inconsistency_C" in stats:
                consis.append(consis_plt(consens,save=False,plot=False))
            if "cophenetic_C" in stats:
                cophe.append(calc_cophenet(consens))
        if "inconsistency_C" in stats:
            stats_tbl.loc["inconsistency_C"]=np.array(consis)
        if "cophenetic_C" in stats:
            stats_tbl.loc["cophenetic_C"]=np.array(cophe)
    return stats_tbl

def rebuild_hex(G, C, data_use, fg_sz=10, min_val=0.1, save=False,title=None):
    """
    Parameters:
        data_use --> "raw", "scaled", or "permuted"
    """
    ## can take both np.array and pd.df as input
    re_const=np.dot(np.array(G),np.array(C))
    real_data=np.array(data_use)
    re_const_arr=re_const.flatten('C')
    real_data_arr=real_data.flatten('C')
    ind_use=np.intersect1d(np.where(re_const_arr>=min_val)[0],np.where(real_data_arr>=min_val)[0])
    hexplot = sns.jointplot(real_data_arr[ind_use],re_const_arr[ind_use],kind="hex", size=fg_sz)
    hexplot.set_axis_labels(xlabel='Real Data', ylabel='Reconstructed Data')
    if title is not None:
        plt.title(title)
    cax = hexplot.fig.add_axes([1, 0.25, 0.04, 0.5])# size and placement of bar [left, bottom, width, height]
    plt.colorbar(cax=cax)
    if save:
        pp = PdfPages(save+'.pdf')
        pp.savefig()
        pp.close()
    return


    
