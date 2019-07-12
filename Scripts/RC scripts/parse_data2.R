library('Matrix')
data_dir="/n/boslfs/LABS/schier_lab/everyone/10X/bushra/Yiqun/"
save_dir="/n/boslfs/LABS/schier_lab/users/yiqunwang/NMF/NMF_for_Bushra/Data/"
stages=c("zf6s","zf10s","zf14s","zf18s","zf20hpf","zf24hpf","zf2dpf")
#stages=c("zf2dpf")

save_permute=T
var_genes=list()
for(stage in stages){
	vargen_name=paste0(stage,"_filtVarGenes.txt")
	var_genes[[stage]]=scan(paste0(data_dir,vargen_name),what = "character",sep="\n")
}

var_genes_use=list()
var_genes_use[["zf6s"]]=unique(c(var_genes[["zf6s"]],var_genes[["zf10s"]],var_genes[["zf14s"]]))
var_genes_use[["zf10s"]]=unique(c(var_genes[["zf6s"]],var_genes[["zf10s"]],var_genes[["zf14s"]],var_genes[["zf18s"]]))
var_genes_use[["zf14s"]]=unique(c(var_genes[["zf6s"]],var_genes[["zf10s"]],var_genes[["zf14s"]],var_genes[["zf18s"]],var_genes[["zf20hpf"]]))
var_genes_use[["zf18s"]]=unique(c(var_genes[["zf10s"]],var_genes[["zf14s"]],var_genes[["zf18s"]],var_genes[["zf20hpf"]],var_genes[["zf24hpf"]]))
var_genes_use[["zf20hpf"]]=unique(c(var_genes[["zf14s"]],var_genes[["zf18s"]],var_genes[["zf20hpf"]],var_genes[["zf24hpf"]],var_genes[["zf2dpf"]]))
var_genes_use[["zf24hpf"]]=unique(c(var_genes[["zf18s"]],var_genes[["zf20hpf"]],var_genes[["zf24hpf"]],var_genes[["zf2dpf"]]))
var_genes_use[["zf2dpf"]]=unique(c(var_genes[["zf20hpf"]],var_genes[["zf24hpf"]],var_genes[["zf2dpf"]]))

#var_genes=c()
#for(stage in stages){
#	vargen_name=paste0(stage,"_filtVarGenes.txt")
#	var_genes=c(var_genes,scan(paste0(data_dir,vargen_name),what = "character",sep="\n"))
#}
#var_genes=unique(var_genes)
#print(paste0(length(var_genes)," variable genes in total."))
for(stage in stages){
  print(paste0("parsing data for ",stage,"..."))
  rds_name=paste0(stage,"_filt_matrix.rds")
  stage_data=readRDS(paste0(data_dir,rds_name))
  #vargen_name=paste0(stage,"_filtVarGenes.txt")
  #var_genes=scan(paste0(data_dir,vargen_name),what = "character",sep="\n")
  var_genes=intersect(var_genes_use[[stage]],rownames(stage_data))
  print(paste0(length(var_genes_use[[stage]])," var genes."))
  print(paste0(length(var_genes)," genes used."))
  write.csv(as.matrix(stage_data)[var_genes,], paste0(save_dir,stage,"_vargene.csv"),quote=F)
  if(save_permute){
    perm_vardata=apply(stage_data[var_genes,],1,function(x) x[sample(ncol(stage_data))])
    perm_vardata=t(perm_vardata)
    #perm_vardata=t(perm_vardata)    colnames(perm_vardata)=colnames(stage_data)
    write.csv(perm_vardata, paste0(save_dir,stage,"_vargene_perm.csv"))
  }
}
