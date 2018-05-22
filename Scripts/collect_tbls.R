args = commandArgs(trailingOnly=TRUE)
if (length(args) > 1) stop("Usage: Rscript collect_tbls.R [DIR_TO_READ_AND_SAVE]")
io_dir <- as.character(args[1])
Ks_dir=list.dirs(io_dir,full.names = TRUE, recursive = FALSE)
Ks_all=unlist(lapply(Ks_dir, function(x) unlist(strsplit(x,"/"))[length(unlist(strsplit(x,"/")))]))
Ks=Ks_all[grep("K=",Ks_all)]
result_obj=list()
if("scaled_data.csv" %in% Ks_all){
  result_obj[["scaled_data"]]=read.csv(paste0(Ks_dir,"/","scaled_data.csv"),row.names = 1)
}
if("permuted_data.csv" %in% Ks_all){
  result_obj[["permuted_data"]]=read.csv(paste0(Ks_dir,"/","permuted_data.csv"),row.names = 1)
}
for(K in Ks){
  print(paste0("Collecting results for ",K,"..."))
  result_obj[[K]]=list()
  reps_dir=list.dirs(paste0(io_dir,"/",K,"/"), full.names=TRUE, recursive = FALSE)
  for(rep_dir in reps_dir){
    rep=unlist(strsplit(rep_dir,"/"))[length(unlist(strsplit(rep_dir,"/")))]
    if(length(grep("rep",rep))>0){
      result_obj[[K]][[rep]]=list()
      files=list.files(rep_dir,full.names=FALSE,recursive = FALSE)
      for(file in files){
        if(length(grep(".csv",file))>0){
          result_obj[[K]][[rep]][[unlist(strsplit(file,"[.]"))[1]]]=read.csv(paste0(rep_dir,"/",file),row.names = 1)
        }
      }
    }
  }
}
save(result_obj,file = paste0(io_dir,"/result_tbls.Robj"))