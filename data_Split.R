
source("setting_files.R")
source("setting_paramts_data.R")

# usible: year.num  split.ts.prop
# need setting: disease.pred.arr  outcome.pred.arr

if(TYPE=="REAL"){
  #------- real --------#
  data_dir = real_data_dir
  split_dir = real_split_dir
  
} else if(TYPE=="SYN"){
  #------ synthetic -------#
  data_dir = syn_data_dir
  split_dir = syn_split_dir
}

if(!dir.exists(split_dir)) dir.create(file.path(split_dir))

#---------------------------------------#

split.separate <- function(in.task, in.y, in.step){
  if(!in.task %in% c('class','regr')){
    cat("    Wrong task:",in.task,"\n")
    return()
  }
  
  if(in.task=='class') {
    delta.type.1="Adel"
    delta.type.2="AdelN"
  } else {
    delta.type.1="Sdel"
    delta.type.2="SdelN"
  }
  data.type="Dt1"

  ix_name.1=paste(data.type,'_',delta.type.1,'_',in.task,'_',in.y,sep='');
  ix_name.2=paste(data.type,'_',delta.type.2,'_',in.task,'_',in.y,sep='');
  split.fname.1 = paste(split_dir,"/",ix_name.1,'_',in.step,'.ts.csv',sep="")
  split.fname.2 = paste(split_dir,"/",ix_name.2,'_',in.step,'.ts.csv',sep="")
  
  if(file.exists(split.fname.1) && file.exists(split.fname.2)) {
    cat("    split() Skip because existed:\n")
    cat("     ",split.fname.1,"\n")
    cat("     ",split.fname.2,"\n")
    return()
  }
  
  cat("   ",split.fname.1," doing \n")
  
  assign("task.type", in.task, envir = .GlobalEnv)
  assign("task.outcome", in.y, envir = .GlobalEnv)
  data.fname = get.fname("data");
  
  if(file.info(data.fname)$size <= 10){
    ts.idxs=c()
  } else {
    data = read.csv(data.fname, header = T)
    seri.res = find.serial.cases(in.step)
    X.idx.list = seri.res$X.idx.list
    Y.vec = seri.res$Y.vec
    
    if(length(Y.vec)==0){
      ts.idxs = c()
    } else {
      pid.vec = data[X.idx.list[[1]], "pid"] # the 1st step
      split.data = as.data.frame(cbind(pid.vec, Y.vec))
      colnames(split.data) = c("pid","y")
      
      avg.per.pid = mean(table(split.data[,"pid"]))
      
      ts.pnum = floor(nrow(split.data)*split.ts.prop/avg.per.pid)
      ts.pidxs = sort(sample(unique(split.data[,"pid"]), ts.pnum, replace = F))
      ts.idxs = which(split.data[,"pid"] %in% ts.pidxs)
    }
  }
  
  write.table(ts.idxs, split.fname.1, sep=",",  col.names=FALSE, row.names=FALSE)
  write.table(ts.idxs, split.fname.2, sep=",",  col.names=FALSE, row.names=FALSE)
}

#-----------------------------------------------------------------------------------------

n_len = year.num - forcast.steps
# max number of records per patient can produce

if(TYPE=="REAL"){
  #------- real --------#
  year.num = t.num # standardize the name 
} else if(TYPE=="SYN"){
  #------ synthetic -------#
  # NULL
}

# classification
for(oc in disease.pred.arr){
  for (st in steps.arr) {
    split.separate('class',in.y=oc, in.step=st)
  }
}

if(do.regr.tag){
  # regression
  for(oc in outcome.pred.arr){
    for (st in steps.arr) {
      split.separate('regr',in.y=oc, in.step=st)
    }
  }
}


