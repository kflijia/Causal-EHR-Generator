source("setting_files.R")
source("setting_paramts_data.R")

if(TYPE=="REAL"){
  #------- real --------#
  data_dir = real_data_dir
  split_dir = real_split_dir
  
} else if(TYPE=="SYN"){
  #------ synthetic -------#
  data_dir = syn_data_dir
  split_dir = syn_split_dir
}

gen.subset.data.class <- function(in.oc){
  whole.data.fname = paste0(whole.wd,"syn_datasets/Dt1_Adel_class_",in.oc,".data.csv")
  if(!file.exists(whole.data.fname) || file.info(whole.data.fname)$size<=10) {
    cat(whole.data.fname," No Data.\n")
    return(FALSE)
  }
  
  whole.data = read.csv(whole.data.fname, header = T)
  
  whole.Ndata.fname = paste0(whole.wd,"syn_datasets/Dt1_AdelN_class_",in.oc,".data.csv")
  whole.Ndata = read.csv(whole.Ndata.fname, header = T)
  
  whole.case.num = nrow(whole.data)
  whole.pat.num = whole.case.num/n_len
  
  # for class only
  assign("task.outcome", in.oc, envir = .GlobalEnv)
  serial.pids.list = vector("list",length(steps.arr))
  ts.idxs.list = vector("list",length(steps.arr))
  
  for (st in steps.arr) {
    st.seri.res = find.serial.cases(st)
    X.idx.list = st.seri.res$X.idx.list
    Y.vec = st.seri.res$Y.vec
    
    if(length(Y.vec)==0){
      ts.idxs.list[[st]] = c()
    } else {
      pid.vec = whole.data[X.idx.list[[1]], "pid"] # the 1st step
      whole.split.data = as.data.frame(cbind(pid.vec, Y.vec))
      colnames(whole.split.data) = c("pid","y")
      serial.pids.list[[match(st,steps.arr)]] = whole.split.data

      #------ split ------
      whole.split.fname = paste0(whole.wd,"syn_splits/Dt1_Adel_class_",in.oc,"_",st,".ts.csv")
      if(!file.exists(whole.split.fname) || file.info(whole.split.fname)$size<=10){
        ts.idxs.list[[st]] = c()
        next;
      }
      whole.ts.idxs = read.csv(whole.split.fname, header = F)
      whole.ts.idxs = whole.ts.idxs[,1]
      ts.idxs.list[[st]] = whole.ts.idxs
      }
    
    } # for st
  
  
  for (s in sample.size.arr[1:(sample.num-1)]) {
    sub.sample.dir = paste("syn_sample_",s,"/", sep="")
    sub.wd = paste(syn_root_dir,sub.sample.dir,sep="")
    if (!dir.exists(sub.wd)) dir.create(file.path(sub.wd))
    sub.wd.dt = paste0(sub.wd,"syn_datasets/")
    if (!dir.exists(sub.wd.dt)) dir.create(file.path(sub.wd.dt))
    sub.wd.sp = paste0(sub.wd,"syn_splits/")
    if (!dir.exists(sub.wd.sp)) dir.create(file.path(sub.wd.sp))
    
    sub.pat.num = ceiling((s/whole.size)*whole.pat.num)
    sub.case.num = sub.pat.num*n_len

    sub.data.fname = paste0(sub.wd.dt,"Dt1_Adel_class_",in.oc,".data.csv")
    sub.Ndata.fname = paste0(sub.wd.dt,"Dt1_AdelN_class_",in.oc,".data.csv")
    if(file.exists(sub.data.fname) && file.exists(sub.Ndata.fname)){
      cat("    subset.data() Skip because existed:\n")
      cat("     ",sub.data.fname,"\n")
      cat("     ",sub.Ndata.fname,"\n")
      sub.data = read.csv(sub.data.fname, header = T);
    } else {
      sub.data = whole.data[1:sub.case.num,,drop=F]
      sub.Ndata = whole.Ndata[1:sub.case.num,,drop=F]
      write.table(sub.data, sub.data.fname, sep=",",  col.names=TRUE, row.names=FALSE)
      write.table(sub.Ndata, sub.Ndata.fname, sep=",",  col.names=TRUE, row.names=FALSE)
    }
    
    max.pat.id = max(sub.data[,1])
    
    for (st in steps.arr) {
      split.fname.1 = paste0(sub.wd,"syn_splits/Dt1_Adel_class_",in.oc,"_",st,".ts.csv")
      split.fname.2 = paste0(sub.wd,"syn_splits/Dt1_AdelN_class_",in.oc,"_",st,".ts.csv")
      if(file.exists(split.fname.1) && file.exists(split.fname.2)){
        cat("    subset.split() Skip because existed:\n")
        cat("     ",split.fname.1,"\n")
        cat("     ",split.fname.2,"\n")
        next;
      }
      
      whole.ts.idxs = ts.idxs.list[[st]]
      
      if(length(whole.ts.idxs)==0) {
        sub.ts.idxs = c()
      } else {
        step.pids = serial.pids.list[[match(st,steps.arr)]][,"pid"]
        sub.step.idxs = which(step.pids<=max.pat.id) # sub pats = 1:max.pat.id
        sub.ts.idxs = whole.ts.idxs[whole.ts.idxs %in% sub.step.idxs]
      }
      write.table(sub.ts.idxs, split.fname.1, sep=",",  col.names=FALSE, row.names=FALSE)
      write.table(sub.ts.idxs, split.fname.2, sep=",",  col.names=FALSE, row.names=FALSE)
    }
  }
  return(TRUE);
}


gen.subset.data.regr <- function(in.oc){
  whole.data.fname = paste0(whole.wd,"syn_datasets/Dt1_Sdel_regr_",in.oc,".data.csv")
  if(!file.exists(whole.data.fname) || file.info(whole.data.fname)$size<=10) {
    cat(whole.data.fname," No Data.\n")
    return(FALSE)
  }
  
  whole.data = read.csv(whole.data.fname, header = T)
  
  whole.Ndata.fname = paste0(whole.wd,"syn_datasets/Dt1_SdelN_regr_",in.oc,".data.csv")
  whole.Ndata = read.csv(whole.Ndata.fname, header = T)
  
  whole.case.num = nrow(whole.data)
  whole.pat.num = whole.case.num/n_len
  
  for (s in sample.size.arr[1:(sample.num-1)]) {
    sub.sample.dir = paste("syn_sample_",s,"/", sep="")
    sub.wd = paste(syn_root_dir,sub.sample.dir,sep="")
    if (!dir.exists(sub.wd)) dir.create(file.path(sub.wd))
    sub.wd.dt = paste0(sub.wd,"syn_datasets/")
    if (!dir.exists(sub.wd.dt)) dir.create(file.path(sub.wd.dt))
    sub.wd.sp = paste0(sub.wd,"syn_splits/")
    if (!dir.exists(sub.wd.sp)) dir.create(file.path(sub.wd.sp))
    
    sub.pat.num = ceiling((s/whole.size)*whole.pat.num)
    sub.case.num = sub.pat.num*n_len
    
    sub.data.fname = paste0(sub.wd.dt,"Dt1_Sdel_regr_",in.oc,".data.csv")
    sub.data = whole.data[1:sub.case.num,,drop=F]
    sub.Ndata.fname = paste0(sub.wd.dt,"Dt1_SdelN_regr_",in.oc,".data.csv")
    sub.Ndata = whole.Ndata[1:sub.case.num,,drop=F]
    write.table(sub.data, sub.data.fname, sep=",",  col.names=TRUE, row.names=FALSE)
    write.table(sub.Ndata, sub.Ndata.fname, sep=",",  col.names=TRUE, row.names=FALSE)
    
    for (st in steps.arr) {
      #------ split ------
      whole.split.fname = paste0(whole.wd,"syn_splits/Dt1_Sdel_regr_",in.oc,"_",st,".ts.csv")
      if(!file.exists(whole.split.fname) || file.info(whole.split.fname)$size<=10){
        sub.ts.idxs = c()
      } else {
        whole.ts.idxs = read.csv(whole.split.fname, header = F)
        whole.ts.idxs = whole.ts.idxs[,1]
        
        sub.split.case.num = sub.pat.num*(n_len-st+1)
        sub.ts.idxs = whole.ts.idxs[whole.ts.idxs %in% (1:sub.split.case.num)]
      }
      
      split.fname.1 = paste0(sub.wd.sp,"Dt1_Sdel_regr_",in.oc,"_",st,".ts.csv")
      split.fname.2 = paste0(sub.wd.sp,"Dt1_SdelN_regr_",in.oc,"_",st,".ts.csv")
      write.table(sub.ts.idxs, split.fname.1, sep=",",  col.names=FALSE, row.names=FALSE)
      write.table(sub.ts.idxs, split.fname.2, sep=",",  col.names=FALSE, row.names=FALSE)
    }
  }
  return(TRUE);
}

#-------------------------------------------------------------------------------------

if(length(sample.size.arr)==0){
  cat("  subset skipped, no need\n",sep="")
  return()
}

n_len = year.num - forcast.steps

sample.num = length(sample.size.arr)
whole.size = sample.size.arr[sample.num]
whole.sample.dir = paste("syn_sample_",whole.size,"/", sep="")
whole.wd = paste(syn_root_dir,whole.sample.dir,sep="")

# classification
assign("task.type", "class", envir = .GlobalEnv)
for(oc in disease.pred.arr){
  sub.res = gen.subset.data.class(in.oc=oc)
  if(sub.res){
    cat("  Outsome=",oc," subset done.\n",sep="")
  } else {
    cat("  Outsome=",oc," subset skipped.\n",sep="")
  }
}

if(do.regr.tag){
  # regression
  assign("task.type", "regr", envir = .GlobalEnv)
  for(oc in outcome.pred.arr){
    sub.res = gen.subset.data.regr(in.oc=oc)
    if(sub.res){
      cat("  Outsome=",oc," subset done.\n",sep="")
    } else {
      cat("  Outsome=",oc," subset skipped.\n",sep="")
    }
  }
}



