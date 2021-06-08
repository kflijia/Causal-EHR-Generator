source("setting_files.R")
source("setting_paramts_data.R")

if(TYPE=="REAL"){
  #------- real --------#
  load(glb_paramt_fname)
  #pat.num
  
  year.num = t.num # consistantly called 'year' from now on 
  
  load(real_data_fname)
  #basic.data.list, interv.data.list, meds.data.list, outcome.data.list
  
  data_dir = real_data_dir
  
  pid.vec = basic.data.list[[1]]$pid

} else if(TYPE=="SYN"){
  #------ synthetic -------#
  
  load(syn_data_fname)
  #basic.data.list, interv.data.list, meds.data.list, outcome.data.list
  
  load(ideal_delta_fname)
  #delta.trj.list: outcomes.num length list, pat.num * year.num matrix delta value
  #delta.outcome.list: year.num length list, pat.num * outcomes.num matrix delta value
  
  data_dir = syn_data_dir
  
  pid.vec = as.numeric(rownames(basic.data.list[[1]]))
}

if(!dir.exists(data_dir)) dir.create(file.path(data_dir))

#---------------------------------------#

compress.meds <- function(in.med.data){
  out.med.data = c()
  for(c in 1:combd.num){
    my.prex = paste0("med_combd",c,"_line")
    my.long.data = in.med.data[,grepl(my.prex,colnames(in.med.data),fixed=T),drop=F]
    my.short.data = apply(my.long.data, 1, sum)
    out.med.data = cbind(out.med.data, my.short.data)
  }
  colnames(out.med.data) = paste0("med_combd",1:combd.num)
  return(out.med.data)
}

scale.mat = as.data.frame(matrix(0, 0, 3))
colnames(scale.mat) = c("y","type","count")

output.csv.data <- function(data.type, delta.type, task.type, the.y){
  # checking validation
  validation.arr = rep(0, 4)
  names(validation.arr) = c("data.type", "delta.type", "task.type", "the.y")
  
  if(! data.type %in% c("Dt1","Dt2")) validation.arr[["data.type"]] = 1
  if(! delta.type %in% c("Sdel","Adel")) validation.arr[["delta.type"]] = 1
  if(! task.type %in% c("regr","class") ) validation.arr[["task.type"]] = 1
  if(! the.y %in% c(outcome.arr, disease.pred.arr)) validation.arr[["the.y"]] = 1
  if(sum(validation.arr)>0) {
    print(validation.arr);
    return();
  }
  
  if(task.type=="regr" && !(the.y %in% outcome.arr) ) validation.arr[c("task.type","the.y")] = 2
  if(task.type=="class" && !(the.y %in% disease.pred.arr) ) validation.arr[c("task.type","the.y")] = 2
  if(sum(validation.arr)>0) {
    print(validation.arr);
    return();
  }
  
  if(task.type=="class" && delta.type=="Sdel" ) validation.arr[c("task.type","delta.type")] = 3
  if(sum(validation.arr)>0) {
    print(validation.arr);
    return();
  }
  
  meds.arr = colnames(meds.data.list[[1]])
  meds.out.arr = paste0("med_combd",1:combd.num)
  
  wt.output.fname = paste(data.type,"_",delta.type,"_",task.type,"_",the.y,".data.csv",sep="")
  wt.output.fname = paste(data_dir,"/",wt.output.fname,sep="")
  wo.output.fname = paste(data.type,"_",delta.type,"N_",task.type,"_",the.y,".data.csv",sep="")
  wo.output.fname = paste(data_dir,"/",wo.output.fname,sep="")
  
  if(file.exists(wt.output.fname) && file.exists(wo.output.fname)) {
    cat("    output.csv.data() Skip because existed:\n")
    cat("     ",wt.output.fname,"\n")
    cat("     ",wo.output.fname,"\n")
    con<-file(wt.output.fname)
    the.count = length(readLines(con))-1
    scale.mat[nrow(scale.mat)+1,] = c(the.y, task.type, the.count)
    assign("scale.mat", scale.mat, envir = .GlobalEnv)
    return()
  }
  
  # start producing data
  
  if(delta.type=="Adel"){
    outcome.trj.mat = c()
    for(the.oc in outcome.arr){
      outcome.trj = delta.trj.list[[the.oc]]
      outcome.trj.mat = cbind(outcome.trj.mat, outcome.trj)
    }
    trj.available.th = 3
  } else { # Sdel, the.y is from outcome.arr
    outcome.trj.mat = delta.trj.list[[the.y]]
    trj.available.th = 2
  }
  
  slt.idxs = which(apply(outcome.trj.mat, 1, function(x){sum(!is.na(x))}) >= trj.available.th)
  slt.count = length(slt.idxs)
  cat("   ",the.y,"number of cases:",slt.count,"\n")
  # record scare, to select y.pred.arr
  scale.mat[nrow(scale.mat)+1,] = c(the.y, task.type, slt.count)
  assign("scale.mat", scale.mat, envir = .GlobalEnv)
  
  
  if(slt.count==0){
    cat("    Warning: ",the.y,"no case! \n")
    write.table(c(), wt.output.fname)
    write.table(c(), wo.output.fname)
    return()
  }
  
  if(data.type=="Dt1") variables.arr = c(outcome.arr, meds.out.arr)
  else variables.arr = c(basic.arr, disease.arr, outcome.arr, meds.out.arr) # Dt2
  if(delta.type=="Adel") { # regression or classification
    delta.arr = paste("delta_",outcome.arr,sep="")
  }
  else { # Sdel, must be regression
    delta.arr = paste("delta_",the.y,sep="") 
  }
  y.arr = "y"
  
  wt.output.data = matrix(0, nrow = length(slt.idxs)*n_len, ncol = length(c(variables.arr,delta.arr,y.arr))+1 )
  colnames(wt.output.data) = c("pid",variables.arr,delta.arr,y.arr)

  # if year.num=30, forecast=5, then each case consume 6 records, 25 cases produced in total
  for(i in 1:n_len) {
    my.idx = (1:length(slt.idxs) - 1) * n_len + i
    
    my.pid = pid.vec[slt.idxs]
    my.pid = as.matrix(my.pid)
    my.bf = outcome.data.list[[i]]$bf[slt.idxs,,drop=F]
    my.bf = as.matrix(my.bf)
    my.meds = meds.data.list[[i]][slt.idxs,meds.arr,drop=F]
    my.meds = as.matrix(my.meds)
    my.out.meds = compress.meds(my.meds)
    
    wt.output.data[my.idx, "pid"] = my.pid
    wt.output.data[my.idx, outcome.arr] = my.bf
    wt.output.data[my.idx, meds.out.arr] = my.out.meds
    
    if(data.type=="Dt2"){
      my.basic = basic.data.list[[i]][slt.idxs,basic.arr,drop=F]
      my.basic = as.matrix(my.basic)
      my.basic[,"age"] = my.basic[,"age"] + (i-1)*years.per.step; 
      my.interv = interv.data.list[[i]][slt.idxs,disease.arr,drop=F] 
      # Notes: interv.data.list not applied pat.time.obs.mat, all observible
      my.interv = as.matrix(my.interv)
      
      wt.output.data[my.idx, basic.arr] = my.basic
      wt.output.data[my.idx, disease.arr] = my.interv
    }
    
    if(task.type=="regr"){
      my.y = outcome.data.list[[i+forcast.steps]]$af[[the.y]][slt.idxs]
      my.y = as.matrix(my.y)
    } else { # class
      my.y = interv.data.list[[i+forcast.steps]][[the.y]][slt.idxs] 
      my.y = as.matrix(my.y)
    }
    
    wt.output.data[my.idx, "y"] = my.y
    #    wo.output.data[my.idx, "y"] = my.y
    
    if(delta.type=="Adel"){
      my.delta = delta.outcome.list[[i]][slt.idxs,outcome.arr, drop=F]
    } else { # Sdel, must be regression, the.y is from outcome.arr
      my.delta = delta.outcome.list[[i]][slt.idxs,the.y, drop=F]
    }
    
    wt.output.data[my.idx, delta.arr] = my.delta
  }
  wt.output.data[is.na(wt.output.data)] = 0  # for REAL
  
  #------- changed here! -------
  wo.output.data = wt.output.data
  wo.output.data[,delta.arr] = 0
  if(randomfill.tag){
    wo.output.data[,delta.arr] = runif(nrow(wo.output.data)*length(delta.arr),0,1) * 0.1
    rdfill.stamp = is.na(wt.output.data) # get same shape of all F
    rdfill.stamp[,delta.arr] = wt.output.data[,delta.arr]==0  # delta default is 0
    wt.output.data[rdfill.stamp] = wo.output.data[rdfill.stamp]
  }
  
  write.table(wt.output.data, wt.output.fname, sep=",",  col.names=TRUE, row.names=FALSE)
  write.table(wo.output.data, wo.output.fname, sep=",",  col.names=TRUE, row.names=FALSE)
  
  
}

slt.scale.4.pred <- function(in.type, in.pred){
  if(is.numeric(pred.slt.num) && pred.slt.num>0){
    tmp.mat = scale.mat[scale.mat[["type"]]==in.type, c("y","count")]
    tmp.idx = order(as.numeric(tmp.mat[["count"]]), decreasing = T)[1:pred.slt.num]
    print(tmp.mat[tmp.idx,]);
    tmp.idx = sort(tmp.idx)
    out.arr = tmp.mat[tmp.idx,"y"]
    if(length(out.arr)==0) cat("    No y slt for tasks!!\n")
    return(out.arr)
  } else {
    return(in.pred)
  }
}

#------------------------------------------------------------------------------------------------

n_len = year.num - forcast.steps
# max number of records per patient can produce

for(oc in disease.pred.arr) {
  output.csv.data(data.type="Dt1", delta.type="Adel", task.type="class", the.y=oc)
}
cat(gen.export.log.flag,"    classification export done. Labels = ");
print(paste(disease.pred.arr, collapse = ", "))
cat("\n")

if(do.regr.tag){
  for(oc in outcome.pred.arr){
    output.csv.data(data.type="Dt1", delta.type="Sdel", task.type="regr", the.y=oc)
  }
  sugg.outcome.pred.arr = slt.scale.4.pred("regr", outcome.pred.arr)
  assign("sugg.outcome.pred.arr", sugg.outcome.pred.arr, envir = .GlobalEnv)
  cat(gen.export.log.flag,"   regression export done. Suggested Labels = ");
  cat(paste(sugg.outcome.pred.arr, collapse = ", "))
  cat("\n")
}


source("functions_checking.R")
cat("Checking Delta Rates:\n")
for(oc in disease.pred.arr) {
  check.delta.rate(data.type="Dt1", delta.type="Adel", task.type="class", the.y=oc)
}
if(do.regr.tag){
  for(oc in outcome.pred.arr) {
    check.delta.rate(data.type="Dt1", delta.type="Sdel", task.type="regr", the.y=oc)
  }
}








