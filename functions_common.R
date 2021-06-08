#------for Prediction Stage only------#
#------And also data Split+Subset! --------
#------also checking... --------

# produce looking-back length for step.num experiments
# mode=0: average interval. mode=1: 1 step interval. mode=2: 5 step interval.
produce.steps <- function(step.num=4,in.mode=0){
  if(in.mode==1) { # 1 step interval
    steps.arr = 1:step.num
  } else if(in.mode==2) { # 5 step interval
    step.base = 4
    step.top = n_len
    step.interval = 5
    steps.arr = c(step.base+1, step.base+1+seq(1:(step.num-1))*step.interval)
    exceed.idx = steps.arr>step.top
    if(sum(exceed.idx)>0){
      steps.arr = c( steps.arr[!exceed.idx], step.top)
    }
  } else { # in.mode==0, by default, take average interval
    step.base = 4
    step.top = round((n_len+1)*(2/3))
    step.interval = (step.top-step.base)/(step.num-1)
    steps.arr = round(c(step.base+1, step.base+1+seq(1:(step.num-2))*step.interval, step.top))
  }
  return(steps.arr)
}

library(dplyr)


get.fname <- function(ftype, D.tag="", model.name=""){
  
  if( (D.tag)!="" && !(D.tag %in% c("wt","wo")) ){
    print("Wrong Delta.tag, wt / wo")
    return("")
  }
  
  if( model.name!="" && !(model.name %in% c("glm","glmer","lm","lmer","gru","lstm")) ){
    print("Wrong model.name, [glm|glmer|lm|lmer|gru|lstm]")
    return("")
  }
  
  data.type = "Dt1"
  
  if (task.type=="class") delta.type="Adel"
  if (task.type=="regr") delta.type="Sdel"
  
  if(!exists("task.type")) task.type="undefined_taskType"
  if(!exists("delta.type")) delta.type="undefined_deltaType"
  if(!exists("task.outcome")) task.outcome="undefined_taskOC"
  if(!exists("task.step")) task.step="undefined_taskStep"
  
  # only using with-delta data file, defined in Data Stage
  data.fname = paste(paste(data.type,delta.type,task.type,task.outcome,sep="_"),"data","csv",sep=".")
  split.fname = paste(paste(data.type,delta.type,task.type,task.outcome,task.step,sep="_"),"ts","csv",sep=".")
  serial.fname = paste(paste(data.type,delta.type,task.type,task.outcome,sep="_"),"serial","rda",sep=".")
  updated.fname = paste(paste(task.type,delta.type,task.outcome,task.step,sep="_"),"updated","csv",sep=".")
  
  output.fname = paste(paste(task.type,D.tag,task.outcome,task.step,model.name,sep="_"),"output","csv",sep=".")
  pvalue.fname = paste(paste(task.type,task.outcome,task.step,model.name,sep="_"),"pvalue","csv",sep=".")
  
  # just for going through
  if(!exists("data_dir")) data_dir="undefined_dataDir";
  if(!exists("split_dir")) split_dir="undefined_splitDir";
  if(!exists("updated_dir")) updated_dir="undefined_updatedDir";
  if(!exists("output_dir")) output_dir="undefined_outputDir"; 
  if(!exists("pvalue_dir")) pvalue_dir="undefined_pvalueDir";
  
  out.fname = case_when(
    ftype=="data"  ~ paste0(data_dir,"/",data.fname),
    ftype=="split" ~ paste0(split_dir,"/",split.fname),
    ftype=="serial" ~ paste0(split_dir,"/",serial.fname), # share same dir
    ftype=="update"  ~ paste0(updated_dir,"/",updated.fname),
    ftype=="output" ~ paste0(output_dir,"/",output.fname),
    ftype=="pvalue" ~ paste0(pvalue_dir,"/",pvalue.fname),
    TRUE ~ ""
  )
  
  if(out.fname==""){
    print("Wrong ftype, [data|split|update|output|pvalue]")
    return("")
  }
  
  return(out.fname)
}

#--------------------------

find.serial.cases <- function(in.step){
  # prerequisite: task.type, task.outcome
  # total.tag return all idxs
  data.fname = get.fname("data");
  serial.fname = get.fname("serial");
  
  serial.list = list()
  my.st.name = paste0("st",in.step);
  if(file.exists(serial.fname)) {
    load(serial.fname);
    exist.steps = names(serial.list)
    if(my.st.name %in% exist.steps){
      return(serial.list[[my.st.name]]);
    }
  }

  data = read.csv(data.fname, header = T)
  devlp.rm.tag = F;
  if(task.type=='class') devlp.rm.tag = T;
  # devlp.rm.tag = F not checking developed.tag

  X.idx.list = vector("list", in.step)
  Y.vec = c()
  for(i in seq(1,nrow(data),n_len)){
    developed.tag = F;
    for(j in i:(i+n_len-2)){
      if(developed.tag) break;
      #if((j+in.step-1)>(i+n_len-1-(forcast.steps-1))) break; # (j+in.step-1) is bordon of X
      if((j+in.step-1)>(i+n_len-1)) break; # equal to forcast.steps=1
      my.X.idx = j:(j+in.step-1)
      my.pid = data[my.X.idx,1];
      if(length(unique(my.pid))>1) {
        cat(pred.log.flag," No unique pid: step=",in.step,", i=",i,", j=",j,"\n",sep="")
        return()
      }
      # for both class and regr
      # X must not containe steps when y has devlp
      last.step.orig.idx = my.X.idx[in.step]-forcast.steps
      if(last.step.orig.idx >= my.X.idx[1]){  # have available y recorded in exported data
        last.step.y = data[last.step.orig.idx,ncol(data)]
        if(last.step.y>=1) break; # has devlp
      }
      
      #my.Y = data[(j+in.step-1+forcast.steps-1),ncol(data)]
      my.Y = data[(j+in.step-1),ncol(data)]
      if(my.Y<=-999) next; # -Inf, marked NA
      if(devlp.rm.tag && my.Y>=1){
        my.Y=1;
        developed.tag = T;
      }
      for (s in 1:in.step) {
        X.idx.list[[s]] = c(X.idx.list[[s]], my.X.idx[s])
      }
      Y.vec = c(Y.vec, my.Y)
    }
  }
  
  serial.list[[my.st.name]] = list(X.idx.list=X.idx.list, Y.vec=Y.vec);
  save(serial.list, file = serial.fname)
  
  return(serial.list[[my.st.name]])
}






