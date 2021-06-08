
source("setting_files.R")
source("setting_paramts_data.R")

load(gen_struct_fname)
load(gen_patient_fname)
#value.init.list: len=labs.num, element=matrix[pat.num * (year.num*2)]
#labs.list: labs.mean.arr[len=labs.num], labs.sd.arr[len=labs.num]
#pat.time.obs.mat: [pat.num * (year.num*2)]
#combd.major.list: len=combd.num, element=array(names: labs, rapid, devlp, crit)
#combd.minor.list: len=combd.num, element=list(names: labs, rapid, devlp, crit), consistent len of all 4 array, may len=0 
#rout.list: rout.major.mat[rout.num * max.step], connect.mat[(combd.num+1) * (combd.num+1)]
#           rout.major.sevrt.mat[rout.num * max.step]
#           rout.minor.sevrt.mat[rout.num * max.step], element=array[may len=0]
#med.list: medicine.mat[combd.num * max.line], element is NULL or list(major.labs[len=1],minor.labs[may len=0],harmful.labs[may len=0])
#          medicine.effect.mat[combd.num * max.line], element is NULL or list(major.effect[len=1],minor.effect[may len=0],harmful.effect[may len=0])
#med.pat.effect.mat: [combd.num * max.line], element is NULL or array[len=pat.num]
#value.disease.list: len=labs.num, element=matrix[pat.num * (year.num*2)]
#pat.rout.list: rout.assigned[len=pat.num]
#               pat.major.rout.sevrt[pat.num * max.step]
#               pat.minor.rout.sevrt[pat.num * max.step], element=array[may len=0]
#               major.altr.devlp.mat[pat.num * max.step], may = NA
#value.treated.list: len=labs.num, element=matrix[pat.num * (year.num*2)]
#pat.dx.mat: [pat.num * year.num], element is NULL(healthy pat) or array(len=combd.num), binary
#pat.med.mat: [pat.num * year.num], element is NULL(healthy pat) or array(len=combd.num), numeric
#treatment.list: len=year.num, element=matrix[pat.num * combd.num] 0=no meds, 1-max.line=medicine line


#-------------------------------------------

medicine.mat=med.list$medicine.mat
#[combd.num * max.line], element is NULL or list(major.labs[len=1],minor.labs[may len=0],harmful.labs[may len=0])
combd.num = nrow(medicine.mat)
max.line = ncol(medicine.mat)
medicine.avlb.mat = !matrix(sapply(medicine.mat, is.null),combd.num,max.line)
meds.arr = c()
meds.max.arr = c()
for (i in 1:combd.num) {
  my.max = sum(medicine.avlb.mat[i,])
  for (j in 1:my.max){
    my.med.name = paste("med_combd",i,"_line",j,sep="")
    meds.arr = c(meds.arr, my.med.name)
  }
  meds.max.arr = c(meds.max.arr, my.max)
}
meds.num = length(meds.arr)


fit.data.list <- function(){
  # ----- randomize basic info -------
  basic.info.mat = as.data.frame(matrix(0, pat.num, length(basic.arr)))
  colnames(basic.info.mat) = basic.arr
  
  age.samples = rnorm(pat.num*5, 0, 1) # random
  age.samples = age.samples[age.samples>=-0.7 & age.samples<=2] # cutoff
  age.samples = round((age.samples-min(age.samples))/(max(age.samples)-min(age.samples))*80)+18
  
  basic.info.mat$age = age.samples[1:pat.num]
  basic.info.mat$male = sample(0:1, pat.num, T)
  basic.info.mat$black = sample(0:1, pat.num, T, prob=c((1-0.026),0.026))
  basic.info.mat$alive = rep(1, pat.num)
  
  basic.data.list=list()
  for(i in 1:year.num){
    basic.data.list[[i]] = basic.info.mat
  }
  
  # --------med + interv(dx) data ----------
  meds.extend <- function(in.meds.mat) {
    # in.meds.mat ncol is combd.num
    out.meds.mat = c()
    for(c in 1:length(meds.max.arr)){
      my.max = meds.max.arr[c]
      my.meds = in.meds.mat[,c]
      my.mat = matrix(0,nrow(in.meds.mat),my.max)
      for(m in 1:my.max){
        my.mat[my.meds>=m, m] = 1
      }
      out.meds.mat = cbind(out.meds.mat, my.mat)
    }
    return(out.meds.mat)
  }
  
  meds.data.list = vector("list",year.num);
  interv.data.list = vector("list",year.num);
  pat.illed = unlist(lapply(pat.med.mat[,1], length))>0 # assigned routine > 0
  for(i in 1:year.num){
    meds.mat = matrix(0, pat.num, combd.num) # healthy pat: all 0
    dx.mat = matrix(0, pat.num, combd.num) # healthy pat: all 0
    for(j in which(pat.illed)){
      meds.mat[j,] = pat.med.mat[[j,i]]
      dx.mat[j,] = pat.dx.mat[[j,i]]
    }
    ext.meds.mat = meds.extend(meds.mat)
    ext.meds.mat = as.data.frame(ext.meds.mat)
    colnames(ext.meds.mat) = meds.arr
    meds.data.list[[i]] = ext.meds.mat
    dx.mat = as.data.frame(dx.mat)
    colnames(dx.mat) = disease.arr   # simple names, use disease.arr
    
    # # unobservible dx
    # #pat.time.obs.mat: [pat.num * (year.num*2)]  T/F
    # pat.year.obs.vec = pat.time.obs.mat[,(i*2-1):(i*2),drop=F]
    # pat.year.obs.vec = apply(pat.year.obs.vec, 1, sum)>0
    # dx.mat[!pat.year.obs.vec,] = NA
    interv.data.list[[i]] = dx.mat
  }
  
  # --------outcome data ----------
  # obs random function
  # randomly flip T/F
  rand.obs <- function(in.obs.vec){
    vec.len = length(in.obs.vec)
    flip.num = round(vec.len*obs.random.rate)
    flip.idxs = sort(sample(1:vec.len, flip.num))
    out.obs.vec = in.obs.vec
    out.obs.vec[flip.idxs] = !out.obs.vec[flip.idxs]
    return(out.obs.vec)
  }
  # outcome 
  outcome.data.list = vector("list",year.num);
  init.value.mat = as.data.frame(matrix(NA,pat.num,labs.num))
  colnames(init.value.mat) = outcome.arr
  for(i in 1:year.num){
    outcome.data.list[[i]]= list(bf=init.value.mat, af=init.value.mat)
  }
  for(j in 1:labs.num){
    # debug
    if(j==3){
      cat("")
    }
    
    labs.value.mat = value.treated.list[[j]]; # pat.num * (year.num*2)
    #labs.obs.vec = pat.labs.obs.mat[,j] # start of time, numeric, maybe NA
    for(i in 1:year.num){
      bf.vec = labs.value.mat[,i*2-1]
      bf.obs.vec = pat.time.obs.mat[,i*2-1]; # T / F
      bf.vec[!bf.obs.vec] = NA # unobservable
      outcome.data.list[[i]]$bf[,j] = bf.vec
      
      af.vec = labs.value.mat[,i*2]
      af.obs.vec = pat.time.obs.mat[,i*2]; # T / F
      af.vec[!af.obs.vec] = NA; # unobservable
      outcome.data.list[[i]]$af[,j] = af.vec
    }
  }
  
  save(basic.data.list, interv.data.list, meds.data.list, outcome.data.list, file = syn_data_fname)
  
}

#-----------------------------------------------------------------------------------------------------

if(file.exists(syn_data_fname)) {
  cat("    fit.data() Skip because existed:\n")
  cat("     ",syn_data_fname,"\n")
  return()
}

fit.data.list()

