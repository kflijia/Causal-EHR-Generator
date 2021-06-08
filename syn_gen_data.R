
source("setting_files.R")
source("setting_paramts_data.R")

#debugSource("functions_syn_gen_data.R")
source("functions_syn_gen_data.R")

# ============= structure data generating/loading ============== 

if(file.exists(gen_struct_fname)){
  # year.num using default value, make sure consistant with loaded structure
  cat("  Loading",gen_struct_fname,"...\n")
  load(gen_struct_fname);
  #labs.list, combd.major.list, combd.minor.list, rout.list, med.list
  labs.num = length(labs.list$labs.mean.arr)
  combd.num = length(combd.major.list)
  rout.mat = rout.list$rout.mat
  rout.num = nrow(rout.mat)
  max.step = ncol(rout.mat)
  
  assign("labs.num", labs.num, envir = .GlobalEnv)
  assign("combd.num", combd.num, envir = .GlobalEnv)
  # cover default setting
  
  tmp.devlp.mat = rout.list$rout.major.devlp.mat
  tmp.devlp.mat[is.na(tmp.devlp.mat)] = 0
  longest.rout.years = max(apply(tmp.devlp.mat, 1, sum))
  
  medicine.mat = med.list$medicine.mat
  
} else {
  cat("  Making",gen_struct_fname,"...\n")
  
  # randomize labs
  labs.list = rand.labs(labs.num)
  
  # randomize routines
  connect.rate = 0.5
  rout.list = rand.routine() # combd.num, connect.rate
  rout.graph = rout.list$out.graph
  rout.mat = rout.list$rout.mat
  rout.num = nrow(rout.mat)
  max.step = ncol(rout.mat)
  print(rout.mat)
  
  # randomize comobidities
  combd.res = rand.combd() #labs.list, rout.list, wide.prop
  labs.list = combd.res[[1]]
  combd.major.list = combd.res[[2]]
  combd.minor.list = combd.res[[3]]
  
  # randomize severty
  rout.list = rand.severity() #rout.list, year.num, combd.minor.list
  rout.mat = rout.list$rout.mat
  rout.major.devlp.mat = rout.list$rout.major.devlp.mat
  rout.minor.devlp.mat = rout.list$rout.minor.devlp.mat
  longest.rout.years = max(apply(rout.major.devlp.mat, 1, max))
  cat("longest routine #years =",longest.rout.years,"\n")
  #rout.list: rout.mat, connect.mat, rout.major.devlp.mat, rout.major.sevrt.mat, rout.minor.sevrt.mat
  
  # randomize medicine
  harmful.param = c(1, 4, 0.2, 0.8)
  names(harmful.param) = c("min.num","max.num","min.prob","max.prob")
  max.line = 5
  med.list = rand.medicine() # labs.num, combd.major.list, combd.minor.list, max.line, harmful.param
  medicine.mat = med.list$medicine.mat
  max.line = ncol(medicine.mat)
  # medicine.mat, medicine.effect.mat
  
  # structure saving
  save(labs.list, combd.major.list, combd.minor.list, rout.list, med.list, file = gen_struct_fname)
}
# ============= patients data generating/loading ============== 

if(file.exists(gen_patient_fname) && file.exists(ideal_delta_fname)) {
  cat("    gen_patient() Skip because existed:\n")
  cat("     ",gen_patient_fname,"\n")
  cat("     ",ideal_delta_fname,"\n")
  load(gen_patient_fname)
  load(ideal_delta_fname)
  
  year.num = ncol(pat.time.obs.mat)/2
  assign("year.num", year.num, envir = .GlobalEnv)
  # cover default setting
  
} else {
  # initialize no disease population
  init.res = init.value()
  # pat.num, labs.list, year.num, wide.prop, obs.prob
  value.init.list = init.res[[1]]
  pat.time.obs.mat = init.res[[2]] # used in get.treatment
  
  # randomize personal severity
  pat.rout.list = rand.pat.severity() # rout.list
  pat.rout.assigned = pat.rout.list$pat.rout.assigned
  pat.start.assigned = pat.rout.list$pat.start.assigned
  pat.start.offset = pat.rout.list$pat.start.offset
  pat.sevrt.rate = pat.rout.list$pat.sevrt.rate
  
  # assigne routine and disease to patients
  disease.res = get.disease()
  # value.init.list, combd.major.list, combd.minor.list, rout.list, pat.rout.list
  value.disease.list = disease.res[[1]]
  maps.list = disease.res[[2]]
  
  # inner patient saving
  #save(value.init.list, pat.time.obs.mat, pat.rout.list, value.disease.list, file = gen_patient_fname)
  
  # randomize personal medicine effect
  special.rate = 0.2
  med.pat.effect.mat = rand.pat.effect() # medicine.mat, pat.num, special.rate
  
  # assigne treatment to patients
  crit.increase.rate = 0.1
  treatment.res = get.treatment()
  # value.init.list, pat.time.obs.mat, combd.major.list, combd.minor.list, rout.list, 
  # labs.list, maps.list, med.list, med.pat.effect.mat, pat.rout.list, crit.increase.rate
  value.treated.list = treatment.res[[1]]
  pat.time.obs.mat = treatment.res[[2]] # updated
  pat.dx.mat = treatment.res[[3]]
  pat.med.mat = treatment.res[[4]]
  treatment.list = treatment.res[[5]]
  #pat.labs.obs.mat = treatment.res[[6]]
  delta.trj.list = treatment.res[[6]]
  delta.outcome.list = treatment.res[[7]]
  
  #---------save patients data & delta here -------#
  
  save(value.init.list, pat.time.obs.mat, pat.rout.list, value.disease.list, maps.list, med.pat.effect.mat,
       value.treated.list, pat.dx.mat, pat.med.mat, treatment.list, #pat.labs.obs.mat,
       file = gen_patient_fname)
  
  save(delta.trj.list, delta.outcome.list, file = ideal_delta_fname)
}

######## check data generation here #######


#source("functions_checking.R")
#display.pat(in.pat = 3, in.labs=NA, onlymajor=F)


