
setwd('D:/Work/2021 Spring/Medcine Reaction final/Experiments/ShareToHaoyu')
#setwd('/home/kumarbio/jiaxx213/workspace/med_react_2021')

# foundamental control !!!! #
#TYPE = "REAL"
TYPE = "SYN"

syn.num = 5
# how many synthetic dataset produced, 
# i.e. number of runs, for averaging results

sample.size.arr = c(500,2500,5000)
sample.size.num = length(sample.size.arr)
whole.size = sample.size.arr[sample.size.num]
#generate the whole.size dataset, than copy subset to be other sample size data

do.regr.tag = F
#whether do the regression task (forecast future labs values) 

if(do.regr.tag){
  # initialized for regr, used in export()
  sugg.pred.outcome.mat = matrix(NA, nrow=pred.slt.num, ncol=sample.size.num);
  colnames(sugg.pred.outcome.mat) = paste("sample",sample.size.arr,sep="_")
  rownames(sugg.pred.outcome.mat) = paste("sugg",1:pred.slt.num,sep="_")
  sugg.pred.outcome.mat[] = NA
}

syn.data <- function(syn.idx, single.ope=NA){
  syn_root_dir = paste("syn_experiments/syn_",syn.idx,"/", sep="")
  assign("syn_root_dir", syn_root_dir, envir = .GlobalEnv)
  if (!dir.exists(syn_root_dir)) {
    dir.create(file.path(syn_root_dir))
  }
  
  # only generate the largest size sample
  assign("pat.num", whole.size, envir = .GlobalEnv)
  
  whole.sample.dir = paste("syn_sample_",whole.size,"/", sep="")
  assign("syn_sample_dir", whole.sample.dir, envir = .GlobalEnv)
  
  whole.wd = paste(syn_root_dir,whole.sample.dir,sep="")
  if (!dir.exists(whole.wd)) dir.create(file.path(whole.wd))
  
  if(is.na(single.ope) || single.ope=="create"){
    cat("[START] Create Dataset: syn_",syn.idx,"\n",sep="")
    # -- create -- #
    source("syn_gen_data.R")
    # -- fit -- #
    source("syn_fit_data.R")
    cat("  syn_",syn.idx," sample ",whole.size," create done.\n",sep="")
    cat("[END] Create Dataset: syn_",syn.idx,"\n",sep="")
    
    source("functions_checking.R")
    check.disease.forecast(whole.wd)
  }
  
  
  if(is.na(single.ope) || single.ope=="export"){
    cat("[START] Export Dataset: syn_",syn.idx,"\n",sep="")

    gen.export.log.flag = paste("  [syn_",syn.idx," | sample ",whole.size,"]",sep="")
    assign("gen.export.log.flag", gen.export.log.flag, envir = .GlobalEnv)

    # -- export -- #
    source("data_Export.R")
    cat("  syn_",syn.idx," sample ",whole.size," export done.\n",sep="")
    if(do.regr.tag){
      #sugg.pred.outcome.mat[,paste("sample_",whole.size,sep="")] = sugg.outcome.pred.arr
      sugg.pred.outcome.mat[] = sugg.outcome.pred.arr # repeat for all size
      print(sugg.pred.outcome.mat)
      write.table(sugg.pred.outcome.mat, syn_suggPreOc_fname, sep=",",  col.names=TRUE, row.names=FALSE)
    }

    cat("[END] Export Dataset: syn_",syn.idx,"\n",sep="")
  }
  
  
  if(is.na(single.ope) || single.ope=="split"){
    source("functions_common.R")
    cat("[START] Split Dataset: syn_",syn.idx,"\n",sep="")
    # -- split -- #
    source("data_Split.R")
    cat("  syn_",syn.idx," sample ",whole.size," split done.\n",sep="")
    cat("[END] Split Dataset: syn_",syn.idx,"\n",sep="")
  }
  
  if(is.na(single.ope) || single.ope=="subset"){
    source("functions_common.R")
    cat("[START] Subset Produce: syn_",syn.idx,"\n",sep="")
    # -- subset -- #
    source("data_Subset.R")
    cat("[END] Subset Produce: syn_",syn.idx,"\n",sep="")
  }
  
}

syn.pred <- function(syn.idx, pred.param.arr){
  # [ FOLDED BY JIA. Apr2020 ]
}

#-----------------------------------------------------------

#== Entrance Here!!==#
#args = commandArgs(trailingOnly=TRUE)
args=c("gen")


args.num = length(args)
if(!args[1] %in% c("gen","pred")){
  cat("Main_Syn.R gen|pred [NULL,create,export,split,subset]|[syn.idx]\n")
  quit()
}

if(args[1]=="gen"){
  if(args.num==1)
    gen.param = c(NA)
  else {
    gen.param = args[2:args.num]
  }
  for (p in gen.param) {
    for (i in 1:syn.num) {
      syn.data(i,p)
    }
  }
}

if(args[1]=="pred"){
  in.syn.idx = as.numeric(args[2])
  in.pred.param.arr = args[3:args.num]
  if(is.na(in.syn.idx) || !in.syn.idx %in% 1:syn.num){
    cat("Wrong Input! Main_Syn.R pred syn.idx=1..",syn.num,"[param of pred_Linear.R] \n",sep="")
    quit()
  }
  syn.pred(in.syn.idx, in.pred.param.arr)
}


