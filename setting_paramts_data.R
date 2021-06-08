

source("functions_basic.R")
source("functions_utils.R")

basic.arr = c("age","male","black","alive")

randomfill.tag = F
# True: all unknown delta filled random noise,
# including all delta in wo.output.data and unknown delta in wt.output.data,
# and the unknown part keeps consistant between wo and wt data
# False: all unknown delta filled 0

# ==============================#


if(TYPE=="REAL") {  #  Real Data
  
  # [ FOLDED BY JIA. Apr2020 ]

} else if(TYPE=="SYN"){  # Synthetic Data
  
  #---- generate -----#
  if(!exists("pat.num")) pat.num = 300  # i.e. sample.size, initialized as 300, controled in Main()
  if(!exists("labs.num")) labs.num = 20  # initialized, may covered by loading old data
  if(!exists("year.num")) year.num = 40  # initialized, may covered by loading old data
  if(!exists("combd.num")) combd.num = 10  # initialized, may covered by loading old data
  #rout.num = 10
  
  wide.prop = 5 
  # used in rand.combd() & init.value(), spread 1/wide.prop mean value on each side
  
  fold.num = 4
  #because pat routine randomly start at any stage, 
  #year.num have estimated 1/fold.num chance to see full routine
  #---- fit & export -----#
  
  disease.arr = paste("combd",1:combd.num,sep="");
  interv.arr = paste(disease.arr, "_ndx", sep="");
  disease.pred.arr = disease.arr[combd.num] # only the last one
  interv.pred.arr = interv.arr; # initialized
  
  outcome.arr = paste("labs",1:labs.num,sep="");
  outcome.pred.arr = outcome.arr; # initialized
  outcome.missing.arr = paste("labs",1:labs.num,"_missing",sep="");
  
  pred.slt.num = 4; # top largest dayaset is selected. 0:disable
  
  forcast.steps = 10; #5
  
  years.per.step = 1; # only used once in export
  
  split.ts.prop = 0.3

  steps.arr = c(5, 10, 15) # must <= (year.num-forecast.num)
  
  obs.prob = 0.65 # probability of observed
  
  if(!exists("gen.export.log.flag")) gen.export.log.flag="  [synthetic]"
}









