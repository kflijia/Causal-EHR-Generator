# Causal-EHR-Generator


***************** Files *******************

[Main_Syn.R]  entrance main()
option: synthetic data generation
	real data export (folded)
	synthetic data prediction tasks (folded)
	real data prediction tasks (folded)

[syn_gen_data.R]  synthetic data generation
[functions_syn_gen_data.R] functions of it

[syn_fit_data.R]  change syn data structure to fit my generalized export portal

[data_Export.R]  data export
[data_Split.R]  data split training/test
[data_Subset.R]  copy subset as smaller sample size data

[setting_files.R]  manage all file paths, global influnce
[setting_paramts_data.R]  parameters setting for data generation, global influnce


[functions_common.R] [functions_basic.R] [functions_utils.R]
useful basic functions 

[functions_checking.R]  all functions for checking the generated data
			including display trajectory figures

[syn_experiments]  the folder of saving all files for syn data experiments



***************** Procedure *******************
Synthetic Data Generation:
stage 1: create, run [syn_gen_data.R]->[syn_fit_data.R]
stage 2: export, export datasets as .csv files
stage 3: split, randomize training/test split and save as .csv files
stage 4: subset, produce smaller sample size datasets



***************** Usage *******************
In Rstudio: edit input parameters on [Main_Syn.R] Line#106
By Linux Command: "Rscript Main_Syn.R xxx "
(Uncomment Main_Syn.R Line#105 and comment Line#106)

Command Example:
Rscript Main_Syn.R gen   # generating all required synthetic data
Rscript Main_Syn.R gen create  #  only do create stage
Rscript Main_Syn.R gen export  #  only do export stage
Rscript Main_Syn.R gen split  #  only do split stage
Rscript Main_Syn.R gen subset  #  only do subset stage

Required Parameters setting:
setwd('xxxxx')  [Main_Syn.R] first line, set your locate
sample.size.arr  [Main_Syn.R] 
pat.num  [setting_paramts_data.R] number of patients, i.e. samples
labs.num  [setting_paramts_data.R] number of labs
year.num  [setting_paramts_data.R] number of observation years
combd.num  [setting_paramts_data.R] number of disease, named as combdX
disease.pred.arr  [setting_paramts_data.R] the diseases your want to predict
forcast.steps  [setting_paramts_data.R] number of future steps you want to predict on
steps.arr  [setting_paramts_data.R] looking-back steps you want in predictions
split.ts.prop  [setting_paramts_data.R] propotion of selecting test samples from whole data
obs.prob  [setting_paramts_data.R] probability of being observed


***************** Display *******************
Uncomment [functions_checking.R] Line#435-#444 and setting it as needed
source("functions_checking.R")
display.pat(in.pat, in.labs=NA, onlymajor=F)

output example:

> display.pat(in.pat = 34, in.labs=NA, onlymajor=F)   # display patient 34

combd rout: ( 1 -> ) 3->4->8->10   # disease trajectory, i.e routine of combds, in () is not under observation yet
major labs: ( 6 -> ) 18->8->10->2   # corresponding major-influened-lab along routine
minor labs: ( [18,1]-> ) [8,7,9]->[7,20,10,25]->[16,2,3,13]->[12,21]   # minor-influened-labs along routine
combd dx years: 1, 15, 11, 26   # the years when diseases diagnosed

-------- develop --------
Labs = 18   # which lab to display, set by in.labs, by default is the first major-influened lab
Major position = 1/4   # the displayed lab is majorly influenced by which disease in routine
Minor position = 3/4  # the displayed lab is minorly influenced by which disease in routine
Harmful position No.1 = 2/4  # which disease's medicine can be harmful to the displayed lab, ordered No.[1..n]

-------- medicine --------
Major med.line=1,2,3,4 position=2,4,20,56/80 Special: TRUE,FALSE,FALSE,FALSE
Minor No.1 med.line=3 position=20/80
Harmful No.1 med.line=1 position=52/80

med.line & position: used medicine lines, and happening positions in the data serial (len(data.serial)=year.num*2)
Special: whether this patient is specially reacted to the medicine, may useless even harmful as contrary
Minor No.1: lines & positions of the 1st medicine that minorly influent the lab
Harmful No.1: lines & positions of the 1st medicine that harmfully influent the lab 

In the figure, filled dots=observed values, unfilled dots=unobserved values






