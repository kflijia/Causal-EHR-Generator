
# --  set on Main() --
#setwd('D:/Work/2021 Spring/Medcine Reaction final/Experiments')
#setwd('/home/kumarbio/jiaxx213/workspace/med_react_2021')

#=============== real fnames setting =================#

if(!exists("real_root_dir")) real_root_dir = "./";  # default
if(!exists("real_sample_dir")) real_sample_dir = "./";  # default

#--------- Data Stage ---------#

data_source_fname = "../../../2019 Fall/Medcine Reaction New/std.meds.2y.rdat"
#std.basic.data.list, std.interv.data.list, std.meds.data.list, 
#std.outcome.data.list, std.latest.outcome.data.list, std.twoends.outcome.data.list

# edit from source file
real_data_fname = "real_data.rda";
#basic.data.list, interv.data.list, meds.data.list, outcome.data.list

status_fname = "real_status.rda";
pairs_slt_fname = "real_pairs_slt.rda"
pairs_mtc_fname = "real_pairs_mtc.rda";
psm_res_fname = "real_psm_res.rda";
psm_parel_dname = "real_psm_parel_res/";
psm_log_fname = "real_psm_output_log.txt";

glb_paramt_fname = "real_glb_paramt.rda";
# pat.num

psm_eval_fname = "real_psm_eval.rda";
psm_anal_fname = "real_psm_anal.rda";

# export
real_data_dir = "real_datasets/"

# split
real_split_dir = "real_splits/"

#--------- Pred Stage ---------#

# update
real_updated_dir = "real_outputs_updated/"

# output
real_output_dir = "real_outputs/"
real_pvalue_dir = "real_pvalues/"




#=============== synthetic fnames setting =================#

# defined in Main()
if(!exists("syn_root_dir")) syn_root_dir = "./";  # default
if(!exists("syn_sample_dir")) syn_sample_dir = "./";  # default

#--------- Data Stage ---------#

gen_struct_fname = paste(syn_root_dir,"syn_generation_struct.rda",sep="")
#labs.list, combd.major.list, combd.minor.list, rout.list, med.list

gen_patient_fname = paste(syn_root_dir,syn_sample_dir,"syn_generation_pat.rda",sep="")
#value.init.list, pat.time.obs.mat, med.pat.effect.mat, value.disease.list, pat.rout.list,
#value.treated.list, pat.dx.mat, pat.med.mat, treatment.list, pat.labs.obs.mat

ideal_delta_fname = paste(syn_root_dir,syn_sample_dir,"syn_ideal_delta.rda",sep="")
#delta.trj.list: outcomes.num length list, p.num * year.num-1 matrix delta value
#delta.outcome.list: year.num-1 length list, p.num * outcomes.num matrix delta value

syn_data_fname = paste(syn_root_dir,syn_sample_dir,"syn_data.rda",sep="");
#basic.data.list, interv.data.list, meds.data.list, outcome.data.list

# export
syn_data_dir = paste(syn_root_dir,syn_sample_dir,"syn_datasets",sep="")
syn_suggPreOc_fname = paste(syn_root_dir,"syn_sugg_preOc.csv",sep="")
#sugg.pred.outcome.mat  matrix(pred.slt.num * syn.num) with header

# split
syn_split_dir = paste(syn_root_dir,syn_sample_dir,"syn_splits",sep="")

#--------- Pred Stage ---------#

# model slt variables for steps
syn_sltVars_dir = paste(syn_root_dir, "syn_slt_vars/", sep="")

# update
syn_updated_dir = paste(syn_root_dir,syn_sample_dir,"syn_outputs_updated",sep="")

# output
syn_output_dir = paste(syn_root_dir,syn_sample_dir,"syn_outputs",sep="")
syn_pvalue_dir = paste(syn_root_dir,syn_sample_dir,"syn_pvalues",sep="")
















