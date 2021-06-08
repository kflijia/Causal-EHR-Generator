source("setting_files.R")
source("setting_paramts_data.R")


#--------------------------------------------------------------------------------
# checking: patient data genaration 


check.disease.forecast <- function(in.wd){
  # Require: [ pat.dx.mat ]
  the.pat.fname = paste0(in.wd,"/syn_generation_pat.rda")
  load(the.pat.fname)
  
  fore.yr.arr = 5:15
  fore.combd.mat.1 = matrix(0, nrow=length(fore.yr.arr), ncol=combd.num)
  colnames(fore.combd.mat.1) = disease.arr
  rownames(fore.combd.mat.1) = fore.yr.arr
  fore.combd.mat.all = fore.combd.mat.1
  
  combd.dx.last.arr = matrix(0, nrow=1, ncol=combd.num)
  colnames(combd.dx.last.arr) = disease.arr
  combd.dx.devlped.arr = combd.dx.last.arr
  
  for (p in 1:pat.num) {
    my.dx.yr.mat = matrix(unlist(pat.dx.mat[p,]), combd.num)
    my.dx.last.arr = apply(my.dx.yr.mat, 1, sum)
    combd.dx.last.arr = combd.dx.last.arr + my.dx.last.arr
    combd.dx.devlped.arr = combd.dx.devlped.arr+as.numeric(my.dx.last.arr>0)
    for (yr in fore.yr.arr) {
      yr.str = as.character(yr)
      my.fore.mat = my.dx.yr.mat[, yr:year.num, drop=F]
      my.case.num.1.vec = apply(my.fore.mat,1,sum)
      fore.combd.mat.1[yr.str, ] = fore.combd.mat.1[yr.str, ] + my.case.num.1.vec
      my.case.num.all = (year.num-yr+1)
      my.case.num.all.vec = rep(my.case.num.all, combd.num)
      fore.combd.mat.all[yr.str, ] = fore.combd.mat.all[yr.str, ] + my.case.num.all.vec
    }
  }
  fore.combd.mat.rate = round(fore.combd.mat.1/fore.combd.mat.all, 3)
  
  combd.dx.arr = (pat.num*year.num - combd.dx.last.arr)/pat.num + 1
  combd.dx.arr[combd.dx.arr>year.num] = 0
  colnames(combd.dx.arr) = colnames(combd.dx.last.arr)
  
  combd.dx.devlped.mat = rbind(combd.dx.devlped.arr, round(combd.dx.devlped.arr/pat.num,3))
  rownames(combd.dx.devlped.mat) = c("Num","Rate")
  
  cat("Positive Case Num:\n")
  print(fore.combd.mat.1)
  cat("\n")
  cat("Positive Case Rate:\n")
  print(fore.combd.mat.rate)
  cat("\n")
  cat("Dx Developed Avg Years:\n")
  print(combd.dx.arr)
  cat("\n")
  cat("Dx Developed Pats:\n")
  print(combd.dx.devlped.mat)
  cat("\n")
}

display.pat <- function(in.pat, in.labs=NA, onlymajor=F){
  data.disp.list = list(value.init.list, value.disease.list, value.treated.list)
  data.disp.name = c("initial","disease","treated")
  data.disp.col = c("blue","red","green")
  data.num = length(data.disp.list) # =3
  
  # ------ grab data ------- #
  
  my.rout.id = pat.rout.list$pat.rout.assigned[in.pat]
  my.rout.start = pat.rout.list$pat.start.assigned[in.pat]
  rout.length = apply(rout.list$rout.mat, 1, function(x){sum(!is.na(x))})
  my.rout.end = rout.length[my.rout.id]
  if(my.rout.id==0){
    cat("No disease.\n")
    return()
  }
  routine.full = rout.list$rout.mat[my.rout.id, ]
  my.rout.combd = routine.full[my.rout.start:my.rout.end]
  
  my.major.labs = unlist(lapply(combd.major.list[my.rout.combd], function(x){x[["labs"]]}))
  my.major.crit = unlist(lapply(combd.major.list[my.rout.combd], function(x){x[["crit"]]}))
  my.minor.labs = unlist(lapply(combd.minor.list[my.rout.combd], function(x){paste(x$labs,collapse=",")}))
  
  my.rout.combd
  my.dx.mat.full = t(matrix(unlist(pat.dx.mat[in.pat,]),combd.num)) # year.num*combd.num
  my.dx.mat = my.dx.mat.full[,my.rout.combd, drop=F]
  my.dx.years = apply(my.dx.mat==0,2,sum)+1
  
  pre.rout.combd = c()
  if(my.rout.start>1) {
    pre.rout.combd = routine.full[1:(my.rout.start-1)]
    pre.major.labs = unlist(lapply(combd.major.list[pre.rout.combd], function(x){x[["labs"]]}))
    pre.minor.labs = unlist(lapply(combd.minor.list[pre.rout.combd], function(x){paste(x$labs,collapse=",")}))
    
    cat("combd rout: (",paste(pre.rout.combd,collape="->"),")",paste(my.rout.combd,collapse="->"),"\n")
    cat("major labs: (",paste(pre.major.labs,collape="->"),")",paste(my.major.labs,collapse="->"),"\n")
    cat("minor labs: (",paste("[",pre.minor.labs,"]",sep="",collape="->"),")",paste("[",my.minor.labs,"]", sep="", collapse = "->"),"\n")
  } else {
    cat("combd rout:",paste(my.rout.combd,collapse="->"),"\n")
    cat("major labs:",paste(my.major.labs,collapse="->"),"\n")
    cat("minor labs:",paste("[",my.minor.labs,"]", sep="", collapse = "->"),"\n")
  }
  cat("combd dx years:",paste(my.dx.years,collapse=", "),"\n")
  
  if(is.na(in.labs)){
    in.labs = my.major.labs[1]
  }
  
  cat("-------- develop --------\n")
  cat("Labs =",in.labs,"\n")
  my.major.position = which(my.major.labs==in.labs)
  my.minor.positions = numeric(0);
  for (c in my.rout.combd) {
    if(in.labs %in% combd.minor.list[[c]]$labs){
      my.minor.positions = c(my.minor.positions, which(my.rout.combd==c))
    }
  }
  my.harmful.positions = numeric(0);
  for (c in my.rout.combd) {
    if(in.labs %in% unlist(lapply(medicine.mat[c,], function(x){x$harmful.labs}))){
      my.harmful.positions = c(my.harmful.positions, which(my.rout.combd==c))
    }
  }
  
  position.exist.tag = rep(F, 3)
  names(position.exist.tag) = c("major","minor","harmful");
  
  my.crit = NA
  if(length(my.major.position)>0) {
    position.exist.tag[["major"]] = T
    my.crit = my.major.crit[my.major.position]
    cat("Major position = ",my.major.position,"/",length(my.rout.combd),"\n",sep="")
  } else {
    cat("Major no position\n")
  }
  if(length(my.minor.positions)>0) {
    position.exist.tag[["minor"]] = T
    for(p in 1:length(my.minor.positions)){
      cat("Minor position No.",p," = ",my.minor.positions[p],"/",length(my.rout.combd),"\n",sep="")
    }
  } else {
    cat("Minor no position\n")
  }
  if(length(my.harmful.positions)>0){
    position.exist.tag[["harmful"]] = T
    for(p in 1:length(my.harmful.positions)){
      cat("Harmful position No.",p," = ",my.harmful.positions[p],"/",length(my.rout.combd),"\n",sep="")
    }
  } else {
    cat("Harmful no position\n")
  }
  
  # ------ display data ------ #
  
  disp.data = c()
  for(i in 1:data.num){
    cur.data = data.disp.list[[i]]
    pat.data = cur.data[[in.labs]][in.pat,]
    disp.data = rbind(disp.data, pat.data)
  }
  disp.start = which(apply(disp.data, 2, function(x){length(unique(x))>1}))[1]
  if(is.na(disp.start)){ # no changing, no start
    disp.start = 1
  }
  disp.cols = disp.start:ncol(disp.data)
  
  disp.data = disp.data[,disp.cols]
  disp.obs = pat.time.obs.mat[in.pat,disp.cols] # T or F
  
  y.edge = (max(disp.data)-min(disp.data))*0.2
  y.top = max(disp.data)+y.edge
  y.bot = min(disp.data)-y.edge
  pch.obs.arr = c(16, 17, 18)
  pch.unobs.arr = c(1, 2, 5)
  treated.names = c("1st Line", "2nd Line", "3rd Line", "4th Line", "5th Line")
  # Create a first line
  my.pch = rep(pch.unobs.arr[1], length(disp.cols))
  my.pch[disp.obs] = pch.obs.arr[1]
  plot(disp.cols, disp.data[1,], type = "b", frame = FALSE, pch = my.pch, 
       col = data.disp.col[1], xlab = "Time Windows", ylab = "Value", ylim = c(y.bot,y.top), 
       main=paste("Value Changing of Lab No.",in.labs,"for Patient No.",in.pat))
  for(i in 2:data.num){
    my.pch = rep(pch.unobs.arr[i], length(disp.cols))
    my.pch[disp.obs] = pch.obs.arr[i]
    lines(disp.cols, disp.data[i,], pch = my.pch, col = data.disp.col[i], type = "b", lty=1)
  }
  
  # -- crit -- #
  if(position.exist.tag[["major"]]){  # !is.na(my.crit)
    lines(c(min(disp.cols), max(disp.cols)), c(my.crit,my.crit), type="l", lty=2, col="black")
  }
  
  # -- legend -- #
  legend.table = data.frame(matrix(0,3,length(data.disp.name)+4, dimnames=list(c("col","lty","tag"), c(data.disp.name,"Major","Minor","Harmful","Criterion"))))
  legend.table["col",] = c(data.disp.col,"magenta","orange","deepskyblue1",'black')
  legend.table["lty",] = c(rep(1,data.num),2,2,2,2)
  legend.table["tag",] = TRUE
  
  if(is.na(my.crit)) legend.table["tag","Criterion"] = F;
  if(onlymajor) {
    legend.table["tag",c("Minor","Harmful")] = F;
    colnames(legend.table)[data.num+1] <- "Treatment   ";
  }
  legend.table = legend.table[c("col","lty"), as.logical(legend.table["tag",]), drop=F]
  
  
  legend("topright", legend=colnames(legend.table), col=as.character(legend.table["col",]),
         lty=as.numeric(legend.table["lty",]), cex=1)
  
  
  cat("-------- medicine --------\n")
  # -- major -- #
  if(position.exist.tag[["major"]]){
    my.major.combd = my.rout.combd[my.major.position]
    my.major.treatment = unlist(lapply(treatment.list, function(x){x[in.pat, my.major.combd]})) # len=year.num
    my.major.treated = which(my.major.treatment>0)
    my.major.treated.pos = my.major.treated * 2
    if(length(my.major.treated.pos)>0){
      my.line = my.major.treatment[my.major.treated]
      my.special = unlist(lapply(med.pat.effect.mat[my.major.combd, my.line], function(x){x[in.pat]<=0.2} ))
      # combd.num*max.line, per len=pat.num
      cat("Major med.line=",paste(my.line,collapse=","),sep="")
      cat(" position=",paste(my.major.treated.pos,collapse=","),"/",year.num*2,sep="")
      cat(" Special: ", paste(my.special,collapse=","),"\n",sep="")
      med.idx = 1
      for(j in my.major.treated.pos){
        lines(c(j,j), c(y.bot,y.top), type="l", lty=2, col="magenta",lwd=2)
        text((j+3.5),(y.bot+(y.top-y.bot)/10),labels = treated.names[med.idx],col="magenta")
        med.idx = med.idx+1
      }
    } else {
      cat("Major no medicine\n",sep="")
    }
  }
  # -- minor -- #
  if(position.exist.tag[["minor"]] && !onlymajor){
    my.minor.treated = c()
    for (m in my.minor.positions) {
      cur.combd = my.rout.combd[m]
      cur.treatment = unlist(lapply(treatment.list, function(x){x[in.pat, cur.combd]})) # len=year.num
      cur.treated = which(cur.treatment>0)
      # for minor, each med line contains the same minor
      if(length(cur.treated)>0){
        my.line = cur.treatment[cur.treated]
        cat("Minor No.",match(m,my.minor.positions)," med.line=",paste(my.line,collapse=","),sep="")
        cat(" position=",paste(cur.treated*2,collapse=","),"/",year.num*2,"\n",sep="")
        my.minor.treated = unique(c(my.minor.treated, cur.treated))
      } else {
        cat("Minor No.",match(m,my.minor.positions)," no medicine\n",sep="")
      }
    }
    if(length(my.minor.treated)>0){
      my.minor.treated.pos = my.minor.treated * 2
      for(m in my.minor.treated.pos){
        lines(c(m,m), c(y.bot,y.top), type="l", lty=2, col="orange",lwd=2)
      }
    }
  }
  # -- harmful -- #
  if(position.exist.tag[["harmful"]] && !onlymajor){
    my.harmful.treated = c()
    for (h in my.harmful.positions) {
      cur.combd = my.rout.combd[h]
      cur.treatment = unlist(lapply(treatment.list, function(x){x[in.pat, cur.combd]})) # len=year.num
      cur.treated = which(cur.treatment>0)
      my.line = cur.treatment[cur.treated]
      # for harmful, each med line contains different harmful
      target.line = which(unlist(lapply(medicine.mat[cur.combd,], function(x){in.labs %in% x$harmful.labs})))
      cur.treated = cur.treated[which(my.line %in% target.line)]
      my.line = intersect(my.line, target.line)
      
      if(length(cur.treated)>0){
        cat("Harmful No.",match(h,my.harmful.positions)," med.line=",paste(my.line,collapse=","),sep="")
        cat(" position=",paste(cur.treated*2,collapse=","),"/",year.num*2,"\n",sep="")
        my.harmful.treated = unique(c(my.harmful.treated, cur.treated))
      }
      else {
        cat("Harmful No.",match(h,my.harmful.positions)," no medicine\n",sep="")
      }
    }
    if(length(my.harmful.treated)>0){
      my.harmful.treated.pos = my.harmful.treated * 2
      for(h in my.harmful.treated.pos){
        lines(c(h,h), c(y.bot,y.top), type="l", lty=2, col="deepskyblue1",lwd=2)
      }
    }
  }
  
}


print.pat.dx = function(in.pat, in.combd.arr){
  out.dx.mat = t(matrix(
    unlist(lapply(pat.dx.mat[in.pat,], function(x){x[in.combd.arr]})), nrow=length(in.combd.arr)
  ))
  colnames(out.dx.mat) = paste("combd", in.combd.arr, sep="")
  rownames(out.dx.mat) = paste("Year", 1:year.num, sep=" ")
  print(out.dx.mat)
}

print.pat.devlp = function(in.pat){
  rout.idx = pat.rout.assigned[in.pat]
  my.rout = rout.mat[rout.idx,]
  my.rout = my.rout[!is.na(my.rout)]
  my.major.labs = unlist(lapply(combd.major.list[my.rout], function(x){x[["labs"]]}))
  my.rout.devlp = rout.major.devlp.mat[rout.idx,]
  my.rout.devlp = my.rout.devlp[!is.na(my.rout.devlp)]
  my.major.effect.arr = final.major.effect.list[[in.pat]]
  
  cat("combd rout:", my.rout, "\n")
  cat("major labs:", my.major.labs, "\n")
  cat("rout devlp:", my.rout.devlp, "\n")
  cat("effect arr:", my.major.effect.arr,"\n")
}


#--------------------------------------------------------------------------------
# checking: delta distribution 

check.delta.rate <- function(data.type, delta.type, task.type, the.y){
  data.fname = paste(data.type,"_",delta.type,"_",task.type,"_",the.y,".data.csv",sep="")
  data.fname = paste(data_dir,"/",data.fname,sep="")
  data.flagstr = paste(task.type,"_",the.y,sep="")
  
  if(!file.exists(data.fname) || file.info(data.fname)$size<10){
    cat("\t",data.flagstr,": No Data.\n",sep="")
    return()
  }
  the.data = read.table(data.fname, sep=",", header=T)
  case.num = nrow(the.data)
  
  col.names = colnames(the.data)
  delta.tags = (substr(col.names, 1, 6)=="delta_")
  delta.mat = the.data[,delta.tags, drop=F]
  
  delta.rate = sum(delta.mat!=0)/(case.num*sum(delta.tags))
  delta.rate.str = paste(round(delta.rate*100, 2),"%",sep="")
  cat("\t",data.flagstr,": ",delta.rate.str," \n",sep="")
  
}

check.rout.distr <- function(data.type, delta.type, task.type, the.y){
  # require: rout.list, pat.rout.list
  # pat.rout.list: pat.rout.assigned, pat.start.assigned, pat.start.offset, pat.sevrt.rate
  
  data.fname = paste(data.type,"_",delta.type,"_",task.type,"_",the.y,".data.csv",sep="")
  data.fname = paste(data_dir,"/",data.fname,sep="")
  data.flagstr = paste(task.type,"_",the.y,sep="")
  
  if(!file.exists(data.fname) || file.info(data.fname)$size<10){
    cat(data.flagstr,": No Data.\n\n",sep="")
    return()
  }
  the.data = read.table(data.fname, sep=",", header=T)
  the.colnames = colnames(the.data)
  the.delta.cols = (substr(the.colnames,1,6)=="delta_")
  the.med.cols = (substr(the.colnames,1,9)=="med_combd")
  
  Rout.Mat = rout.list$rout.mat
  Pat.Rout.Assigned = pat.rout.list$pat.rout.assigned
  Pat.Start.Assigned = pat.rout.list$pat.start.assigned
  
  pid.vec = the.data[,"pid"]
  pid.arr = unique(pid.vec)
  pat.rout.vec = Pat.Rout.Assigned[pid.vec]
  pat.rout.arr = Pat.Rout.Assigned[pid.arr]
  pat.start.arr = Pat.Start.Assigned[pid.arr]
  
  pid.num = length(pid.arr)
  pat.med.mat = matrix(0,pid.num,sum(the.med.cols))
  pat.delta.mat = matrix(0, pid.num, n_len)
  for (p in 1:pid.num) {
    cur.med.blc = the.data[pid.vec==(pid.arr[p]), the.med.cols]
    cur.med = as.numeric(apply(cur.med.blc, 2, sum)>0)
    pat.med.mat[p, ] = cur.med
    
    cur.delta.blc = the.data[pid.vec==(pid.arr[p]), the.delta.cols]
    cur.delta = apply((cur.delta.blc!=0), 1, sum)
    pat.delta.mat[p, ] = as.numeric(cur.delta>0)
  }
  
  rout.stats = table(pat.rout.arr)
  rout.stats = rout.stats[order(rout.stats, decreasing = T)]
  rout.arr = as.numeric(names(rout.stats))
  rout.arr = setdiff(rout.arr, "0")
  cat(data.flagstr,": \n",sep="")
  for(r in rout.arr){
    my.pat.idxs = (pat.rout.arr==r)
    my.pat.num = sum(my.pat.idxs)
    
    
    my.pat.data.idxs = (pat.rout.vec==r)
    
    my.case.num = sum(my.pat.data.idxs)
    cat("[Rout_",r,"] Pat.num=",my.pat.num,", Case.num=",my.case.num,". \n",sep="")
    
    cat("             Routine:")
    my.rout = Rout.Mat[r,]
    my.rout.len = sum(!is.na(my.rout))
    my.rout = my.rout[1:my.rout.len]
    cat(" [", paste(my.rout, collapse="]->["), "]\n", sep="")
    
    cat("          Starts.num:")
    my.start.counts = rep(0, my.rout.len)
    names(my.start.counts) = as.character(1:my.rout.len)
    my.start.stats = table(pat.start.arr[my.pat.idxs])
    my.start.counts[names(my.start.stats)] = my.start.stats
    my.start.counts = as.numeric(my.start.counts)
    counts.1 = my.start.counts[nchar(as.character(my.rout))==1]
    counts.2 = my.start.counts[nchar(as.character(my.rout))==2]
    counts.str1 = paste0(" (", paste(counts.1, collapse=")->("), ")->(")
    counts.str2 = paste0(paste(counts.2, collapse=") ->("), ")")
    cat(counts.str1, counts.str2, "\n", sep="")
    
    cat("       Medicined.num:")
    my.med.pnums = apply(pat.med.mat[my.pat.idxs,,drop=F], 2, sum)
    my.rout.med.pnums = my.med.pnums[my.rout]
    cat(" [", paste(my.rout.med.pnums, collapse="]->["), "]\n", sep="")
    
    cat("           Delta.num:")
    yr.seg.len = 5
    my.delta.pnums = apply(pat.delta.mat[my.pat.idxs,,drop=F], 2, sum)
    my.delta.pnum.mat = matrix(0, yr.seg.len, ceiling(length(my.delta.pnums)/yr.seg.len))
    my.delta.pnum.mat[1:length(my.delta.pnums)] = my.delta.pnums
    my.delta.seg.pnums = apply(my.delta.pnum.mat, 2, sum)
    cat(" [", paste(my.delta.seg.pnums, collapse="]->["), "]\n", sep="")
    
  }
}

#--------------------------------------------------------------------------------


data_dir = "syn_experiments/syn_1/syn_sample_5000/syn_datasets/" # setting!

gen_struct_fname = "syn_experiments/syn_1/syn_generation_struct.rda" # setting!
gen_patient_fname = "syn_experiments/syn_1/syn_sample_5000/syn_generation_pat.rda" # setting!

load(gen_struct_fname)
load(gen_patient_fname)

n_len = year.num - forcast.steps
# max number of records per patient can produce

#cat("Checking Delta Rates:\n")
#for(oc in disease.pred.arr) {
#  check.delta.rate(data.type="Dt1", delta.type="Adel", task.type="class", the.y=oc)
#}
#for(oc in outcome.pred.arr) {
#  check.delta.rate(data.type="Dt1", delta.type="Sdel", task.type="regr", the.y=oc)
#}

#check.rout.distr(data.type="Dt1", delta.type="Adel", task.type="class", the.y="combd15")
#data.type="Dt1"; delta.type="Adel"; task.type="class"; the.y="combd15"




#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# checking: delta information 

# [ FOLDED BY JIA. Apr2020 ]















