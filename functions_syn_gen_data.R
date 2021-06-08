

rand.labs <- function(labs.num){
  labs.mean.arr = round(runif(labs.num)*145) + 5
  labs.sd.arr = round(runif(labs.num,0.1,0.5),2)#round(runif(labs.num)*0.5+0.5, 2)
  labs.list=list(labs.mean.arr=labs.mean.arr, labs.sd.arr=labs.sd.arr)
  cat("    rand.labs() done.\n")
  return(labs.list)
}

init.value <- function(){
# Require: [pat.num, labs.list, year.num, wide.prop, obs.prob]
  labs.mean.arr = labs.list$labs.mean.arr
  labs.sd.arr = labs.list$labs.sd.arr
  
  value.init.list = vector("list",labs.num)
  for(i in 1:labs.num){
    value.init.list[[i]] = matrix(NA, pat.num, year.num*2) # each year has 2 measures, bf + af
    my.mean = labs.mean.arr[i]
    my.sd = labs.sd.arr[i]
    my.wide = my.mean/wide.prop
    
    # random value pool
    my.value.pool = rnorm(round(pat.num*year.num*2*1.2), mean = 0, sd = my.sd) 
    my.value.bottom = quantile(my.value.pool, 0.05)
    my.value.upper = quantile(my.value.pool, 0.95)
    my.value.pool = my.value.pool[my.value.pool>=my.value.bottom & my.value.pool<=my.value.upper]
    
    # select values
    my.value.vec = my.value.pool[1:(pat.num*year.num*2)]
    my.value.vec = (my.value.vec-min(my.value.vec))/(max(my.value.vec)-min(my.value.vec)) # standardize
    my.value.vec = my.value.vec * my.wide * 2 + my.mean - my.wide # move to mean center
    value.init.list[[i]][,] = my.value.vec
  }
  
  # randomize observable years
  pat.time.obs.mat = matrix(FALSE, nrow=pat.num, ncol=year.num*2)
  obs.sample.num = ceiling(length(pat.time.obs.mat)*obs.prob)
  obs.samples = sample(length(pat.time.obs.mat), obs.sample.num)
  pat.time.obs.mat[obs.samples] = TRUE
  # pat.time.obs.mat will be changed in get.treatment, the year has med will have T on bf
  
  cat("    init.value() done.\n")
  return(list(value.init.list, pat.time.obs.mat))
}

library("purrr")
library("bnlearn")
rand.routine <- function(){
# Require: [combd.num, connect.rate]
  if(combd.num>20){
    cat("    combd.num must <= 20\n")  # control complexity
    return()
  }
  if(!is.numeric(connect.rate) || connect.rate<=0 || connect.rate>=1){
    cat("    connect.rate must in (0,1)\n")
    return()
  }
  
  # want a long graph
  glength.min = ceiling(combd.num*0.5)
  gwidth.max = floor(combd.num*0.4)
  
  #---------------------------------------------
  
  merge.graph <- function(graph.a, graph.b){
    if(length(graph.a)==0) {
      nodes.num.b = length(graph.b$nodes)
      graph.ab = graph.b[c("nodes","arcs")]
      graph.ab$covers = c(nrow(graph.b$arcs), nodes.num.b*(nodes.num.b-1)/2)
      return(graph.ab)
    }
    
    nodes.a = names(graph.a$nodes); nodes.num.a = length(nodes.a)
    nodes.b = names(graph.b$nodes); nodes.num.b = length(nodes.b)
    nodes.ovlp = intersect(nodes.a, nodes.b)
    nodes.new = setdiff(nodes.b, nodes.a)
    
    merge.nodes = graph.a$nodes
    for(n in nodes.ovlp){
      merge.nodes[[n]] = Map(c, merge.nodes[[n]], graph.b$nodes[[n]])
      merge.nodes[[n]]$mb = unique(merge.nodes[[n]]$mb)
      merge.nodes[[n]]$nbr = unique(merge.nodes[[n]]$nbr)
      merge.nodes[[n]]$parents = unique(merge.nodes[[n]]$parents)
      merge.nodes[[n]]$children = unique(merge.nodes[[n]]$children)
    }
    merge.nodes = append(merge.nodes, graph.b$nodes[nodes.new])
    merge.arcs = unique(rbind(graph.a$arcs, graph.b$arcs))
    
    if(length(graph.a$covers)>0) {
      covers.a = graph.a$covers
    } else {
      covers.a = c(nrow(graph.a$arcs), nodes.num.a*(nodes.num.a-1)/2)
    }
    if(length(graph.b$covers)>0) {
      covers.b = graph.b$covers
    } else {
      covers.b = c(nrow(graph.b$arcs), nodes.num.b*(nodes.num.b-1)/2)
    }
    merge.covers = covers.a + covers.b
    names(merge.covers) = c("revealed","fullcover")

    graph.ab = list(nodes=merge.nodes, arcs=merge.arcs, covers=merge.covers)
    return(graph.ab)
  }
  
  random.long.graph <- function(in.nodes, in.num, in.method, in.prob){
    part.size = 5 # num of nodes in each part
    ovlp.size = 2 # num of nodes overlaped between parts
    
    left.nodes = in.nodes
    the.graph.list = vector("list",in.num)
    while(length(left.nodes)>0){
      left.size = length(left.nodes)
      if(left.size<=part.size){
        cur.nodes = left.nodes
        left.nodes = c()
      } else {
        cur.nodes = left.nodes[1:part.size]
        left.nodes = left.nodes[(part.size-ovlp.size+1):left.size]
      }
      cur.graph.list = random.graph(nodes=cur.nodes, num=in.num, method=in.method, prob=in.prob);
      # learning, nodes, arcs
      for (i in 1:in.num) {
        the.graph.list[[i]] = merge.graph(the.graph.list[[i]], cur.graph.list[[i]])
        # nodes, arcs, covers(revealed, fullcover)
      }
    }
    
    the.graph.list
  }
  
  
  # get graph length/width, entrance/exit width, connect rate
  # maybe: max in/out.degree
  get.shape <- function(in.graph){
    nodes.list = in.graph$nodes;
    nodes.arr = names(nodes.list);
    nodes.num = length(nodes.arr);
    entrance.tag = unlist(lapply(nodes.list, function(x){length(x$parents)==0}))
    entrance.arr = nodes.arr[entrance.tag]
    exit.tag = unlist(lapply(nodes.list, function(x){length(x$children)==0}))
    exit.arr = nodes.arr[exit.tag]
    
    broad.level.arr = data.frame(matrix(0,nrow=1,ncol=nodes.num,dimnames = list(NULL, nodes.arr)))
    depth.level.arr = broad.level.arr
    # leave entrances to be 0
    cur.level.nodes = entrance.arr # using depth-first
    cur.level = 0;
    while(length(cur.level.nodes)>0){
      next.level.nodes = unique(unlist(lapply(nodes.list[cur.level.nodes], function(x){x$children})))
      depth.level.arr[next.level.nodes] <- (cur.level+1)
      cur.level.nodes <- next.level.nodes
      
      rm.tag = broad.level.arr[next.level.nodes] > 0
      next.level.nodes = next.level.nodes[!rm.tag]
      broad.level.arr[next.level.nodes] <- (cur.level+1)
      cur.level <- cur.level + 1
    }
    
    graph.length = max(depth.level.arr);
    graph.width =  max(table(as.numeric(broad.level.arr)));
    
    covers.arr = in.graph$covers
    covers.revealed = covers.arr[1]
    covers.fullcover = covers.arr[2]
    actual.connect.rate = as.numeric(round(covers.revealed/covers.fullcover, 3))
    
    # check connected
    connected.tag.arr = data.frame(matrix(F,nrow=1,ncol=nodes.num,dimnames = list(NULL, nodes.arr)))
    blkt.border = entrance.arr[1] # blanket border, start from one entrance
    while(length(blkt.border)>0){
      connected.tag.arr[blkt.border] = T
      blkt.border = unique(unlist(lapply(nodes.list[blkt.border], function(x){x$mb})))
      rm.tag = connected.tag.arr[blkt.border]
      blkt.border = blkt.border[!rm.tag]
    }
    connected.tag = all(as.matrix(connected.tag.arr))
    
    return(list(glength = graph.length, gwidth = graph.width, nodes.depth = depth.level.arr,
                entr.arr = entrance.arr, exit.arr = exit.arr, 
                cnt.tag = connected.tag, cnt.rate = actual.connect.rate))
  }
  
  eval.shape <- function(in.shape){
    quali.val = 0
    if(!in.shape$cnt.tag) return(quali.val)
    if(in.shape$cnt.rate < connect.rate) return(quali.val)
    if(in.shape$cnt.rate > (connect.rate+0.1)) return(quali.val)
    
    shape.ratio = round(in.shape$glength/in.shape$gwidth, 2)
    quali.tag = (in.shape$glength >= glength.min) && (in.shape$gwidth <= gwidth.max)
    if(quali.tag) quali.val = shape.ratio
    return(quali.val)
  }
  
  get.rout<- function(in.graph, in.shape){
    nodes.list = in.graph$nodes
    nodes.arr = names(nodes.list);
    nodes.num = length(nodes.arr);
    arcs.mat = in.graph$arcs
    entr.arr = in.shape$entr.arr
    exit.arr = in.shape$exit.arr
    exit.num = length(exit.arr)
    nodes.depth = in.shape$nodes.depth
    glength = in.shape$glength
    
    out.graph = in.graph[c("nodes", "arcs")]
    if(exit.num>1){
      # alter to sigle exit
      exit = exit.arr[exit.num]
      exit.depth = nodes.depth[[exit]]
      for(e in exit.arr[1:(exit.num-1)]){
        arcs.mat = rbind(arcs.mat, c(e,exit));
        nodes.list[[e]]$mb = unique(c(nodes.list[[e]]$mb, exit))
        nodes.list[[e]]$nbr = unique(c(nodes.list[[e]]$nbr, exit))
        nodes.list[[e]]$children = unique(c(nodes.list[[e]]$children, exit))
        if(exit.depth < nodes.depth[[e]]+1) exit.depth=nodes.depth[[e]]+1;
      }
      out.graph$nodes = nodes.list
      out.graph$arcs = arcs.mat
      nodes.depth[[exit]] = exit.depth
      glength = max(nodes.depth)
    }
    out.graph$nodes.depth = nodes.depth
    out.graph$glength = glength
    # nodes, arcs, nodes.depth, glength
    
    connect.mat = matrix(NA,nodes.num,nodes.num)
    colnames(connect.mat) = nodes.arr
    rownames(connect.mat) = nodes.arr
    connect.mat[upper.tri(connect.mat)] = 0
    for(r in 1:nrow(arcs.mat)){
      connect.mat[arcs.mat[r,"from"], arcs.mat[r,"to"]] = 1
    }
    
    rout.len = glength+1
    rout.mat = matrix(NA,nrow=0,ncol=rout.len)
    next.step <- function(prev.rout, cur.node, rout.mat){
      next.nodes = nodes.list[[cur.node]]$children
      if(length(next.nodes)==0){
        out.rout = c(prev.rout, cur.node)
        out.rout = match(out.rout, nodes.arr)
        out.rout = c(out.rout, rep(0,rout.len-length(out.rout)))
        rout.mat = rbind(rout.mat, out.rout)
      } else {
        for(n in next.nodes){
          rout.mat = next.step(c(prev.rout, cur.node), n, rout.mat)
        }
      }
      return(rout.mat)
    }
    for(e in entr.arr){
      rout.mat = next.step(c(),e,rout.mat)
    }
    rownames(rout.mat) <- NULL
    return(list(out.graph=out.graph, rout.mat=rout.mat, connect.mat=connect.mat))
  }
  
  #------------------------
  
  try.batch.size = 10; # number of graph per try
  
  the.graph = list()
  the.shape = list()
  
  while(length(the.graph)==0){
    candi.list = random.long.graph(in.nodes = LETTERS[1:combd.num], in.num = try.batch.size, 
                                   in.method = "ordered", in.prob=connect.rate);
    candi.shape.list = lapply(candi.list, get.shape)
    candi.quali.arr = unlist(lapply(candi.shape.list, eval.shape))
    if(any(candi.quali.arr>0)){
      slt.candi.idx = which.max(candi.quali.arr)
      the.graph = candi.list[[slt.candi.idx]]
      # nodes, arcs
      the.shape = candi.shape.list[[slt.candi.idx]]
      # glength, gwidth, nodes.depth, entr.arr, exit.arr, cnt.tag, cnt.rate
    }
  }
  
  out.list = get.rout(the.graph, the.shape)
  out.graph = out.list$out.graph # nodes, arcs, nodes.depth, glength
  rout.mat = out.list$rout.mat
  connect.mat = out.list$connect.mat
  rout.mat[rout.mat==0] = NA
  
  rout.num.max = 20
  if(nrow(rout.mat)>=rout.num.max){
    slt.step = (nrow(rout.mat)-1)/(rout.num.max-1)
    slt.idxs = round(seq(1, nrow(rout.mat), slt.step))
    rout.lens = apply(!is.na(rout.mat),1,sum)
    max.steps = max(rout.lens[slt.idxs])
    rout.mat = rout.mat[slt.idxs,1:max.steps,drop=F]
  }
  
  cat("    Num of Routines =", nrow(rout.mat))
  cat(", Max Step =",out.graph$glength,"\n")
  rout.list = list(out.graph=out.graph, rout.mat=rout.mat, connect.mat=connect.mat)
  
  cat("    rand.routine() done.\n")
  return(rout.list)
}

rand.combd <- function(){
# Require: [labs.list, rout.list, wide.prop]
  labs.mean.arr = labs.list$labs.mean.arr
  labs.sd.arr = labs.list$labs.sd.arr
  labs.num = length(labs.mean.arr)
  rout.mat = rout.list$rout.mat
  rout.num = nrow(rout.mat)
  max.step = ncol(rout.mat)
  connect.mat = rout.list$connect.mat
  combd.num = ncol(connect.mat)
  if(labs.num<=combd.num){
    cat("    labs.num must > combd.num\n")
    return()
  }
  
  # fix critirian for each lab
  labs.crit.arr = rep(0, labs.num)
  for(i in 1:labs.num){
    my.lab.mean = labs.mean.arr[i]
    my.lab.sd = labs.sd.arr[i]
    
    my.sample.wide = my.lab.mean/wide.prop
    my.samples = rnorm(10000, mean = 0, sd = my.lab.sd)
    my.bottom = quantile(my.samples, 0.05)
    my.upper = quantile(my.samples, 0.95)
    
    my.samples = (my.samples-my.bottom)/(my.upper-my.bottom) # standardize
    my.samples = my.samples * my.sample.wide * 2 + my.lab.mean - my.sample.wide
    my.crit = round(max(my.samples))
    labs.crit.arr[i] = my.crit
  }
  labs.list$labs.crit.arr=labs.crit.arr
  
  # assign major impact
  combd.major.labs = sample(1:labs.num, combd.num, F)
  combd.major.list = vector("list", combd.num)
  for(i in 1:combd.num){
    my.major = combd.major.labs[i]
    my.value = c(my.major, labs.mean.arr[my.major], labs.crit.arr[my.major])
    names(my.value) = c("labs","mean","crit")
    combd.major.list[[i]] = my.value
  }
  # assign minor impact
  other.labs = 1:labs.num
  other.labs = other.labs[-combd.major.labs]
  other.labs.assigned.combd = sample(combd.num, length(other.labs), replace = T)
  combd.minor.list = vector("list", combd.num)
  for (i in 1:combd.num) {
    # all major labs just after i in combd routs belong to minors
    cuz.idxs = which(t(rout.mat)==i)
    cuz.idxs = cuz.idxs[cuz.idxs %% max.step > 0] # rm the last step
    my.minor.arr = t(rout.mat)[cuz.idxs+1] # combds
    my.minor.arr = unique(my.minor.arr[!is.na(my.minor.arr)])
    my.minor.arr = combd.major.labs[my.minor.arr] # combd -> labs
    my.minor.others = other.labs[other.labs.assigned.combd==i]
    my.minor.arr = c(my.minor.arr, my.minor.others)
    if(length(my.minor.arr)==0){
      #combd.minor.list[[i]]=list(labs=integer(0), mean=numeric(0), crit=numeric(0))
      #next;
      my.minor.arr = c(sample(combd.major.labs, 1), sample(other.labs,1))
      # guarantee minors exist
    }
    my.mean.arr = labs.mean.arr[my.minor.arr]
    my.crit.arr = labs.crit.arr[my.minor.arr]
    my.value.list = list(labs=my.minor.arr, mean=my.mean.arr, crit=my.crit.arr)
    combd.minor.list[[i]] = my.value.list
  }
  
  cat("    rand.combd() done.\n")
  return(list(labs.list, combd.major.list, combd.minor.list))
  
}


# =========== Disease ============== #

rand.devlp.segs <- function(rand.size, devlp.type="major"){
  if(devlp.type!="major" && devlp.type!="minor"){
    cat("Wrong devlp type: major or minor.\n")
    return()
  }
  
  devlp.rand.var = 0.6 # degree of flat
  devlp.rand.bott = -0.7 # cut off bottom, how bias toward low number of years
  devlp.rand.top = 2 # cut off top, how long the tailis
  devlp.rand.wide = 22 # wide of [min # of years, max # of years]
  devlp.rand.offset = 1 #2 # the lowest possible # of years
  
  element.num=5000
  
  devlp.samples = rnorm(element.num*5,0,devlp.rand.var) # random
  devlp.samples = devlp.samples[devlp.samples>=devlp.rand.bott & devlp.samples<=devlp.rand.top] # cutoff
  devlp.samples = round((devlp.samples-min(devlp.samples))/(max(devlp.samples)-min(devlp.samples))*devlp.rand.wide)+devlp.rand.offset # standersize and move
  # around 1~23, median 8 
  
  shrink.rate = 1.5
  major.devlp.samples = round((devlp.samples-devlp.rand.offset)/shrink.rate + devlp.rand.offset)
  # around 1~16, median 6
  minor.devlp.samples = devlp.samples+10
  # around 11~33, median 18
  
  if(devlp.type=="major"){
    out.samples = sample(major.devlp.samples, rand.size, T)
  } else {
    out.samples = sample(minor.devlp.samples, rand.size, T)
  }
  return(out.samples)
}

sevrt.2.devlp <- function(in.devlp.arr, in.sevrt){
# apply the sevrt rate onto the devlp
# sevrt: positive = more illed, negative = less illed
  if(!is.numeric(in.sevrt) || !(length(in.sevrt)==1 || length(in.sevrt)==length(in.devlp.arr))) {
    cat("Invalid sevrt input:", in.sevrt, "\n")
    return()
  }
  #if(any(!is.numeric(in.devlp.arr) && !is.infinite(in.devlp.arr)) || any(in.devlp.arr<1)) {
  if(any(is.na(in.devlp.arr)) || !is.numeric(in.devlp.arr) || any(in.devlp.arr<1)) {
    cat("Invalid devlp input:", paste(in.devlp.arr,collapse=", "), "\n")
    return()
  }
  
  # in.sevrt may be array
  if(length(in.sevrt)==1) {
    in.sevrt.arr = rep(in.sevrt, length(in.devlp.arr))
  } else {
    in.sevrt.arr = in.sevrt
  }
  
  the.min.arr = round(in.devlp.arr/2); the.min.arr[the.min.arr<1]=1
  the.max.arr = round(in.devlp.arr*2)
  
  # initially assume in.sevrt>0
  the.botm = the.min.arr
  the.up = in.devlp.arr
  neg.tags = (in.sevrt.arr<=0)
  the.botm[neg.tags] = in.devlp.arr[neg.tags]
  the.up[neg.tags] = the.max.arr[neg.tags]

  out.devlp.arr = in.devlp.arr - (the.up - the.botm)*in.sevrt.arr
  out.devlp.arr[is.na(out.devlp.arr)] = Inf # in.devlp.arr may contains Inf leads to NaN
  out.devlp.arr
}

merge.devlp <- function(in.devlp.1.arr, in.devlp.2.arr){
  if(length(in.devlp.1.arr)!=length(in.devlp.2.arr)){
    cat("Different merge.devlp input length: [", paste0(in.devlp.1.arr,collapse=","), "] [",
        paste0(in.devlp.2.arr,collapse=","), "]\n",sep="")
    return()
  }
  in.stack = rbind(in.devlp.1.arr, in.devlp.2.arr);
  
  if(any(is.na(in.stack)) || !is.numeric(in.stack) || any(in.stack<1) ) {
    cat("Invalid merge.devlp input: [", paste0(in.devlp.1.arr,collapse=","), "] [",
        paste0(in.devlp.2.arr,collapse=","), "]\n",sep="")
    return()
  }
  
  devlp.length = length(in.devlp.1.arr)
  new.devlp.arr = rep(Inf, devlp.length)
  for(i in 1:devlp.length){
    cur.pair = in.stack[,i]
    high.devlp = max(cur.pair)
    low.devlp = min(cur.pair)
    if(is.infinite(high.devlp)){
      new.devlp.arr[i] = low.devlp;
      next
    }
    high.devlp.ebase = 2^(1/high.devlp)
    # devlp is on years! so the ebase is on years, not times
    added.height = (high.devlp.ebase^low.devlp)/2
    # the less sever devlp is (halfly=0.5 shrunk) added on more sever one
    new.height = (2+added.height)
    new.ebase = new.height^(1/low.devlp)
    new.devlp = log(2, base = new.ebase)
    new.devlp.arr[i] = new.devlp
  }

  new.devlp.arr
}


rand.extd = 3 # used in devlp.2.labs(), get.disease(), get.treatment()

devlp.2.labs <- function(in.maps, start.yr, end.yr){
# only calculate upper part values of time-series 
# bottom part is rep(my.mean, time.len)
  devlp.map = in.maps$devlp.map # devlp is on years, not times
  mean.map = in.maps$mean.map
  crit.map = in.maps$crit.map
  
  if(length(unique(c(length(devlp.map),length(mean.map),length(crit.map))))>1 
     || length(devlp.map)!=labs.num){
    cat("    devlp map wrong!\n")
    return()
  }
  if(start.yr<1 || end.yr>(year.num+longest.rout.years+rand.extd)){
    cat("    devlp.2.labs yr wrong!\n")
    return()
  }
  
  first.yr = 1
  offset.steps = (start.yr-first.yr)*2
  time.seq = c((first.yr-0.5), seq(first.yr,end.yr, 0.5))
  the.values.mat=matrix(0, nrow=labs.num, ncol=length(time.seq)) # need to be cut offset.steps
  for (l in 1:labs.num) {
    my.devlp = devlp.map[l]
    if(is.infinite(my.devlp)) next;
    
    my.mean = mean.map[l]
    my.crit = crit.map[l]
    my.ebase = 2^(1/my.devlp)
    my.value.seris = ((my.ebase^time.seq)-1)*(my.crit-my.mean) # 1~2 => 0~1
    
    # turn point
    above.time.idx = which((my.value.seris+my.mean-my.crit)>0) # if observed longer than devlp
    if(length(above.time.idx)>0){
      above.time.seq = time.seq[above.time.idx]
      turn.time.point = min(above.time.seq)-0.5
      turn.point = match(turn.time.point, time.seq)
      turn.base.value = my.value.seris[turn.point]
      turned.time.seq = turn.time.point - (above.time.seq - turn.time.point) # decreasing
      turned.value.seris = ((my.ebase^turned.time.seq)-1)*(my.crit-my.mean) # decreasing
      #turned.value.seris = (turn.base.value - turned.value.seris)*0.5 + turn.base.value # increasing and shrink half
      turned.value.seris = (turn.base.value - turned.value.seris) + turn.base.value # increasing and shrink half
      my.value.seris[above.time.idx] = turned.value.seris
    }
    the.values.mat[l,] = my.value.seris
  }
  
  out.values.mat = the.values.mat[,(offset.steps+1):length(time.seq),drop=F]
  
  return(out.values.mat)
}

#--------------------------------------------------------------------

# lengend: step=routine stage, year=real year, time=time stage (year*2)

rand.severity <- function(){
# Require: [rout.list, year.num, combd.minor.list]
  
  # randomize severity along routines
  rout.mat=rout.list$rout.mat
  rout.num=nrow(rout.mat)
  max.step=ncol(rout.mat)
  element.num = rout.num*max.step
  rout.major.devlp.mat=matrix(0, rout.num, max.step) # initial as 0 instead of Inf
  rout.minor.devlp.mat=matrix(vector("list",element.num), rout.num, max.step)
  
  # for rout.major.devlp.mat
  arcs.mat = cbind(as.numeric(rout.mat[,1:(max.step-1)]), as.numeric(rout.mat[,2:max.step]))
  arcs.mat = unique(arcs.mat)
  enter.nodes = unique(rout.mat[,1])
  arcs.mat = rbind(arcs.mat, cbind(rep(0,length(enter.nodes)), enter.nodes))
  arcs.mat = arcs.mat[apply(is.na(arcs.mat), 1, sum)==0,,drop=F]
  arcs.identities = apply(arcs.mat, 1, function(x){paste0(x,collapse="_")})
  arcs.mat = arcs.mat[order(arcs.identities),]
  arcs.identities = sort(arcs.identities)
  
  arcs.num = nrow(arcs.mat)
  arcs.weights = rand.devlp.segs(arcs.num, "major")
  arcs.mat = cbind(arcs.mat,arcs.weights)
  colnames(arcs.mat) = c("start","end", "weight")
  
  # start from 1st step
  my.identities = apply(cbind(rep(0,rout.num),rout.mat[,1]), 
                        1, function(x){paste0(x,collapse="_")})
  my.weights = arcs.mat[match(my.identities, arcs.identities), "weight"]
  rout.major.devlp.mat[,1] = my.weights
  for (i in 2:max.step) {
    my.arcs = rout.mat[,(i-1):i,drop=F]
    keep.idx = !apply(is.na(my.arcs),1,any)
    my.arcs = my.arcs[keep.idx,,drop=F]
    my.identities = apply(my.arcs, 1, function(x){paste0(x,collapse="_")})
    my.weights = arcs.mat[match(my.identities, arcs.identities),"weight"]
    rout.major.devlp.mat[keep.idx,i] = my.weights
  }
  
  # for rout.minor.devlp.mat
  for(i in 1:rout.num){
    for(j in 1:max.step){
      my.combd = rout.mat[i,j]
      if(is.na(my.combd)){
        next;
      }
      my.minors = combd.minor.list[[my.combd]]$labs
      my.minors.devlp = rand.devlp.segs(length(my.minors), "minor")
      rout.minor.devlp.mat[[i,j]] = my.minors.devlp
    }
  }
  
  # change devlp.segs to accumulated devlp
  for (i in 1:rout.num){
    for (j in 2:max.step) {
      my.combd = rout.mat[i,j]
      if(is.na(my.combd)){
        next;
      }
      prev.devlp = rout.major.devlp.mat[i,(j-1)]
      rout.major.devlp.mat[i,j] = rout.major.devlp.mat[i,j]+prev.devlp
      rout.minor.devlp.mat[[i,j]] = rout.minor.devlp.mat[[i,j]]+prev.devlp
    }
  }
  
  
  # save
  rout.list$rout.major.devlp.mat = rout.major.devlp.mat
  rout.list$rout.minor.devlp.mat = rout.minor.devlp.mat
  
  cat("    rand.severity() done.\n")
  return(rout.list)
}

rand.pat.severity <- function(){ 
# Require: [rout.list]
  rout.major.devlp.mat = rout.list$rout.major.devlp.mat
  
  # randomize routines
  pat.rout.assigned = sample(0:rout.num, pat.num, T)
  
  # prepare for rand start point
  # rout.major.devlp.mat is accumulated, and initialized by 0
  rout.devlp.ynum.arr = apply(rout.major.devlp.mat, 1, max) 
  rout.devlp.left.mat = matrix(0, nrow(rout.major.devlp.mat), ncol(rout.major.devlp.mat))
  for (m in 1:max.step) {
    rout.devlp.left.mat[,m] = rout.devlp.ynum.arr - rout.major.devlp.mat[,m]
  }
  latest.start.arr = apply(rout.devlp.left.mat>=year.num, 1, sum)
  latest.start.arr[latest.start.arr==0] = 1
  
  # randomize start point 
  pat.start.assigned = rep(0, pat.num)
  for (i in 1:rout.num) {
    my.latest = latest.start.arr[i]
    my.pat.idx = (pat.rout.assigned==i)
    pat.start.assigned[my.pat.idx] = sample(1:my.latest, sum(my.pat.idx), replace = T)
  }

  # randomize offset years
  # to solve the 0 rout.assigned
  rout.major.devlp.mat = rbind(rep(0,max.step),rout.major.devlp.mat) # (rout.num+1) * max.step
  
  pat.devlp.accum.mat = rout.major.devlp.mat[(pat.rout.assigned+1),] # pat.num * max.step
  pat.devlp.accum.mat = cbind(rep(0,pat.num), pat.devlp.accum.mat) # pat.num * (max.step+1)
  pat.start.t.idxs = (1:pat.num-1)*ncol(pat.devlp.accum.mat)+pat.start.assigned+1
  max.values = t(pat.devlp.accum.mat)[pat.start.t.idxs]
  pat.devlp.accum.mat = cbind(rep(0,pat.num), pat.devlp.accum.mat) # pat.num * (max.step+2)
  pat.start.t.idxs = (1:pat.num-1)*ncol(pat.devlp.accum.mat)+pat.start.assigned+1
  min.values = t(pat.devlp.accum.mat)[pat.start.t.idxs]

  diff.values = max.values - min.values
  prop.rand = runif(pat.num,0,1)
  pat.start.offset = min.values+(diff.values*prop.rand)
  
  # randomize rate
  pat.sevrt.rate = rnorm(pat.num, mean=0, sd=0.1)
  pat.rout.list = list(pat.rout.assigned=pat.rout.assigned, pat.start.assigned=pat.start.assigned, 
                       pat.start.offset=pat.start.offset, pat.sevrt.rate=pat.sevrt.rate)
  
  cat("    rand.pat.severity() done.\n")
  return(pat.rout.list)
}


get.disease <- function(){
# Require: [value.init.list, combd.major.list, combd.minor.list, rout.list, pat.rout.list]
  
  #value.init.list: len=labs.num, pat.num * year.num*2
  #combd.major.list: len=combd.num, array names: labs, mean, crit
  #combd.minor.list: len=combd.num, list names: labs, mean, crit
  #rout.list: rout.mat(rout.num*max.step), connect.mat(combd.num+1*combd.num+1), 
  #           rout.major.devlp.mat(rout.num*max.step), rout.minor.devlp.mat(rout.num*max.step)
  #pat.rout.list: pat.rout.assigned, pat.start.assigned, pat.start.offset, pat.sevrt.rate

  rout.major.devlp.mat = rout.list$rout.major.devlp.mat # rout.num * max.step, number
  rout.minor.devlp.mat = rout.list$rout.minor.devlp.mat # rout.num * max.step, array
  
  pat.rout.assigned = pat.rout.list$pat.rout.assigned
  pat.start.offset = pat.rout.list$pat.start.offset
  pat.sevrt.rate = pat.rout.list$pat.sevrt.rate
  
  # start
  value.disease.list = value.init.list
  rout.lens = apply(rout.mat, 1, function(x){sum(!is.na(x))})
  maps.list = list() # len = 3
  for(i in 1:3){
    maps.list[[i]] = matrix(NA, nrow = pat.num, ncol = labs.num)
  }
  names(maps.list) = c("devlp.map", "mean.map", "crit.map")
  
  for(i in 1:pat.num){
    #cat(i," ");
    my.rout.id = pat.rout.assigned[i]
    if(my.rout.id==0){
      next;
    }
    my.rout = rout.mat[my.rout.id,];
    my.rout.offset = pat.start.offset[i] # years offset
    my.rout.end = rout.lens[my.rout.id]
    
    my.pat.rate = pat.sevrt.rate[i]
    # creat devlp map
    my.devlp.map = rep(Inf, labs.num) # base number=Inf, means 0 increase
    my.mean.map = rep(NA, labs.num)
    my.crit.map = rep(NA, labs.num)
    for (j in 1:my.rout.end) { # calculate full rout, then cut off the offset
      the.combd = my.rout[j]
      
      # major
      the.major = combd.major.list[[the.combd]]
      the.major.labs = the.major[["labs"]]
      the.major.mean = the.major[["mean"]]
      the.major.crit = the.major[["crit"]]
      
      my.std.devlp = rout.major.devlp.mat[my.rout.id,j]
      my.new.devlp = sevrt.2.devlp(my.std.devlp, my.pat.rate)
      my.exist.devlp = my.devlp.map[the.major.labs]
      
      my.devlp = merge.devlp(my.exist.devlp, my.new.devlp) # may be Inf
      my.devlp.map[the.major.labs] = my.devlp
      my.mean.map[the.major.labs] = the.major.mean
      my.crit.map[the.major.labs] = the.major.crit
      
      # minors
      the.minor = combd.minor.list[[the.combd]]
      if(length(the.minor)==0) next;
      the.minor.labs = the.minor[["labs"]]
      the.minor.mean = the.minor[["mean"]]
      the.minor.crit = the.minor[["crit"]]
      
      my.std.devlp = rout.minor.devlp.mat[[my.rout.id,j]]
      my.new.devlp = sevrt.2.devlp(my.std.devlp, my.pat.rate)
      my.exist.devlp = my.devlp.map[the.minor.labs]
      
      my.devlp = merge.devlp(my.exist.devlp, my.new.devlp) # may be Inf
      my.devlp.map[the.minor.labs] = my.devlp
      my.mean.map[the.minor.labs] = the.major.mean
      my.crit.map[the.minor.labs] = the.major.crit
    }
    my.maps = list(devlp.map=my.devlp.map, mean.map=my.mean.map, crit.map=my.crit.map)
    
    # then cut offset
    labs.value.mat = devlp.2.labs(in.maps=my.maps, start.yr=1, end.yr=year.num+floor(my.rout.offset)+rand.extd) 
    # labs.num * ((year.num+floor(my.rout.offset)+rand.extd)*2)
    offset.times = floor(my.rout.offset*2)
    labs.value.mat = labs.value.mat[,(offset.times+1):(offset.times+year.num*2),drop=F]

    for (l in 1:labs.num) {
      cur.values = value.disease.list[[l]][i,]
      cur.values = cur.values + labs.value.mat[l,]
      value.disease.list[[l]][i,] = cur.values
    }
    maps.list$devlp.map[i,] = my.devlp.map
    maps.list$mean.map[i,] = my.mean.map
    maps.list$crit.map[i,] = my.crit.map
  } # i = 1:pat.num
  
  #cat("\n")
  cat("    get.disease() done.\n")
  return(list(value.disease.list, maps.list))
}


# ============ Treatment ============= #

rand.medicine <- function(){
# Require: [labs.num, combd.major.list, combd.minor.list, max.line, harmful.param]
  
  # randomize medications for each comobidity
  # max.line: the max number of medicine for each comobidity
  # harmful.param: max.num: max number of labs harmful to. max.prob: max probability of being harmful
  #                min.num: min number of labs harmful to. min.prob: max probability of being harmful
  
  #combd.major.list: len=combd.num, array names: labs, mean, crit
  #combd.minor.list: len=combd.num, list names: labs, mean, crit

  harmful.min.num = harmful.param[["min.num"]]
  harmful.max.num = harmful.param[["max.num"]]
  harmful.min.prob = harmful.param[["min.prob"]]
  harmful.max.prob = harmful.param[["max.prob"]]
  
  # major, minor and harmful
  assigned.lines = round(runif(combd.num, 2, max.line))
  max.line = max(assigned.lines)
  medicine.mat = matrix(vector("list",combd.num*max.line), nrow = combd.num, ncol = max.line)
  assigned.harmful.nums = round(seq(harmful.min.num, harmful.max.num, (harmful.max.num-harmful.min.num)/(max.line-1)))
  assigned.harmful.probs = round(seq(harmful.min.prob, harmful.max.prob, (harmful.max.prob-harmful.min.prob)/(max.line-1)), 2)
  # how many times extend the severity
  medicine.effect.mat = medicine.mat
  effect.range.arr = seq(0.35, 0.8, (0.8-0.35)/max.line)
  for(i in 1:combd.num){
    my.major.labs = combd.major.list[[i]][["labs"]] # len=1
    my.minor.labs = combd.minor.list[[i]][["labs"]] # len>=1
    my.harmful.cands = 1:labs.num
    my.harmful.cands = my.harmful.cands[-c(my.major.labs, my.minor.labs)]
    my.line.num = assigned.lines[i]
    my.harmful.orders = sort(sample(1:max.line, my.line.num))
    for(j in 1:my.line.num){
      cur.line = my.harmful.orders[j]
      cur.harmful.num = assigned.harmful.nums[cur.line]
      cur.harmful.prob = assigned.harmful.probs[cur.line]
      the.harmful.num = sample(c(0,1:cur.harmful.num), 1, prob = c((1-cur.harmful.prob),rep(cur.harmful.prob,cur.harmful.num)))
      if(the.harmful.num==0){
        the.harmful.labs = c()
        the.harmful.effect = c()
      } else {
        the.harmful.labs = sample(my.harmful.cands, the.harmful.num)
        #the.harmful.effect = runif(the.harmful.num, -1.1, -1.02) # harmful random range
        the.harmful.effect = runif(the.harmful.num, -0.3, -0.05) # harmful random range
      }
      if(length(my.minor.labs)==0){
        the.minor.effect = c()
      } else {
        the.minor.effect = runif(length(my.minor.labs), 0.05, 0.34) # minor random range
      }
      medicine.mat[[i,j]] = list(major.labs=my.major.labs, minor.labs=my.minor.labs, harmful.labs=the.harmful.labs)
      the.major.effect = runif(1,effect.range.arr[j],effect.range.arr[j+1])
      medicine.effect.mat[[i,j]] = list(major.effect=the.major.effect, minor.effect=the.minor.effect, harmful.effect=the.harmful.effect)
    }
  }
  cat("    rand.medicine() done.\n")
  return(list(medicine.mat=medicine.mat, medicine.effect.mat=medicine.effect.mat))
}

rand.pat.effect <- function(){
# Require: [ medicine.mat, pat.num, special.rate,
# randomize patients medicine effects, only applicable on major ]
  
  #medicine.mat: combd.num*max.line. NULL or list(major.labs,minor.labs,harmful.labs)
  
  combd.num = nrow(medicine.mat)
  max.line = ncol(medicine.mat)
  medicine.avlb.mat = !matrix(sapply(medicine.mat, is.null),combd.num,max.line)
  
  # randomize effectiveness
  library(fGarch)
  med.pat.effect.mat = matrix(vector("list",combd.num*max.line),combd.num,max.line)
  skewed.samples <- rnbinom(10000, 10, 0.5)
  skewed.samples = -(skewed.samples - median(skewed.samples))
  skewed.samples = skewed.samples - min(skewed.samples)
  skewed.samples = skewed.samples/mean(skewed.samples)
  # 0 ~ 1.5, median=1
  for(i in 1:combd.num){
    for(j in 1:max.line){
      if(medicine.avlb.mat[i,j]){
        pat.effect = sample(skewed.samples, pat.num, T)
        special.idxs = sample(1:pat.num, round(pat.num*special.rate))
        pat.effect[special.idxs] = - ((abs(pat.effect[special.idxs]-1)) * 0.9) # lower than 0.9
        med.pat.effect.mat[[i,j]] = pat.effect
      }
    }
  }
  
  cat("    rand.pat.effect() done.\n")
  return(med.pat.effect.mat)
}


get.treatment <- function(){
# Require: [value.init.list, pat.time.obs.mat,   combd.major.list, combd.minor.list, 
# rout.list, labs.list, maps.list, med.list, med.pat.effect.mat, pat.rout.list, crit.increase.rate]
  
  # assige medicine to patients
  # medicine can propost routines
  # medicine.effect.mat: contains [major,minor,harmful] effect*3 for each combd+line
  # med.pat.effect.mat: contains 1 rate per patient for each combd+line
  
  #value.init.list: len=labs.num, pat.num * year.num*2
  #combd.major.list: len=combd.num, array names: labs, mean, crit
  #combd.minor.list: len=combd.num, list names: labs, mean, crit
  #rout.list: rout.mat(rout.num*max.step), connect.mat(combd.num+1*combd.num+1), rout.major.devlp.mat(rout.num*max.step), rout.minor.devlp.mat(rout.num*max.step)
  #labs.list: labs.mean.arr, labs.sd.arr, labs.crit.arr. (len=labs.num)
  #maps.list: devlp.map(pat.num*labs.num), mean.map(pat.num*labs.num), crit.map(pat.num*labs.num)
  #medicine.mat: combd.num*max.line. NULL or list(major.labs,minor.labs,harmful.labs)
  #medicine.effect.mat: combd.num*max.line. 
  #  - NULL or list(major.effect[0.35, 0.8] increasing, minor.effect[0.05, 0.34], harmful.effect[-0.3, -0.05])
  #med.pat.effect.mat: combd.num*max.line. vector len=pat.num
  #pat.rout.assigned: len=pat.num
  
  direct.impact.tag = F
  # medicine whether directly impact labs

  labs.mean.arr = labs.list$labs.mean.arr
  labs.crit.arr = labs.list$labs.crit.arr
  
  rout.length = apply(rout.mat, 1, function(x){sum(!is.na(x))})
  
  pat.rout.assigned = pat.rout.list$pat.rout.assigned
  pat.start.assigned = pat.rout.list$pat.start.assigned
  pat.start.offset = pat.rout.list$pat.start.offset
  
  medicine.mat = med.list$medicine.mat
  medicine.effect.mat = med.list$medicine.effect.mat
  max.line = ncol(medicine.mat)
  medicine.avlb.mat = !matrix(sapply(medicine.mat, is.null),combd.num,max.line)
  
  # start - initialization
  value.treated.list = value.init.list
  
  pat.med.mat = matrix(vector("list",pat.num*year.num), pat.num, year.num) 
  # element:array(len=combd.num) or NULL, numeric 
  pat.dx.mat = matrix(vector("list",pat.num*year.num), pat.num, year.num) 
  # element:array(len=combd.num) or NULL, binary
  treatment.list = vector("list", year.num) # pat.num*combd.num each
  for (t in 1:length(treatment.list)) {
    treatment.list[[t]] = matrix(0, pat.num, combd.num)
  }

  ## patient X labs matrix, observable time for each labs
  #pat.labs.obs.mat = matrix(NA, nrow=pat.num, ncol=labs.num)
  
  # initialize delta data structure
  # delta defaulted filled 0
  outcome.arr = paste("labs",1:labs.num,sep="")
  delta.outcome.list = vector("list",year.num-1)
  for(t in 1:year.num){  # 1:year.num delta may happen in year=1
    delta.outcome.list[[t]] = matrix(0,pat.num,length(outcome.arr))
    colnames(delta.outcome.list[[t]]) = outcome.arr
  }
  delta.trj.list = vector("list",length(outcome.arr))
  names(delta.trj.list) = outcome.arr
  for(oc in 1:labs.num){
    delta.mat = matrix(0, pat.num, year.num) # 1:year.num
    delta.trj.list[[oc]] = delta.mat
  }
  
  update.devlp <- function(in.devlp, in.med.effect, in.pat.effect, effect.return=F){
    if(length(in.devlp)==0) return(in.devlp)
    cur.devlp = in.devlp
    the.effect = in.med.effect * in.pat.effect
    the.sevrt = -the.effect
    out.devlp = sevrt.2.devlp(cur.devlp, the.sevrt)
    
    if(effect.return){
      return(list(out.devlp, the.effect))
    } else {
      return(out.devlp)
    }
  }
  
  impact.labs <- function(in.labs.value.mat, in.labs, in.med.effect, in.pat.effect){
    labs.len = length(in.labs)
    time.len = ncol(in.labs.value.mat)
    if(labs.len!=length(in.med.effect)) return()
    if(!is.numeric(in.pat.effect) || length(in.pat.effect)!=1) return()
    if(labs.len==0) return(in.labs.value.mat)
    
    cur.lab.value = in.labs.value.mat[in.labs,,drop=F]
    random.offset.rate = runif(1,-0.3,0.3)
    random.offset = apply(cur.lab.value, 1, mean)*random.offset.rate
    cur.lab.value = cur.lab.value + matrix(random.offset, nrow=labs.len, ncol=time.len)

    the.effect = (in.med.effect * in.pat.effect) # may be negative
    # the.effect=1 means shrink 0.5 at most
    the.slope.mat = matrix(1, nrow=labs.len, ncol=time.len)
    for(l in 1:labs.len){
      the.slope.mat[l,] = 1 - (1-(1:time.len)/time.len)*(the.effect[l]*0.5)
    }
    the.slope.mat[the.slope.mat<=0] = 0.1 # just in case
    cur.lab.value = cur.lab.value * the.slope.mat
    
    in.labs.value.mat[in.labs,] = cur.lab.value
    return(in.labs.value.mat)
  }
  
  for(i in 1:pat.num){
    # initialize dx + med
    pat.dx.mat[i,1:year.num] = lapply(pat.dx.mat[i,1:year.num], function(x){rep(0,combd.num)})
    pat.med.mat[i,1:year.num] = lapply(pat.med.mat[i,1:year.num], function(x){rep(0,combd.num)})
    
    my.rout.id = pat.rout.assigned[i]
    if(my.rout.id==0){
      next;
    }
    my.rout = rout.mat[my.rout.id,];
    my.rout.start = pat.start.assigned[i]
    my.rout.offset = pat.start.offset[i]
    my.rout.offset = floor(my.rout.offset)
    my.rout.end = rout.length[my.rout.id]
    my.rout.major.labs = rep(0,my.rout.end)
    for (j in 1:my.rout.end) { # fill full routine
      my.rout.major.labs[j] = combd.major.list[[my.rout[j]]][["labs"]]
    }
    
    # pat.time.obs.mat is from initial function
    # most is af, contains some bf
    # change to all bf, 1 per year
    my.obs.time = which(pat.time.obs.mat[i,])# by year*2
    my.obs.time = my.obs.time[my.obs.time%%2==0]
    my.obs.time = my.obs.time-1 # change to bf
    # times of offset
    offset.times = my.rout.offset*2
    my.obs.time = my.obs.time + offset.times

    # devlp map initalized
    my.devlp.map = maps.list$devlp.map[i,]
    my.mean.map = maps.list$mean.map[i,]
    my.crit.map = maps.list$crit.map[i,]
    my.maps = list(devlp.map=my.devlp.map, mean.map=my.mean.map, crit.map=my.crit.map)
    # refill offset
    labs.value.mat = devlp.2.labs(in.maps=my.maps, start.yr=1, end.yr=year.num+my.rout.offset+rand.extd) 
    # labs.num * ((year.num+my.rout.offset+rand.extd)*2)
    
    for (j in my.rout.start:my.rout.end) { # my routine
      the.combd = my.rout[j]
      the.med.avlb = sum(medicine.avlb.mat[the.combd,])
      the.obs = my.obs.time
      the.med.line = 0 # illed and need the 1st line medcine
      while(length(the.obs)>0){
        the.major.labs = combd.major.list[[the.combd]][["labs"]]
        the.major.crit = combd.major.list[[the.combd]][["crit"]]*((1+crit.increase.rate)^the.med.line)
        the.values = labs.value.mat[the.major.labs,]
        the.values = the.values[the.obs]
        the.values = the.values + value.init.list[[the.major.labs]][i,(the.obs-offset.times)]
        
        keep.tag = the.values>=the.major.crit
        the.obs=the.obs[keep.tag]
        if(length(the.obs)==0) break;
        
        cur.obs.time = the.obs[1] # by year*2, must be the bf stage
        cur.obs.year = (cur.obs.time+1)/2 
        real.obs.time = cur.obs.time - offset.times
        real.obs.year = (real.obs.time+1)/2 
        if(the.med.line==0){ # the 1st time observed, get dx
          pat.dx.mat[i,real.obs.year:year.num] = lapply(pat.dx.mat[i,real.obs.year:year.num], function(x){x[the.combd]=1;x})
        }
        if(the.med.line+1>the.med.avlb){ # all possible medicine have been used
          break;
        }
        the.med.line = the.med.line+1
        # Here is the certain point of delta produced
        pat.med.mat[i,real.obs.year:year.num] = lapply(pat.med.mat[i,real.obs.year:year.num], function(x){x[the.combd]=the.med.line;x})
        treatment.list[[real.obs.year]][i,the.combd] = the.med.line
        # update pat.time.obs.mat, make the bf = T
        pat.time.obs.mat[i, real.obs.year*2-1] = T
        
        # update devlp.map
        med.labs.list = medicine.mat[[the.combd, the.med.line]]
        med.major.labs = med.labs.list$major.labs
        med.minor.labs = med.labs.list$minor.labs
        med.harmful.labs = med.labs.list$harmful.labs
        
        med.effect.list = medicine.effect.mat[[the.combd, the.med.line]]
        med.major.effect = med.effect.list$major.effect
        med.minor.effect = med.effect.list$minor.effect
        med.harmful.effect = med.effect.list$harmful.effect
        
        pat.major.effect = med.pat.effect.mat[[the.combd, the.med.line]][i]
        #pat.minor.effect = runif(length(med.minor.labs), med.major.effect*0.5, med.major.effect*1.5) 
        # randomized personal minor effect
        #pat.harmful.effect = runif(length(med.harmful.labs), med.major.effect*0.5, med.major.effect*1.5)
        # randomized personal harmful effect

        # major effect
        # assert: med.major.labs == the.major.labs
        my.devlp.map[med.major.labs] = update.devlp(in.devlp=my.devlp.map[med.major.labs], in.med.effect=med.major.effect,
                                        in.pat.effect=pat.major.effect, effect.return=F)
        # minor effect
        my.devlp.map[med.minor.labs] = update.devlp(in.devlp=my.devlp.map[med.minor.labs], in.med.effect=med.minor.effect,
                                    in.pat.effect=1, effect.return=F)
        # harmful effect
        my.devlp.map[med.harmful.labs] = update.devlp(in.devlp=my.devlp.map[med.harmful.labs], in.med.effect=med.harmful.effect,
                                    in.pat.effect=1, effect.return=F)

        
        my.maps = list(devlp.map=my.devlp.map, mean.map=my.mean.map, crit.map=my.crit.map)
        cur.labs.value = devlp.2.labs(in.maps=my.maps, start.yr=cur.obs.year, end.yr=year.num+my.rout.offset+rand.extd)
        
        # if(i==3 && the.major.labs==6){
        #   cat("debug\n")
        # }
        
        
        # medicine impact labs directly
        if(direct.impact.tag){
          cur.labs.value = impact.labs(cur.labs.value, in.labs=med.major.labs, 
                                       in.med.effect=med.major.effect, in.pat.effect=pat.major.effect)
          cur.labs.value = impact.labs(cur.labs.value, in.labs=med.minor.labs, 
                                      in.med.effect=med.minor.effect, in.pat.effect=1)
          cur.labs.value = impact.labs(cur.labs.value, in.labs=med.harmful.labs, 
                                       in.med.effect=med.harmful.effect, in.pat.effect=1)
        }
        
        # change starts from the time after cur.obs.time
        labs.value.mat[,(cur.obs.time+1):ncol(labs.value.mat)] = cur.labs.value[,2:ncol(cur.labs.value)]
        cur.delta.vec = cur.labs.value[,2] # labs.num length
        # record delta: maybe major/minor time-conflict among different combd, addup to form combined effect
        exist.delta.vec = delta.outcome.list[[real.obs.year]][i, ]
        delta.outcome.list[[real.obs.year]][i, ] = exist.delta.vec+cur.delta.vec # ncol=labs.num
        
        # next obs time point
        the.obs = the.obs[-c(1)]
      } # end while
      
      for (l in 1:labs.num) {
        cur.values = value.init.list[[l]][i,]
        cur.values = cur.values + labs.value.mat[l,(offset.times+1):(offset.times+year.num*2)]
        value.treated.list[[l]][i,] = cur.values
      }
      
    } # j=my.rout.start:my.rout.end
    
  }# i=1:pat.num
  
  # synchronize delta data structure
  for (t in 1:year.num) {
    delta.outcome = delta.outcome.list[[t]] # pat.num*labs.num
    for (oc in outcome.arr) {
      delta.trj.list[[oc]][,t] = delta.outcome[,oc]
    }
  }
  
  cat("    get.treatment() done.\n")
  return(list(value.treated.list, pat.time.obs.mat, pat.dx.mat, pat.med.mat, 
              treatment.list, #pat.labs.obs.mat,
              delta.trj.list, delta.outcome.list))
}







