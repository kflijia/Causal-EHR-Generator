
split.tr.ts = function(in.num, tr.prop){
  tr.idxs = sort(sample(in.num, tr.prop*in.num));
  ts.idxs = which(! 1:in.num %in% tr.idxs);
  return(list(tr.idxs, ts.idxs));
}

split.tr.tags = function(in.num, tr.prop){
  tr.idxs = sort(sample(in.num, tr.prop*in.num));
  tr.tags = 1:in.num %in% tr.idxs;
  tr.tags
}


get.mse<-function(yhat,y){
  mse = sum((yhat-y)^2) * (1/length(yhat))
  return(mse)
}


get.auc <- function(yhat, y, get.tag.tag=F) {
  
  y.len = length(y);
  y.posnum = sum(y);
  y.negnum = sum(y==0);
  
  auc.tag = 1;
  if (y.len>1 && y.posnum>0 && y.negnum>0) {
    R1 <- sum(rank(yhat)[y==1]);
    n1 <- length(which(y==1));
    U1 <- R1 - n1*(n1+1)/2;
    the.auc = U1/(n1*(length(y)-n1));
  } else {
    the.auc = mean(abs(yhat - !y));
    auc.tag = 0;
  }
  
  if (get.tag.tag){
    return(list(the.auc, auc.tag));
  } else {
    return(the.auc);
  }
}


get.glm = function(tr.list, ts.list, terms.arr=c(), baseline.mod=NULL, crossvalid.num=0){
  tr.X = as.matrix(tr.list$X);
  tr.Y = as.matrix(tr.list$Y);
  tr.data = cbind(tr.X, tr.Y);
  tr.num = nrow(tr.X);
  
  ts.X = as.matrix(ts.list$X);
  ts.Y = as.matrix(ts.list$Y);
  ts.data = cbind(ts.X, ts.Y);
  
  if(length(terms.arr)==0){ # all terms
    the.formula = as.formula("label~.");
    if(!is.null(baseline.mod)){ # have offset
      the.formula = as.formula("label~.-o");
    }
  } else { # selected terms
    the.formula = as.formula(paste("label~",paste(terms.arr,collapse="+"),sep=""));
  }
  
  crossvalid.tag=F; # no crossvalidation
  if(is.numeric(crossvalid.num) && crossvalid.num>0){
    crossvalid.tag=T;
    split.idxs = rep(1:crossvalid.num,ceiling(tr.num/crossvalid.num))[1:tr.num];
  }
  
  if(!crossvalid.tag) { # no crossvalidation
    if(is.null(baseline.mod)){ # no offset
      glm.mod = glm(the.formula, family = "binomial", data = data.frame(tr.data))
      ts.yhatm <- predict(glm.mod, as.data.frame(ts.X), type="response")
    } else { # have offset
      tr.offset <- predict(baseline.mod, data.frame(tr.data), type="link")
      ts.offset <- predict(baseline.mod, data.frame(ts.data), type="link")
      glm.mod = glm(the.formula, family = "binomial", data = data.frame(tr.data, o=tr.offset), offset=o)
      ts.yhatm <- predict(glm.mod, data.frame(ts.X, o=ts.offset), type="response")
    }
  } else { # crossvalidation
    vl.auc.arr = c();
    glm.mod.list = list();
    for(i in 1:crossvalid.num){
      my.vl.data = tr.data[split.idxs==i,,drop=F];
      my.tr.data = tr.data[split.idxs!=i,,drop=F];
      if(is.null(baseline.mod)){ # no offset
        my.glm.mod = glm(the.formula, family = "binomial", data = data.frame(my.tr.data))
        my.vl.yhatm <- predict(my.glm.mod, data.frame(my.vl.data), type="response")
      } else { # have offset
        my.tr.offset <- predict(baseline.mod, data.frame(my.tr.data), type="link")
        my.vl.offset <- predict(baseline.mod, data.frame(my.vl.data), type="link")
        my.glm.mod = glm(the.formula, family = "binomial", data = data.frame(my.tr.data, o=my.tr.offset), offset=o)
        my.vl.yhatm <- predict(my.glm.mod, data.frame(my.vl.data, o=my.vl.offset), type="response")
      }
      my.vl.auc = get.auc(my.vl.yhatm, my.vl.data[,ncol(my.vl.data)]);
      glm.mod.list[[i]] = my.glm.mod;
      vl.auc.arr[i] = my.vl.auc;
    }
    glm.mod = glm.mod.list[[which.max(vl.auc.arr)]];
    if(is.null(baseline.mod)){ # no offset
      ts.yhatm <- predict(glm.mod, as.data.frame(ts.X), type="response")
    } else { # have offset
      ts.offset <- predict(baseline.mod, data.frame(ts.data), type="link")
      ts.yhatm <- predict(glm.mod, data.frame(ts.X, o=ts.offset), type="response")
    }
  }
  
  tr.yhatm = glm.mod$fitted.values;
  tr.auc = get.auc(tr.yhatm, tr.Y);
  coeffs = as.matrix(data.frame(glm.mod$coefficients));
  tr.dev = glm.mod$deviance
  tr.nulldev = glm.mod$null.deviance;
  tr.dev.ratio = 1 - tr.dev/tr.nulldev;
  tr.resid = resid(glm.mod);
  tr.res = list(tr.auc=tr.auc, tr.yhatm=tr.yhatm, tr.dev=tr.dev, tr.nulldev=tr.nulldev, tr.dev.ratio=tr.dev.ratio, tr.resid=tr.resid);
  # not returned
  
  ts.auc = get.auc(ts.yhatm, ts.Y);
  a = ifelse(ts.Y==0, 1, ts.Y/ts.yhatm);
  b = ifelse(ts.Y==1, 1, (1-ts.Y)/(1-ts.yhatm));
  ts.resid = ((ts.Y-ts.yhatm)/abs(ts.Y-ts.yhatm)) * 
    ( 2*ts.Y*log(a)+2*(1-ts.Y)*log(b) )^0.5;
  ts.dev = sum(ts.resid^2);
  #temp.res = glm(label~., family = "binomial", data = as.data.frame(ts.data));
  #ts.nulldev = temp.res$null.deviance;
  #ts.dev.ratio = 1 - ts.dev/ts.nulldev;
  #ts.res = list(ts.auc=ts.auc, ts.yhatm=ts.yhatm, ts.dev=ts.dev, ts.nulldev=ts.nulldev, ts.dev.ratio=ts.dev.ratio, ts.resid=ts.resid);
  ts.res = list(ts.auc=ts.auc, ts.yhatm=ts.yhatm, ts.dev=ts.dev, ts.resid=ts.resid);
  
  return(list(coeffs, ts.res, glm.mod));
}

get.glmnet = function(tr.list, ts.list){
  tr.X = as.matrix(tr.list$X);
  tr.Y = as.matrix(tr.list$Y);
  
  ts.X = as.matrix(ts.list$X);
  ts.Y = as.matrix(ts.list$Y);
  
  glmnet.mod <- glmnet(tr.X, tr.Y, family="binomial");
  yhatm <- predict(glmnet.mod, tr.X, type="response"); # matrix:length(ts.sub)*#lambdas
  my.aucs = c();
  for (j in 1:ncol(yhatm)){
    my.aucs = c(my.aucs, get.auc(yhatm[,j], tr.Y));
  }
  aucs <- my.aucs;
  lambdas <- glmnet.mod$lambda;
  if (length(lambdas)<=10) {
    #l = lambdas[which.max(aucs)];
    l.idx = which.max(aucs);
  } else {
    my.loess = loess(aucs~lambdas);
    if (all(is.na(my.loess[[2]]))) {
      #l = lambdas[which.max(aucs)];
      l.idx = which.max(aucs);
    } else {
      #l <- lambdas[which.max(predict(my.loess, lambdas))];
      l.idx <- which.max(predict(my.loess, lambdas))
    }
  }
  
  coeffs = coef(glmnet.mod)[,l.idx,drop=F];
  coeffs = as.matrix(coeffs);
  tr.auc = aucs[l.idx];
  l = lambdas[l.idx];
  ts.fit <- predict(glmnet.mod, ts.X, type="response", s=l);
  ts.auc = get.auc(ts.fit, ts.Y);
  
  linear.predict = (cbind(rep(1,nrow(ts.X)),ts.X) %*% coeffs);
  
  return(list(coeffs, tr.auc, ts.auc, l, ts.fit, linear.predict));
}


boost.fit <- function(x, y, beta.init=rep(0, ncol(x)+1), maxiter=10000, gamma=1, epsilon=1e-10) {
  intercept <- beta.init[1]
  beta      <- beta.init[2:(ncol(x)+1)]
  
  u <- intercept + as.vector(x %*% beta);
  p <- exp(u) / (1+exp(u));
  prev.nll <- sum(y*u - log(1+exp(u)));  # log likelihood
  cat("Initial NLL: ", -prev.nll, "\n")
  
  for(iter in 1:maxiter) {
    z <- -(y-p)
    cors <- apply(x, 2, cov, y=z);
    idx  <- which.max(abs(cors));
    
    beta[idx] <- beta[idx] - gamma * sum(z*x[,idx])/sum(x[,idx]*x[,idx])
    u <- intercept + as.vector(x %*% beta)
    p <- exp(u) / (1+exp(u))
    nll.coef <- sum(y*u - log(1+exp(u)))
    z <- -(y-p);
    
    intercept <- intercept - gamma * sum(z)/nrow(x)
    u <- intercept + as.vector(x %*% beta)
    p <- exp(u) / (1+exp(u))
    
    nll <- sum(y*u - log(1+exp(u)))
    if(iter < 10 | iter < 100 & iter %% 10 == 0 | iter < 1000 & iter %% 100 == 0 | iter %% 1000 == 0) { 
      #cat("iter=", iter, "nll=", -nll, ifelse(nll > prev.nll, "(*)", ""), "\n")
    }
    if(nll - prev.nll < epsilon) { 
      cat("iter=", iter, "nll=", -nll, "nll.coef=", nll.coef, "prev.nll=", prev.nll, "... converged.\n");
      break; 
    }
    prev.nll <- nll;                
  }
  
  return(list(intercept, beta, p, u));
}


get.boost = function(tr.list, ts.list, in.beta.init, in.gamma = 1){
  tr.X = as.matrix(tr.list$X);
  tr.Y = as.vector(tr.list$Y);
  
  ts.X = as.matrix(ts.list$X);
  ts.Y = as.vector(ts.list$Y);
  
  fit.res = boost.fit(tr.X, tr.Y, beta.init=in.beta.init, gamma=in.gamma);
  intcpt = fit.res[[1]];
  beta = fit.res[[2]];
  tr.fit = fit.res[[3]];
  
  ts.linear.fit <- intcpt + as.vector(ts.X %*% beta);
  ts.fit <- exp(ts.linear.fit) / (1+exp(ts.linear.fit));
  
  tr.auc = get.auc(tr.fit, tr.Y);
  ts.auc = get.auc(ts.fit, ts.Y);
  
  coeffs = c(intcpt, beta);
  
  # tr not returned 
  
  a = ifelse(ts.Y==0, 1, ts.Y/ts.fit);
  b = ifelse(ts.Y==1, 1, (1-ts.Y)/(1-ts.fit));
  ts.resid = ((ts.Y-ts.fit)/abs(ts.Y-ts.fit)) * 
    ( 2*ts.Y*log(a)+2*(1-ts.Y)*log(b) )^0.5;
  ts.dev = sum(ts.resid^2);
  ts.res = list(ts.auc=ts.auc, ts.yhatm=ts.fit, ts.dev=ts.dev, ts.resid=ts.resid);
  
  return(list(coeffs, ts.res));
}


