# ------------------------------------------------------
# auc
#
auc <- function(yhat, y, na.rm=FALSE) {
  
  if(na.rm) {
    rows <- which(!is.na(yhat) & !is.na(y))
    yhat <- yhat[rows]
    y <- y[rows]
  }
  
  R1 <- sum(rank(yhat)[y==1]);
  n1 <- length(which(y==1));
  U1 <- R1 - n1*(n1+1)/2;
  U1/(n1*(length(y)-n1));
}

# ------------------------------------------------------
# loessPlot
#
loessPlot <- function(x, y, span=.75, add=FALSE, col="red", lwd=2, maxrows=2000, xlim=NULL, ylim=NULL, 
                      main=NULL, xlab=NULL, ylab=NULL, pointcol=NULL, pointpch=NULL, ...) {
  
  if(length(x) > maxrows) { rows <- sample(1:length(x), maxrows, replace=FALSE); }
  else { rows <- 1:length(x); }
  
  fit <- loess(y[rows]~x[rows], span=span);
  
  if(!add) { plot(x, y, xlim=xlim, ylim=ylim, main=main, col=pointcol, pch=pointpch, xlab=xlab, ylab=ylab); }
  
  xs <- seq(min(x), max(x), length.out=maxrows)
  lines(xs, predict(fit, xs), lwd=lwd, col=col, ...); 
}


# ------------------------------------------------------
# backwards elimination
#

library(survival)

# coxph.backelim
#
# Performs backwards elimination.
# x: model matrix for the predictors
# s: Surv object having columns 'time' and 'status'
#
# Returns the model matrix of the resultant model.
# The model can be built as
# coxph(Surv(time, status)~., xres)
# where xres is the resultant model matrix.
#
coxph.backelim <- function(x, s, alpha=.05, offset=NULL, ...) {
  
  if(nrow(x) != nrow(s)) { stop("Mismatch in number of rows: x:", nrow(x), "s:", nrow(s)); }
  
  cols <- colnames(x);
  
  while(TRUE) {
    
    # Build the model
    if(!is.null(offset)) {
      frm <- paste("Surv(time, status)~.", "+offset(", offset, ")-", offset);
    } else {
      frm <- "Surv(time, status)~.";
    }
    fit <- coxph(as.formula(frm), data.frame(x[,cols], time=s[,"time"], status=s[,"status"]), ...)
    pvals <- summary(fit)$coef[,5]
    rn <- rownames(summary(fit)$coef);
    
    # Remove NAs
    naidxs <- which(is.na(pvals));
    if(length(naidxs) > 0) {
      cat("Removing ", rn[naidxs], "for producing NA coefficients\n");
      cols <- setdiff(cols, rn[naidxs]);
      next;
    }
    
    # Find the highest p-value
    if(max(pvals,na.rm=T) < alpha) { break; }
    max.pval.idx <- which.max(pvals);
    cat("Removing ", rn[max.pval.idx], "with p-value", pvals[max.pval.idx], "\n");
    cols <- setdiff(cols, rn[max.pval.idx]);
    
  }
  
  cols
}


lm.backelim <- function(x, y, alpha=.05, verbose=0, ...) {
  
  if(nrow(x) != length(y)) { stop("Mismatch in number of rows: x:", nrow(x), "y:", length(y)); }
  
  cols <- colnames(x);
  
  while(TRUE) {
    
    # Build the model
    fit <- lm(y~., data.frame(x[,cols,drop=FALSE], y=y), ...);
    pvals <- summary(fit)$coef[,4]
    rn <- rownames(summary(fit)$coef);
    pvals[1] <- 0;   # Forcefully keep the intercept
    
    # Remove NAs
    nacols <- setdiff(cols, rn)
    if(length(nacols) > 0) {
      if(verbose > 0) { cat("Removing ", nacols, "for producing NA coefficients\n"); }
      cols <- setdiff(cols, nacols);
      next;
    }
    
    # Find the highest p-value
    if(max(pvals,na.rm=T) < alpha) { break; }
    max.pval.idx <- which.max(pvals);
    if(verbose > 0) { cat("Removing ", rn[max.pval.idx], "with p-value", pvals[max.pval.idx], "\n"); }
    cols <- setdiff(cols, rn[max.pval.idx]);
    
    if(length(cols) < 1) { break; }
    
  }
  
  cols
}


glm.backelim <- function(x, y, alpha=.05, family="binomial", offset=NULL, verbose=0, ...) {
  
  if(nrow(x) != length(y)) { stop("Mismatch in number of rows: x:", nrow(x), "y:", length(y)); }
  if(!is.null(offset) && length(offset) != length(y)) { stop("Mismatch in number of rows: offset:", length(offset), "y:", length(y)); }
  
  cols <- colnames(x);
  
  # Remove the cols with zero variability
  sds <- apply(x, 2, sd, na.rm=TRUE);
  
  while(TRUE) {
    
    # Build the model
    if(!is.null(offset)) {
      fit <- glm(y~.+offset(o)-o, data.frame(x[,cols,drop=FALSE], y=y, o=offset), family=family, ...);
    } else { 
      fit <- glm(y~., data.frame(x[,cols,drop=FALSE], y=y), family=family, ...);
    }
    pvals <- summary(fit)$coef[,4]
    rn <- rownames(summary(fit)$coef);
    pvals[1] <- 0;   # Forcefully keep the intercept
    
    # Remove NAs
    nacols <- names(which(is.na(coef(fit))));
    if(length(nacols) > 0) {
      if(verbose > 0) { cat("Removing ", nacols, "for producing NA coefficients\n"); }
      cols <- setdiff(cols, nacols);
      next;
    }
    
    # Find the highest p-value
    if(max(pvals,na.rm=T) < alpha) { break; }
    max.pval.idx <- which.max(pvals);
    if(verbose > 0) { cat("Removing ", rn[max.pval.idx], "with p-value", pvals[max.pval.idx], "\n"); }
    cols <- setdiff(cols, rn[max.pval.idx]);
    
  }
  
  cols
}

# ---------------------------------------------------------------------------
# ml.fit
#
# Builds a GBM model for a data set. 
# The function determines the optimal number of trees (through a 30% leave-out validation set)
# and builds a model for the entire data set using the optimal number of trees.
#
# x, y:         data set. 'y' is binary encoded as 0/1. (NA's are ignored)
# tree.seq:     try these numbers of trees to find the optimal number.
#
# Returns a GBM object with three additional slots:
# auc:          validation AUC at the optimal number of trees
# aucs:         AUC for all numbers of trees in tree.seq
# n.trees:      the optimal number of trees.
#
ml.fit <- function(x, y, tree.seq=c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000)) {
  
  pos <- which(y==1);
  neg <- which(y==0);
  #	cat("DEBUG: ml.fit: y: pos:", length(pos), "neg:", length(neg), "NA:", length(is.na(y)), "\n");
  
  # Create training and validation sets
  pos.tr <- sample(pos, .7*length(pos));  # replace=FALSE
  neg.tr <- sample(neg, .7*length(neg));
  pos.vd <- setdiff(pos, pos.tr);
  neg.vd <- setdiff(neg, neg.tr);
  tr <- c(pos.tr, neg.tr);
  vd <- c(pos.vd, neg.vd);
  rows <- c(tr, vd);
  
  # Training
  fit <- gbm(y~., data=data.frame(y=y[tr], x[tr,]), distribution="bernoulli", n.trees=max(tree.seq));
  
  # Select the optimal number of trees
  yhatm <- predict(fit, data.frame(x[vd,]), type="response", n.tree=tree.seq);
  aucs <- apply(yhatm, 2, auc, y=y[vd]);
  idx <- which.max(aucs);
  n.trees <- tree.seq[idx];
  #	cat("DEBUG: ml.fit: aucs:", aucs, "\n");
  
  # Rebuild the model on the entire data set with the optimal number of trees
  fit <- gbm(y~., data=data.frame(y=y[rows], x[rows,]), distribution="bernoulli", n.trees=n.trees);
  fit$n.trees <- n.trees;
  fit$aucs <- aucs;
  fit$auc <- aucs[idx];
  fit
}

