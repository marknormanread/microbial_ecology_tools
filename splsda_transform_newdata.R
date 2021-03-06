#############################################################################################################
# Adapted from the mixOmics package github distribution. 
#
# Aiming to use sPLS-DA as a feature selection scheme, wanting to then build classifiers outside of the mixomics 
# framework, we need a way of transforming new data (e.g. test data in a cross validation loop) into the sPLS 
# representation constructed on the training data. 
#
# Mixomics does not explictly expose this information, it is instead folded within a broader "predict" framework. 
# We don't want the predict functionality, but DO want the transformation of new data that prediction relies on. 
# This functionality is contained within the mixOmics `predict.mint.block.pls.R` module. 
#
# The present module is that, modified only to cut out prediction and return the transformed data. 
#
# Mark N. Read, 2019
#############################################################################################################

#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 27-05-2015
# last modified: 05-10-2017
#
# Copyright (C) 2015
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################

# ========================================================================================================
# This function makes a prediction of a 'newdata' by using the results of 'object'.
# Depending on the class of the object (mint).(block).(s)pls(da) (16 classes so far), the input data is different
# and the preparation of the data is different - different scaling for instance if object="mint...."
# However, the prediction formula is the same for all classes, thus only one code
# ========================================================================================================

#predict.mixOmics <-
#predict.pls <-  predict.spls<- predict.plsda <- predict.splsda <-
#predict.mint.pls <- predict.mint.spls <- predict.mint.plsda <- predict.mint.splsda <-
#predict.block.pls <- predict.block.spls <- predict.block.plsda <- predict.block.splsda <-
#predict.mint.block.pls <- predict.mint.block.spls <- predict.mint.block.plsda <- predict.mint.block.splsda <- #predict.sgcca <-

# note FR: 27/05/16: the way I dealt with missing block to predict will probably not work for mint.block analysis


# object: any object (mint).(block).(s)pls(da)
# newdata: data to predict, same form as what's in object$X
# study.test: optional, used if a mint object.  grouping factor indicating which samples of `newdata' are from the same study
# dist: prediction distance to use
# multilevel: if there was a multilevel in object, one should be used for `newdata'
# newdata.scale: hidden parameter. used in tune/perf functions. pass `newdata' that has been scale already (gain in computational time)
# misdata.all: hidden parameter. used in tune/perf functions. any missing values in object$X
# is.na.X: hidden parameter. used in tune/perf functions. where are the missing values in object$X
# is.na.newdata: hidden parameter. used in tune/perf functions. where are the missing values in `newdata'
# noAveragePredict: hidden parameter. used in tune/perf functions. no calculation of Average Prediction (gain in computational time)
library(mixOmics)

splsda_transform <-
  function(object, newdata,study.test,dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), multilevel = NULL, ...)
  {
    time=FALSE  # TODO delete?
    
    # pass R check
    newdata.scale = misdata.all = is.na.X = is.na.newdata = noAveragePredict = NULL
    
    # input parameter: noAveragePredict=> no averagePredict calculation, used in tune.block.splsda
    
    if(is(object, c("rgcca","sparse.rgcca")))
      stop("no prediction for RGCCA methods")
    
    #check on dist
    if (!any(dist %in% c("all", "max.dist", "centroids.dist", "mahalanobis.dist")) & is(object, "DA"))
      stop("ERROR : choose one of the four following modes: 'all', 'max.dist', 'centroids.dist' or 'mahalanobis.dist'")
    #end general checks
    
    ncomp = object$ncomp
    
    if(is(object, "DA")) # a DA analysis (mint).(block).(s)plsda
    {
      #if DA analysis, the unmap Y is in ind.mat
      Y.factor=object$Y
      Y=object$ind.mat
    }else{
      #if not DA, Y is in object$Y
      Y=object$Y
      if(is.null(Y)) # block analysis
      {
        Y=object$X[[object$indY]]
      }
    }
    q=ncol(Y)
    
    if(time) time1 = proc.time()
    #print(names(list(...)))
    
    if(hasArg(newdata.scale)) # if an input `newdata.scale' is given, we use that one and we don't scale the data and don't do the logratio or multilevel
    {
      newdata = list(...)$newdata.scale
      object$logratio = NULL
      multilevel = NULL
    }
    
    mint.object = c("mint.pls", "mint.spls", "mint.plsda", "mint.splsda")
    block.object = c("block.pls", "block.spls", "block.plsda", "block.splsda")
    ### if the object is a block, the input newdata is different, we check newdata, make sure it's a list and check newdata/X
    if(!is(object, block.object)) # not a block (pls/spls/plsda/splsda/mint...)
    {
      p=ncol(object$X)
      if(is.list(object$X))
        stop("Something is wrong, object$X should be a matrix and it appears to be a list") #this should never happen/intern check
      
      if(is.list(newdata) & !is.data.frame(newdata))
        stop("'newdata' must be a numeric matrix")
      
      # deal with near.zero.var in object, to remove the same variable in newdata as in object$X (already removed in object$X)
      if(length(object$nzv$Position) > 0)
        newdata = newdata[, -object$nzv$Position,drop=FALSE]
      
      if(all.equal(colnames(newdata),colnames(object$X))!=TRUE)
        stop("'newdata' must include all the variables of 'object$X'")
      
      #not a block, the input newdata should be a matrix
      if (length(dim(newdata)) == 2) {
        if (ncol(newdata) != p)
          stop("'newdata' must be a numeric matrix with ncol = ", p,
               " or a vector of length = ", p, ".")
      }
      
      if (length(dim(newdata)) == 0) {
        if (length(newdata) != p)
          stop("'newdata' must be a numeric matrix with ncol = ", p,
               " or a vector of length = ", p, ".")
        dim(newdata) = c(1, p)
      }
      
      #check col/rownames of newdata
      check=Check.entry.single(newdata, ncomp,q=1)
      newdata=check$X
      
      if(length(rownames(newdata))==0) rownames(newdata)=1:nrow(newdata)
      if(max(table(rownames(newdata)))>1) stop('samples should have a unique identifier/rowname')
      
      # we transform everything in lists
      X=list(X=object$X)
      object$X=X
      newdata=list(newdata=newdata)
      
      object$indY=2
      ind.match = 1
      
    }else{
      
      # a block, newdata should be a list, each blocks should have the same number of samples
      if(!is.list(newdata))
        stop("'newdata' should be a list")
      
      if(!is.null(object$indY) && !is(object, "DA")) #if DA, object$X is already without Y
      {
        X = object$X[-object$indY]
      }else{
        X=object$X
      }
      object$X=X
      
      p = lapply(X, ncol)
      
      #matching newdata and X
      
      # error if no names on keepX or not matching the names in X
      if(length(unique(names(newdata)))!=length(newdata) | sum(is.na(match(names(newdata),names(X)))) > 0)
        stop("Each entry of 'newdata' must have a unique name corresponding to a block of 'X'")
      
      # error if not same number of samples in each block of newdata
      if (length(unique(sapply(newdata,nrow))) != 1)
        stop("All entries of 'newdata' must have the same number of rows")
      
      # I want to match keepX to X by names
      ind.match = match(names(X), names(newdata))
      if(any(is.na(ind.match)))
        warning("Some blocks are missing in 'newdata'; the prediction is based on the following blocks only: ", paste(names(X)[!is.na(ind.match)],collapse=", "))
      
      #replace missing blocks by 0
      newdataA = list()
      for (q in 1:length(X))
      {
        
        if (!is.na(ind.match[q])) # means there is a newdata with the same name as X[q] #(q <= length(newdata))
        {
          newdataA[[q]] = newdata[[ind.match[q]]]
        }else{
          newdataA[[q]] = matrix(0,nrow = unique(sapply(newdata,nrow)), ncol = ncol(X[[q]]), dimnames = list(rownames(newdata[[1]]), colnames(X[[q]]))) # if missing, replaced by a 0 matrix of correct size
        }
      }
      names(newdataA) = names(X)
      newdata = newdataA
      
      # deal with near.zero.var in object, to remove the same variable in newdata as in object$X (already removed in object$X)
      if(!is.null(object$nzv))
      {
        newdata = lapply(1:(length(object$nzv)-1),function(x){if(length(object$nzv[[x]]$Position>0)) {newdata[[x]][, -object$nzv[[x]]$Position,drop=FALSE]}else{newdata[[x]]}})
      }
      if(length(newdata)!=length(object$X)) stop("'newdata' must have as many blocks as 'object$X'")
      
      names(newdata)=names(X)
      
      if (any(lapply(newdata, function(x){length(dim(x))}) != 2)) {
        if (any(unlist(lapply(newdata, ncol)) != unlist(p)))
          stop("'newdata' must be a list with ", length(p), " numeric matrix and ncol respectively equal to ", paste(p, collapse = ", "), ".")
      }
      
      if (any(lapply(newdata, function(x){length(dim(x))}) == 0)) {
        if (any(unlist(lapply(newdata, ncol)) != unlist(p)))
          stop("'newdata' must be a list with ", length(p), " numeric matrix and ncol respectively equal to ", paste(p, collapse = ", "), ".")
        dim(newdata) = c(1, p) #don't understand that, Benoit?
      }
      
      #check dimnames and ncomp per block of A
      for(q in 1:length(newdata))
      {
        check=Check.entry.single(newdata[[q]], ncomp[q],q=q)
        newdata[[q]]=check$X
      }
      names(newdata)=names(X)
      
      #check that newdata and X have the same variables
      if(all.equal(lapply(newdata,colnames),lapply(X,colnames))!=TRUE)
        stop("Each 'newdata[[i]]' must include all the variables of 'object$X[[i]]'")
      
      #need to reorder variates and loadings to put 'Y' in last
      if(!is.null(object$indY))
      {
        indY=object$indY
        object$variates=c(object$variates[-indY],object$variates[indY])
        object$loadings=c(object$loadings[-indY],object$loadings[indY])
      }
      
      
    }
    
    # logratio and multilevel transform if necessary
    if (!is.null(object$logratio))
      newdata = lapply(newdata, logratio.transfo, logratio = object$logratio)
    
    if(!is.null(multilevel))
      newdata = lapply(newdata, withinVariation, design = data.frame(multilevel))
    
    p = lapply(X, ncol)
    q = ncol(Y)
    J = length(X) #at this stage we have a list of blocks
    variatesX = object$variates[-(J + 1)];  # Coordinates of instances in spls space. 
    loadingsX = object$loadings[-(J + 1)]
    
    scale = object$scale # X and Y are both mean centered by groups and if scale=TRUE they are scaled by groups
    
    
    ### if the object is not a mint analysis, the input study.test is missing and we can go faster to scale the data
    
    # we only scale the data if the input `newdata.scale' is missing. Only use in tune functions where we do not need to scale for every prediction as we predict on the same newdata for a grid of keepX
    if(!hasArg(newdata.scale))
    {
      if(!is(object, mint.object))#| nlevels(factor(object$study))<=1) #not a mint object or just one level in the study
      {   # not a mint (pls/spls/plsda/splsda/block...)
        
        # scale newdata if just one study
        if (!is.null(attr(X[[1]], "scaled:center")))
          newdata[which(!is.na(ind.match))] = lapply(which(!is.na(ind.match)), function(x){sweep(newdata[[x]], 2, STATS = attr(X[[x]], "scaled:center"))})
        if (scale)
          newdata[which(!is.na(ind.match))] = lapply(which(!is.na(ind.match)), function(x){sweep(newdata[[x]], 2, FUN = "/", STATS = attr(X[[x]], "scaled:scale"))})
        
        means.Y = matrix(attr(Y, "scaled:center"),nrow=nrow(newdata[[1]]),ncol=q,byrow=TRUE);
        if (scale)
        {sigma.Y = matrix(attr(Y, "scaled:scale"),nrow=nrow(newdata[[1]]),ncol=q,byrow=TRUE)}else{sigma.Y=matrix(1,nrow=nrow(newdata[[1]]),ncol=q)}
        concat.newdata=newdata
        names(concat.newdata)=names(X)
        
      }else{
        # a mint analysis
        
        #check study.test
        if(missing(study.test))
        {
          if(nlevels(object$study)==1)
          {study.test=factor(rep(1,nrow(newdata[[1]])))}else{
            stop("'study.test' is missing")}
        }else{
          study.test=as.factor(study.test)
        }
        if (any(unlist(lapply(newdata, nrow)) != length(study.test)))
          stop(paste0("'study' must be a factor of length ",nrow(newdata[[1]]),"."))
        
        #scale newdata if more than one study. If some levels of study.test are included in study.learn, the means/sigmas of study.learn are used to scale
        M=nlevels(study.test)
        study.learn=factor(object$study)
        names.study.learn=levels(study.learn);names.study.test=levels(study.test)
        match.study=match(names.study.test,names.study.learn)
        match.study.indice=which(!is.na(match.study))
        
        newdata.list.study = lapply(newdata,study_split,study.test) #a list of lists. [[j]][[m]]: j is for the blocks, m is for the levels(study)
        
        # each study is normalized, depending on how the normalization was performed on the learning set (center/scale)
        newdata.list.study.scale.temp =NULL#vector("list",length=J) #list(list())
        concat.newdata  = vector("list",length=J)
        for(j in 1:J) # loop on the blocks
        {
          for (m in 1:M) #loop on the groups of study.test
          {
            # case 1: sample test is part of the learning data set
            if(m%in%match.study.indice) #if some study of study.test were already in the learning set
            {
              
              if(scale==TRUE)
              {
                if(nlevels(object$study)>1)
                {
                  newdata.list.study.scale.temp=scale(newdata.list.study[[j]][[m]], center = attr(X[[j]],paste0("means:", levels(study.test)[m])), scale = attr(X[[j]],paste0("sigma:", levels(study.test)[m])))
                }else{#if just one level in object$study, the normalisation attributes are not named the same
                  newdata.list.study.scale.temp=scale(newdata.list.study[[j]][[m]], center = attr(X[[j]],"scaled:center"), scale = attr(X[[j]],"scaled:scale"))
                }
              }
              
              if(scale==FALSE)
              {
                if(nlevels(object$study)>1)
                {
                  newdata.list.study.scale.temp=scale(newdata.list.study[[j]][[m]], center = attr(X[[j]],paste0("means:", levels(study.test)[m])),scale=FALSE)
                }else{#if just one level in object$study, the normalisation attributes are not named the same
                  newdata.list.study.scale.temp=scale(newdata.list.study[[j]][[m]], center = attr(X[[j]],"scaled:center"),scale=FALSE)
                }
                
                
              }
              
            }else{
              # case 2: sample test is a new study # a new study which was not in the learning set
              newdata.list.study.scale.temp=scale(newdata.list.study[[j]][[m]],center=TRUE,scale=scale)
            }
            #concatenation of the scaled data
            concat.newdata[[j]] = rbind(concat.newdata[[j]], unlist(newdata.list.study.scale.temp))#[[j]][[m]]))
          }
        }
        
        # we now need to reorganise concat.newdata as newdata. Indeed, the concatenation was done without taking care of the order of the samples in newdata
        for(j in 1:J) # loop on the blocks
        {
          indice.match=match(rownames(newdata[[j]]),rownames(concat.newdata[[j]]))
          #match indice
          concat.newdata[[j]]=concat.newdata[[j]][indice.match,]
          concat.newdata[[j]][which(is.na(concat.newdata[[j]]))]=0 # taking care of the NA due to normalisation: put to 0 so they don't influence the product below (in Y.hat), reviens au meme que de supprimer les colonnes sans variances au depart.
          
        }
        names(concat.newdata)=names(X)
        
        means.Y=matrix(0,nrow=nrow(concat.newdata[[1]]),ncol=q)
        sigma.Y=matrix(1,nrow=nrow(concat.newdata[[1]]),ncol=q)
        
        #loop on the blocks to define means.Y and sigma.Y for mint analysis
        for(m in 1:M)
        {
          if(m%in%match.study.indice) #if some study of study.test were already in the learning set
          {
            if(nlevels(object$study)>1)
            {
              means.Y[which(study.test%in%levels(study.learn)[match.study[m]]),]=matrix(attr(Y,paste0("means:", levels(study.test)[m])),nrow=length(which(study.test%in%levels(study.learn)[match.study[m]])),ncol=q,byrow=TRUE)
            }else{#if just one level in object$study, the normalisation attributes are not named the same
              means.Y[which(study.test%in%levels(study.learn)[match.study[m]]),]=matrix(attr(Y,"scaled:center"),nrow=length(which(study.test%in%levels(study.learn)[match.study[m]])),ncol=q,byrow=TRUE)
            }
            
            
            if(scale==TRUE)
            {
              # I want to build a vector with sigma.Y for each group
              if(nlevels(object$study)>1)
              {
                sigma.Y[which(study.test%in%levels(study.learn)[match.study[m]]),]=matrix(attr(Y,paste0("sigma:", levels(study.test)[m])),nrow=length(which(study.test%in%levels(study.learn)[match.study[m]])),ncol=q,byrow=TRUE)
              }else{#if just one level in object$study, the normalisation attributes are not named the same
                sigma.Y[which(study.test%in%levels(study.learn)[match.study[m]]),]=matrix(attr(Y,"scaled:scale"),nrow=length(which(study.test%in%levels(study.learn)[match.study[m]])),ncol=q,byrow=TRUE)
              }
              
              
            }
          }
        }
        
      }
      ### end if object is a mint analysis
    } else {
      means.Y = matrix(attr(Y, "scaled:center"),nrow=nrow(newdata[[1]]),ncol=q,byrow=TRUE);
      if (scale)
      {sigma.Y = matrix(attr(Y, "scaled:scale"),nrow=nrow(newdata[[1]]),ncol=q,byrow=TRUE)}else{sigma.Y=matrix(1,nrow=nrow(newdata[[1]]),ncol=q)}
      concat.newdata=newdata
      names(concat.newdata)=names(X)
    }

    
    ### at this stage we have
    # X         # list of blocks
    # Y         # observation
    # newdata   #list of blocks for the prediction, same length as A, scaled
    
    ###  replace NA by 0 in the training data
    if( (hasArg(misdata.all) & hasArg(is.na.X)) && any(list(...)$misdata.all)) # ind.na.X: all blocks except Y
    {    # if misdata.all and ind.na.X are provided, we don't calculate the is.na(X) as it takes time. Used in tune functions.
      for(j in c(1:J)[list(...)$misdata.all])
        X[[j]][list(...)$is.na.X[[j]]]=0 # faster than using replace
    } else {
      # replace NA by 0
      X = lapply(X,function(x)
      {
        if (anyNA(x)){
          ind = is.na(x)
          x[ind] = 0
        }
        x
      })
    }
    
    ###  replace NA by 0 in the test data
    
    if( (hasArg(misdata.all) & hasArg(is.na.newdata)) && any(list(...)$misdata.all)) # ind.na.newdata: all blocks of newdata
    {
      # if misdata.all and ind.na.X are provided, we don't calculate the is.na(X) as it takes time. Used in tune functions.
      concat.newdata = lapply(1:J, function(q){replace(concat.newdata[[q]], list(...)$is.na.newdata[[q]], 0)})
      
    } else {
      # replace NA by 0
      concat.newdata = lapply(concat.newdata,function(x)
      {
        if (anyNA(x)){
          ind = is.na(x)
          x[ind] = 0
        }
        x
      })
    }
    
    if (any(sapply(concat.newdata, anyNA)))
      stop("Some missing values are present in the test data")
    
    
    # replace NA by 0 in Y
    Y[is.na(Y)] = 0
    
    # -----------------------
    #       prediction
    # -----------------------
    
    B.hat = t.pred = Y.hat = list()
    # J=1 for splsda case. 
    for (i in 1 : J)
    {
      Pmat = Cmat = Wmat = NULL
      
      ### Start estimation using formula Y = XW(P'W)C (+ Yr, residuals on Y) See page 136 La regression PLS Theorie et pratique Tenenhaus
      # Estimation matrix W, P and C
      Pmat = crossprod(X[[i]], variatesX[[i]])
      Cmat = crossprod(Y, variatesX[[i]])
      Wmat = loadingsX[[i]]
      
      # Prediction Y.hat, B.hat and t.pred
      Ypred = lapply(1 : ncomp[i], function(x){concat.newdata[[i]] %*% Wmat[, 1:x] %*% solve(t(Pmat[, 1:x]) %*% Wmat[, 1:x]) %*% t(Cmat)[1:x, ]})
      Ypred = sapply(Ypred, function(x){x*sigma.Y + means.Y}, simplify = "array")
      
      Y.hat[[i]] = array(Ypred, c(nrow(newdata[[i]]), ncol(Y), ncomp[i])) # in case one observation and only one Y, we need array() to keep it an array with a third dimension being ncomp
      
      t.pred[[i]] = concat.newdata[[i]] %*% Wmat %*% solve(t(Pmat) %*% Wmat)
      t.pred[[i]] = matrix(data = sapply(1:ncol(t.pred[[i]]),
                                         function(x) {t.pred[[i]][, x] * apply(variatesX[[i]], 2,
                                                                               function(y){(norm(y, type = "2"))^2})[x]}), nrow = nrow(concat.newdata[[i]]), ncol = ncol(t.pred[[i]]))
      
      B.hat[[i]] = sapply(1 : ncomp[i], function(x){Wmat[, 1:x] %*% solve(t(Pmat[, 1:x]) %*% Wmat[, 1:x]) %*% t(Cmat)[1:x, ]}, simplify = "array")
      ### End estimation using formula Y = XW(P'W)C (+ Yr, residuals on Y) See page 136 La regression PLS Theorie et pratique Tenenhaus
      
      rownames(t.pred[[i]]) = rownames(newdata[[i]])
      colnames(t.pred[[i]]) = paste("dim", c(1:ncomp[i]), sep = " ")
      rownames(Y.hat[[i]]) = rownames(newdata[[i]])
      colnames(Y.hat[[i]]) = colnames(Y)
      dimnames(Y.hat[[i]])[[3]]=paste("dim", c(1:ncomp[i]), sep = " ")
      rownames(B.hat[[i]]) = colnames(newdata[[i]])
      colnames(B.hat[[i]]) = colnames(Y)
      dimnames(B.hat[[i]])[[3]]=paste("dim", c(1:ncomp[i]), sep = " ")
      
    }

    #-- valeurs sortantes --#
    names(Y.hat) = names(t.pred) = names(B.hat) = names(object$X)

    out = list(Y.hat=Y.hat, t.pred=t.pred, B.hat=B.hat, 
               # Manually verified that this is the new data transformed into the the new PCA space. 
               # Checked by building an splsda object on the training dataset, and then RE-transforming the same 
               # training data through the splsda object with this function.
               newdata.transformed=t.pred) 
    
    class(out) = paste("predict")
    
    return(out)
  }

# --------------------------------------
# Check.entry.single: check input parameter for a single matrix that comes from a list of matrix
# --------------------------------------
# X: a single numeric matrix
# ncomp: the number of components to include in the model
# q: position of X in a list of matrix, as used in horizontal analysis (e.g. block.pls)
Check.entry.single = function(X,  ncomp, q)
{
  
  #-- validation des arguments --#
  if (length(dim(X)) != 2)
    stop(paste0("'X[[", q, "]]' must be a numeric matrix."))
  
  if(! any(class(X) %in% "matrix"))
    X = as.matrix(X)
  
  if (!is.numeric(X))
    stop(paste0("'X[[", q, "]]'  must be a numeric matrix."))
  
  N = nrow(X)
  P = ncol(X)
  
  if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop(paste0("invalid number of variates 'ncomp' for matrix 'X[[", q, "]]'."))
  
  ncomp = round(ncomp)
  
  # add colnames and rownames if missing
  X.names = dimnames(X)[[2]]
  if (is.null(X.names))
  {
    X.names = paste("X", 1:P, sep = "")
    dimnames(X)[[2]] = X.names
  }
  
  ind.names = dimnames(X)[[1]]
  if (is.null(ind.names))
  {
    ind.names = 1:N
    rownames(X)  = ind.names
  }
  
  if (length(unique(rownames(X))) != nrow(X))
    stop("samples should have a unique identifier/rowname")
  if (length(unique(X.names)) != P)
    stop("Unique indentifier is needed for the columns of X")
  
  return(list(X=X, ncomp=ncomp, X.names=X.names, ind.names=ind.names))
}
