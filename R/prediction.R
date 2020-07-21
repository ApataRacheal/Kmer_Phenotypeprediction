
#' extract_best_model
#'
#' This function extract the best lasso model from cross-validation results
#' @param cv.res list containing cross-validation results (produced by lasso_cv() or clusterlasso_cv() ) 
#' @param modsel.criterion prediction performance criterion to select best model (balanced.accuracy.best or auc) (defaut = balanced.accuracy.best)
#' @param best.eps AUC tolerance considered to flag best model (defaut = 1)
#' @param method string defining the method considered ("lasso" or "clusterlasso" - default = "lasso")
#' @export
#' @examples
#' extract_best_model()
extract_best_model = function(cv.res, modsel.criterion = "balanced.accuracy.best", best.eps = 1, method = "lasso"){

    # check model selection crition is accuracy.macro.best or auc #
    #-------------------------------------------------------------#
    if( !(modsel.criterion %in% c("balanced.accuracy.best","auc")) ){
        stop("can only consider balanced.accuracy.best or auc as model selection criterion")
    }

    # check method name is valid #
    #----------------------------#
    if( !(method %in% c("lasso","clusterlasso")) ){
        stop("method can only be defined as lasso or clusterlasso")
    }

    # extract average perf across cv repetitions #
    #--------------------------------------------#
	tt = sapply(cv.res$cv.res, function(x){x$perf$perf[,modsel.criterion]}) 
	p.mu = apply(tt, 1, mean)
	p.sd = apply(tt, 1, sd)

    # flag best model #
    #-----------------#
    ind.best = which(p.mu >= (max(p.mu)-best.eps) )
    ind.best = ind.best[1]

    # extract model coefficients #
    #----------------------------#
    # extract coefficients
    beta = cv.res$global.model$global.fit$beta[,ind.best]
    ind.active = which(abs(beta)>0)
    beta = beta[ind.active]
    # extract intercept
    intercept = cv.res$global.model$global.fit$a0[ind.best]

    # extract "best" decision thresholds #
    #------------------------------------#
    thresh.best = sapply(cv.res$cv.res, function(x){x$perf$perf$thresh.best})
    thresh = round( mean(thresh.best[ind.best,]), digits = 3)

    # create results #
    #----------------#
    res = list("indices" = ind.active, "beta" = beta, "intercept" = intercept, "threshold" = thresh)

    # if clusterLasso : extract cluster ids #
    #---------------------------------------#
    if(method == "clusterlasso"){
        res[["cluster.ids"]] = cv.res$global.model$cluster.ids
    }

    # return results #
    #----------------#
    return(res)
}

#' predict_clusterlasso
#'
#' This function makes predictions from the model selected by the extract_best_model function
#' @param X n.samples x n.features matrix 
#' @param model (best) lasso / clusterlasso model
#' @param method string defining the method considered ("lasso" or "clusterlasso" - default = "lasso")
#' @param clust.summary function used to summarize clusters in clusterlasso (default "mean")
#' @export
#' @examples
#' predict_clusterlasso()
predict_clustlasso = function(X, model, method = "lasso", clust.summary = "mean"){

    # check method name is valid #
    #----------------------------#
    if( !(method %in% c("lasso","clusterlasso")) ){
        stop("method can only be defined as lasso or clusterlasso")
    }

    # build clusters #
    #----------------#
    if(method == "clusterlasso"){
        Xclust = buildClusters(X, model$cluster.ids, clust.summary = clust.summary)
    }else{
        Xclust = X
    }

	# extract model features #
    #------------------------#
    # check X made of a single instance or not
    if(is.null(nrow(Xclust))){
    	Xs = Xclust[model$indices]
        Xs = matrix(Xs, nrow = 1) 
    }else{
        Xs = Xclust[,model$indices] 
    }

	# compute score #
    #---------------#
	score = apply(Xs, 1, function(x){ sum(x*model$beta)})
	score = score + model$intercept

	# compute proba #
    #---------------#
	probs = exp(score)/(1+exp(score))

	# make predictions #
    #------------------#
	preds = rep(-1, nrow(Xs))
	preds[probs >= model$thresh] = 1

	# return results #
    #----------------#
	res = list("preds"=preds, "probs" = probs)
	return(res) 
}
