
#' compute_perf
#'
#' This function computes several performance indicators for binary classification 
#' @param preds [n_sample x n_models] matrix containing binary (-1/+1) predictions 
#' @param probs [n_sample x n_models] matrix containing posterior probability p(y=1|x) (more generally : score associated to positive class) 
#' @param y reference labels
#' @export
#' @examples
#' compute_perf()
compute_perf = function(preds, probs, y){
	# check y values belong to {-1,+1}
	if( !identical(sort(unique(y)),c(-1,1)) ){
		stop("reference labels must belong to {-1;1}")
	}

	# format if prediction from a single run
	if(is.null(nrow(preds))){
		preds = matrix(preds, ncol = 1)
		probs = matrix(probs, ncol = 1)
	}

	#-------------------#
	# intialize results #
	#-------------------#
	res = list()
	#------------------------------------------------------------------#
	# compute performance indicators at default (prediction) threshold #
	#------------------------------------------------------------------#
	# accuracy, sensi, speci, balanced-accuracy
	accuracy = round(100 * apply(preds, 2, function(x){mean(x==y)}), digits = 1)
	sensi = round(100 * apply(preds, 2, function(x){mean(x[y==1] == 1)}), digits = 1)
	speci = round(100 * apply(preds, 2, function(x){mean(x[y==-1] == -1)}), digits = 1)
	balanced.accuracy = round(0.5*(sensi+speci), digits = 1)
	# precision, recall, f1
	rec = sensi
	prec = round(100 * apply(preds, 2, function(x){mean( y[x==1] == 1)}), digits = 1)
	prec[is.nan(prec)] = 0 # fix NaNs (occurs when all predictions are negative)
	f1 = round(2*(rec*prec)/(rec+prec), digits = 1)
	f1[is.nan(f1)] = 0 # fix NaNs (occurs when all predictions are negative)
	# auc
	auc = apply(probs, 2, function(x){
			    p = prediction(x, y, label.ordering = c(-1,1))
			    p2 = performance(p, measure = "auc")
			    return(p2@y.values[[1]])
			})
	auc = round(100*auc, digits = 1)
	# merge
	res$perf = data.frame("accuracy"=accuracy, "sensi"=sensi, "speci"=speci, "balanced.accuracy"=balanced.accuracy, "auc"=auc, "recall"=rec, "precision"=prec, "f1"=f1)
	#---------------------------------------#
	# compute roc & precision/recall curves #
	#---------------------------------------#
	roc = apply(probs, 2, function(x){
			    p = prediction(x, y, label.ordering = c(-1,1))
			    p2 = performance(p, measure = "tpr", x.measure = "fpr")
			    return(p2)
			})
	pr = apply(probs, 2, function(x){
			    p = prediction(x, y, label.ordering = c(-1,1))
			    p2 = performance(p, measure = "rec", x.measure = "prec")
		    return(p2)
		})
	#----------------#
	# return results #
	#----------------#
	res$roc.curves = roc
	res$pr.curves = pr
	return(res)
}


#' stratified_cv_folds
#'
#' This function defines stratified cross-validation folds
#' @param y reference labels
#' @param n.folds number of folds
#' @export
#' @examples
#' stratified_cv_folds()
stratified_cv_folds = function(y, n.folds){
	# split indices by class
	ind.by.class = split(c(1:length(y)), y)
	# make a permutation
		# NB : following line does not work when 1 instance per class (can be the case for pop. structure)
		#ind.by.class = lapply(ind.by.class, sample)
	for(i in 1:length(ind.by.class)){
		if(length(ind.by.class[[i]])  > 1){
			ind.by.class[[i]] = sample(ind.by.class[[i]])
		}
	}
	# initialize cv fold vector
	cv.fold = rep(0, length(y))
	# cut into folds
	for(i in 1:length(ind.by.class)){
		# get indices
		ind = ind.by.class[[i]]
		# cut into n.folds
		# get majority part
		n.by.fold = rep(floor(length(ind)/n.folds), n.folds)
		# distribute modulo part
		n.remain = length(ind) %% n.folds
		if(n.remain > 0){
			tmp = sample(c(1:n.folds), n.remain)
			n.by.fold[tmp] = n.by.fold[tmp] + 1
		}
		# build fold vector
		fold = rep(c(1:n.folds), times = n.by.fold)
		# store in cv.fold variable
		cv.fold[ind] = fold
	}
	# return
	return(cv.fold)
}


#' sparse_cor 
#'
#' This function computes a correlation matrix between the columns of a sparse matrix, or of a pair of sparse matrices. 
#' @param x input matrix
#' @param y second input matrix (optional)
#' @param return.cov boolean flag to return covariance matrix in addition to correlation matrix (default : FALSE)
#' @export
#' @examples
#' sparse_cor()
sparse_cor <- function(x, y = NULL, return.cov = FALSE){
	# NB : code taken and adapted from https://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r
	if(is.null(y)){
	   	n <- nrow(x)
   		cMeans <- Matrix::colMeans(x)
   	 	covmat <- (as.matrix(Matrix::crossprod(x)) - n*Matrix::tcrossprod(cMeans))/(n-1)
   		sdvec <- sqrt(diag(covmat)) 
   	 	cormat <- covmat/Matrix::tcrossprod(sdvec)
   	}else{
   		n <- nrow(x)
    	n2 = nrow(y)
    	if(n != n2){ 
    		stop("matrix dimensions don't match")
    	}	
	    cMeans_x <- Matrix::colMeans(x)
   	 	cMeans_y <- Matrix::colMeans(y)
	    cSd_x = apply(x, 2, sd)
   	 	cSd_y = apply(y, 2, sd)
		covmat <- (  as.matrix(Matrix::crossprod(x,y)) - n*Matrix::tcrossprod(cMeans_x,cMeans_y)  )/(n-1)
    	cormat <- covmat/Matrix::tcrossprod(cSd_x, cSd_y)
   	}
    if(return.cov){
	    return( list(cor=cormat,cov=covmat) )
	 }else{
	    return(cormat)
	}
}


#' get_perf_best_cutoff 
#'
#' This function extract the predictive performance (sensi/speci/accuracy/macro-accuracy) at optimal cutpoints on ROC curves obtained by cross-validation. 
#' @param cv.res cross-validation results (list returned by lasso or cluster-lasso cross-validation fuctions)
#' @export
#' @examples
#' get_perf_best_cutoff()
get_perf_best_cutoff <- function(cv.res){
	# process each repetition of the cross-validation result object
	for(i in 1:length(cv.res$cv.res)){
		# initialize results
		cv.res$cv.res[[i]]$perf$perf$sensi.best = 0
		cv.res$cv.res[[i]]$perf$perf$speci.best = 0
		cv.res$cv.res[[i]]$perf$perf$accuracy.best = 0
		cv.res$cv.res[[i]]$perf$perf$balanced.accuracy.best = 0
		cv.res$cv.res[[i]]$perf$perf$thresh.best = 0
		# extract number of strains
		nr = sum(cv.res$cv.res[[i]]$y.ref == 1)			
		ns = sum(cv.res$cv.res[[i]]$y.ref == -1)
		n = length(cv.res$cv.res[[i]]$y.ref)
		# process each roc curve
		for(j in 1:length(cv.res$cv.res[[i]]$perf$roc.curves)){
			# extract optimal cutoff
			rc = cv.res$cv.res[[i]]$perf$roc.curves[[j]]
			d = rc@x.values[[1]]^2 + (1-rc@y.values[[1]])^2
			ind = floor(mean(which(d==min(d))))
			# compute indicators
			se =  rc@y.values[[1]][ind]
			sp = 1 - rc@x.values[[1]][ind]
			thresh = mean(rc@alpha.values[[1]][ind],rc@alpha.values[[1]][ind+1])
			# store 
			cv.res$cv.res[[i]]$perf$perf$sensi.best[j] = round(100*se, digits = 1)
			cv.res$cv.res[[i]]$perf$perf$speci.best[j] = round(100*sp, digits = 1)
			cv.res$cv.res[[i]]$perf$perf$balanced.accuracy.best[j] = round(100*0.5*(sp+se), digits = 1)			
			cv.res$cv.res[[i]]$perf$perf$accuracy.best[j] = round( 100*(se*nr + sp*ns)/n, digits = 1)		
			cv.res$cv.res[[i]]$perf$perf$thresh.best[j] = thresh
		}
	}
	# return results
	return(cv.res)
}

