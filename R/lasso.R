
#' lasso_cv
#'
#' This function carries out Lasso-penalized logistic regression cross-validation using glmnet 
#' @param X [n_sample x n_features] feature matrix 
#' @param y [n_sample] vector containing reference labels (+1/-1) 
#' @param subgroup [n_sample] vector containing a "subgroup" information (e.g., bacterial lineages), to be stratified among the folds (optional)
#' @param n.lambda number of regularization values considered (default 200)
#' @param n.folds number of folds (default 10)
#' @param n.repeat number of repetitions of the cv process (default 3)
#' @param seed random seed (default 27)
#' @param verbose show progress (default TRUE)
#' @param alpha elasticnet parameter (default alpha=1, i.e., lasso)
#' @export
#' @examples
#' lasso_cv()
lasso_cv = function(X, y, subgroup = NULL, n.lambda = 200, n.folds = 10, n.repeat = 3, seed = 27, verbose = TRUE, alpha = 1){

	# check y values belong to {-1,+1}
	if( !identical(sort(unique(y)),c(-1,1)) ){
		stop("reference labels must belong to {-1;1}")
	}
	# set random seed
	if(!is.null(seed)){
		set.seed(seed)
	}
	#-----------------------------------------------#
	# fit global model to get grid of lambda to use #
	#-----------------------------------------------#
	if(verbose){
		cat("**** fitting model on entire dataset ****\n")
	}
		# fit model #
		#-----------#
	global.fit = glmnet(X, y, family = "binomial", nlambda = n.lambda, standardize = FALSE, alpha = alpha)	# NB : never standardize in glmnet, done beforehand
	lambda.grid = global.fit$lambda
		# extract support #
		#-----------------#
	n.feats = rep(0, length(lambda.grid))
	supports = list()
	for(i in 1:length(lambda.grid)){
		supports[[i]] = which(abs(global.fit$beta[,i]) > 0)
		n.feats[i] = length(supports[[i]])
	}
		# create data frame
	model.stats = list("supports"=supports, "n.feats"=n.feats, "lambda.grid"=lambda.grid)

	#----------------------------#
	# create global model object #
	#----------------------------#
	global.model = list("global.fit"=global.fit, "model.stats" = model.stats)

	#----------------------#
	# run cross-validation #
	#----------------------#
	cv.res = list()
	# carry out repetitions
	for(j in 1:n.repeat){
		if(verbose){
			cat("*** processing repeat", j, "out of", n.repeat, "***\n")
		}
		#--------------------#
		# initialize results #
		#--------------------#
		res.tmp = list()
		res.tmp$preds = matrix(0, nrow = nrow(X), ncol = length(lambda.grid))
		res.tmp$probs = res.tmp$preds
		# define folds
		if(is.null(subgroup)){
			cv.folds = stratified_cv_folds(y, n.folds)
		}else{
			subgroup_pheno = factor( paste(as.character(y), as.character(subgroup), sep = "_") )
			cv.folds = stratified_cv_folds(subgroup_pheno, n.folds)
		}
		# store
		res.tmp$cv.folds = cv.folds
		res.tmp$y.ref = y
		#-------------------#
		# process each fold #
		#-------------------#
		for(i in 1:n.folds){
			if(verbose){
				cat("\t- processing fold", i, "out of", n.folds, "\n")
			}
			ind.train = which(cv.folds != i)
			ind.test = which(cv.folds == i)
			# fit model using appropriate grid of lambda
			fold.fit = glmnet(X[ind.train,], y[ind.train], family = "binomial" , lambda = lambda.grid, standardize = FALSE, alpha = alpha)
			# make predictions
				# extract predictions
			probs = predict(fold.fit, X[ind.test,], type = "response")
			preds = predict(fold.fit, X[ind.test,], type = "class")
				# check dimensions match (numerical issues of lasso)
			if(ncol(probs) == ncol(res.tmp$probs)){
				res.tmp$probs[ind.test,] = probs
				res.tmp$preds[ind.test,] = preds
			}else{
				# if numerical issue : fix matrix
				probs.fix = matrix(0, nrow = nrow(probs), ncol = ncol(res.tmp$probs))
				preds.fix = probs.fix
				# store first columns
				probs.fix[,1:ncol(probs)] = probs
				preds.fix[,1:ncol(preds)] = preds
				# fix missing columns with last column obtained (smallest lambda)
				for(k in seq(ncol(probs)+1,ncol(probs.fix)) ){
					probs.fix[,k] = probs[,ncol(probs)]
					preds.fix[,k] = preds[,ncol(probs)]
				}
				res.tmp$probs[ind.test,] = probs.fix
				res.tmp$preds[ind.test,] = preds.fix
			}
		}
		#-----------------------------------------------#
		# compute performance for binary classification #
		#-----------------------------------------------#
		res.tmp$perf = compute_perf(res.tmp$preds, res.tmp$probs, res.tmp$y.ref)
		#-------#
		# store #
		#-------#
		cv.res[[j]] = res.tmp
	}

	#----------------#
	# format results #
	#----------------#
	res = list("cv.res" = cv.res, "global.model" = global.model)

	#------------------------------------------------------#
	# compute perf at optimal cut-points on the ROC curve  #
	#------------------------------------------------------#
	res = get_perf_best_cutoff(res)

	#----------------#
	# return results #
	#----------------#
	return(res)
}

