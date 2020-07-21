#' clusterlasso_cv
#'
#' This function carries out ClusterLasso-penalized logistic regression cross-validation using glmnet 
#' @param X [n_sample x n_features] feature matrix 
#' @param y [n_sample] vector containing reference labels (+1/-1) 
#' @param subgroup [n_sample] vector containing a "subgroup" information (e.g., bacterial lineages), to be stratified among the folds (optional)
#' @param n.lambda number of regularization values considered (default 200)
#' @param n.folds number of folds (default 10)
#' @param n.repeat number of repetitions of the cv process (default 3)
#' @param seed random seed (default 27)
#' @param verbose show progress (default TRUE)
#' @param screen.thresh correlation threshold considered to screen variables (default 0.9)
#' @param clust.thresh correlation threshold considered to cluster variables (default 0.9)
#' @param clust.summary function used to summarize clusters (default "mean")
#' @export
#' @examples
#' clusterlasso_cv()
clusterlasso_cv = function(X, y, subgroup = NULL, n.lambda = 200, n.folds = 10, n.repeat = 3, seed = 27, verbose = TRUE, screen.thresh = 0.9, clust.thresh = 0.9, clust.summary = "mean"){
	# check y values belong to {-1,+1}
	if( !identical(sort(unique(y)),c(-1,1)) ){
		stop("reference labels must belong to {-1;1}")
	}
	# set random seed
	if(!is.null(seed)){
		set.seed(seed)
	}

	#----------------------------------------------------------------#
	# fit model on entire dataset to get final model and lambda grid #
	#----------------------------------------------------------------#
	if(verbose){
		cat("**** fitting model on entire dataset ****\n")
	}
		# fit model of individual variables #
		#-----------------------------------#
	if(verbose){
		cat("\t-fitting lasso on individual variables\n")
	}
	global.fit = glmnet(X, y, family = "binomial", nlambda = n.lambda, standardize = FALSE)	# NB : never standardize in glmnet, done beforehand
	lambda.grid.global = global.fit$lambda
		# extract support #
		#-----------------#
	n.feats = rep(0, length(lambda.grid.global))
	supports = list()
	for(i in 1:length(lambda.grid.global)){   
		supports[[i]] = which(abs(global.fit$beta[,i]) > 0)
		n.feats[i] = length(supports[[i]])
	}
		# create data frame
	model_global.stats = list("supports"=supports, "n.feats"=n.feats, "lambda.grid"=lambda.grid.global)   

		# compute clustered dataset #
		#--------------------------#
	if(verbose){
		cat("\t-clustering variables\n")
	}
		# cluster features
	cluster.ids = learnClusters(X = X, y = y, screen.thresh = screen.thresh, clust.thresh = clust.thresh, lambda.grid = lambda.grid.global, verbose = verbose)
		# transform matrices
	Xclust = buildClusters(X, cluster.ids, clust.summary = clust.summary)
		# fit model on clustered variables #
		#----------------------------------#
	if(verbose){
		cat("\t-fitting model on clustered variables\n")
	}
	cluster.fit = glmnet(Xclust, y, family = "binomial", nlambda = n.lambda, standardize = FALSE)	# NB : never standardize in glmnet, done beforehand
	lambda.grid.clust = cluster.fit$lambda
		# extract support #
		#-----------------#
	n.feats = rep(0, length(lambda.grid.clust))
	supports = list()
	for(i in 1:length(lambda.grid.clust)){
		supports[[i]] = which(abs(cluster.fit$beta[,i]) > 0)
		n.feats[i] = length(supports[[i]])
	}
		# create data frame
	model_clust.stats = list("supports"=supports, "n.feats"=n.feats, "lambda.grid"=lambda.grid.clust)  

	#----------------------#
	# store data in a list #
	#----------------------#
	global.model = list("global.fit"=cluster.fit, "model.stats"=model_clust.stats, "cluster.ids"=cluster.ids, "global.fit.individual"=global.fit, "model.stats.individual" = model_global.stats)

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
		res.tmp$preds = matrix(0, nrow = nrow(X), ncol = length(lambda.grid.clust))  
		res.tmp$probs = res.tmp$preds
		#--------------#	
		# define folds #
		#--------------#
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
			# identify clusters
			cluster.ids.fold = learnClusters(X = X[ind.train,], y = y[ind.train], screen.thresh = screen.thresh, clust.thresh = clust.thresh, lambda.grid = lambda.grid.global, verbose = verbose)
			# transform matrices
			Xclust.train = buildClusters(X[ind.train,], cluster.ids.fold, clust.summary = clust.summary)
			Xclust.test  = buildClusters(X[ind.test,], cluster.ids.fold, clust.summary = clust.summary)
			# fit model using appropriate grid of lambda
			fold.fit = glmnet(Xclust.train, y[ind.train], family = "binomial" , lambda = lambda.grid.clust, standardize = FALSE)
			# make predictions
				# extract predictions
			probs = predict(fold.fit, Xclust.test, type = "response")
			preds = predict(fold.fit, Xclust.test, type = "class")
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




#' learnClusters
#'
#' This function identifies cluster involved in the ClusterLasso method
#' @param X [n_sample x n_features] feature matrix 
#' @param y [n_sample] vector containing reference labels (+1/-1) 
#' @param screen.thresh correlation threshold considered to screen variables (default 0.9)
#' @param clust.thresh correlation threshold considered to cluster variables (default 0.9)
#' @param lambda.grid use specified grid of regularization values (default NULL : computed from data)
#' @param n.lambda number of regularization values considered (default 200)
#' @param max.size maximum size to consider to build entire Nactive x Nfeature correlation matrix (default 250000)
#' @param chunk.size chunk size to consider when Nactive x Nfeature correlation matrix is computed by chunk (default 10000)
#' @param verbose show progress (default TRUE)
#' @export
#' @examples
#' learnClusters()
learnClusters = function(X, y, screen.thresh = 0.9, clust.thresh = 0.9, lambda.grid = NULL, n.lambda = 200, max.size = 250000, chunk.size = 10000, verbose = TRUE){

	#-----------------#
	# fit lasso model #
	#-----------------#
	if(is.null(lambda.grid)){
		fit = glmnet(X, y, family = "binomial", nlambda = n.lambda, standardize = FALSE)
	}else{
		fit = glmnet(X, y, family = "binomial", lambda = lambda.grid, standardize = FALSE)
	}
	beta = fit$beta

	#-------------------------#
	# extract active features #
	#-------------------------#
	m = apply(abs(beta), 1, max) 
	ind.active = which(m > 0)
	# extract max coefficients observed across path
	beta.active = beta[ind.active,]
	beta.active.max = apply(abs(beta.active), 1, max)
	# extact feature matrix
	X.active = X[,ind.active]

	#-----------------------------------------------------#
	# identify variables correlated with active variables #
	#-----------------------------------------------------#
	if(verbose){
		cat("\t\t-flagging variables correlated to active\n")
	}
	# if matrix not too large : compute full correlation matrix
	if(ncol(X) < max.size){
		G = sparse_cor(X.active, X)
		G = round(G, digits = 5) # WARNING : ROUND MATRIX FOR NUMERICAL ISSUES (1 != almost 1...)
		# flag variables
		ind.ext = apply(G, 1, function(x){which(x>screen.thresh)})
		ind.ext = unique(unlist(ind.ext))
	}else{ # otherwise : compute by chunk
		chunk.id = seq(ncol(X)) %/% chunk.size
		chunk.id = chunk.id + 1 # start at 1 instead of 0
		ind.ext = c()
		for(i in 1:max(chunk.id)){
			if(verbose){
				cat("\t\t\t- processing chunk no", i, "out of", max(chunk.id), "...")
				start.time = Sys.time()
			}
			ind.chunk = which(chunk.id == i)
			Gtmp = sparse_cor(X.active, X[,ind.chunk])
			Gtmp = round(Gtmp, digits = 5) # WARNING : ROUND MATRIX FOR NUMERICAL ISSUES (1 != almost 1...)
			ind.tmp = apply(Gtmp, 1, function(x){which(x>screen.thresh)})
			ind.tmp = ind.chunk[ unique(unlist(ind.tmp)) ]
			ind.ext = unique(c(ind.ext,ind.tmp))
			if(verbose){
				end.time = Sys.time()
				cat(" took", difftime(end.time, start.time, units = "mins"), "minutes\n")
			}
		}
	}
	# extract maximum beta observed on regularization path (--> for visualisation)
	beta.ext.max = rep(0, length(ind.ext))
	ind.match = match(ind.active, ind.ext)
	beta.ext.max[ind.match] = beta.active.max

	#-----------------------------------#
	# carry out hierarchical clustering #
	#-----------------------------------#
	# compute distance matrix
	G = sparse_cor(X[,ind.ext])
	D = as.dist(1 - abs(G))
	# carry out hierarchical clustering
	hc = hclust(D, method = "single")
	# extract clusters
	clusters = cutree(hc, h = 1-clust.thresh)

	#----------------#
	# return results #
	#----------------#
	res = list("clusters"=clusters, "ind.active" = ind.active, "ind.extended" = ind.ext, "beta.active" = beta.active, "beta.active.max" = beta.active.max, "beta.extended.max"=beta.ext.max)
	return(res)
}


#' buildClusters
#'
#' This function "builds" clusters involved in the ClusterLasso method, i.e., summarizes features of each cluster into a single variable
#' @param X [n_sample x n_features] feature matrix 
#' @param cluster.ids cluster ids returned by function learnClusters (entire list returned)
#' @param clust.summary function used to summarize clusters (default "mean")
#' @export
#' @examples
#' buildClusters()
buildClusters = function(X, cluster.ids, clust.summary = "mean"){
	# get summarizing function
	clust.sum = match.fun(clust.summary) 
	# extract cluster ids
	clust = cluster.ids$clusters
	n.clust = max(clust)
	ind.by.clust = split(seq(length(clust)), factor(clust))
	# merge variables (specific treatment if made of a single instance)
	if(is.null(nrow(X))){
		X.ext = X[cluster.ids$ind.extended]
		X.clust = sapply(ind.by.clust, function(x){ if(length(x)==1){tt = X.ext[x]}else{ tt = clust.sum(X.ext[x])}; return(tt)})
	}else{
		X.ext = X[,cluster.ids$ind.extended]
		X.clust = sapply(ind.by.clust, function(x){ if(length(x)==1){tt = X.ext[,x]}else{ tt = apply(X.ext[,x], 1, clust.sum)}; return(tt)})		
	}

	# return results
	return(X.clust)
}

