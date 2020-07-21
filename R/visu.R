

#' show_cv_overall
#'
#' This function generates a figure (made of 4 individual plots) summarizing overall cross-validation results
#' @param cv.res list containing cross-validation results (produced by lasso_cv() or clusterlasso_cv() ) 
#' @param plot.id string containing a label to consider for the plots (typically, the drug corresponding to the cv.res object)  
#' @param modsel.criterion prediction performance criterion to select best model (balanced.accuracy.best or auc) (defaut = balanced.accuracy.best)
#' @param best.eps AUC tolerance considered to flag best model (defaut = 1)
#' @export
#' @examples
#' show_cv_overall()

show_cv_overall = function(cv.res, plot.id = "cv-results", modsel.criterion = "balanced.accuracy.best", best.eps = 1){

    # check model selection crition is accuracy.macro.best or auc #
    #-------------------------------------------------------------#
    if( !(modsel.criterion %in% c("balanced.accuracy.best","auc")) ){
        stop("can only consider balanced.accuracy.best or auc as model selection criterion")
    }
    # specify list of performance indicators #
    #----------------------------------------#
    perf.crit = c("sensi.best","speci.best","accuracy.best", "balanced.accuracy.best","auc")

    # define colors for performance indicators #
    #------------------------------------------#
    tmp = brewer.pal(10, "Paired")
    cols.perf = c(tmp[c(5,6,1,2,4)])    # sensi, speci, accuracy, balanced.accuracy, auc
    names(cols.perf) = perf.crit

    # extract average perf across cv repetitions #
    #--------------------------------------------#
    p.mu = c()
    p.sd = c()
    for(perf in perf.crit){
        tt = sapply(cv.res$cv.res, function(x){x$perf$perf[,perf]}) 
        p.mu = cbind(p.mu, apply(tt, 1, mean))
        p.sd = cbind(p.sd, apply(tt, 1, sd))
    }
    colnames(p.mu) = perf.crit
    p.mu = round(p.mu, digits = 1)

    # flag best model #
    #-----------------#
    ind.best = which(p.mu[,modsel.criterion] >= (max(p.mu[,modsel.criterion])-best.eps) )
    ind.best = ind.best[1]

    # plot 1 :  perf vs lambda #
    #--------------------------#
    plot(p.mu[,1], type = "l", col = cols.perf[1], lwd = 2, ylim = range(p.mu), xlab = "lambda index", ylab = "performance", main = paste(plot.id, ": average performance\nacross CV repetitions"))
    grid()
    for(i in 1:ncol(p.mu)){
        lines(p.mu[,i], type = "l", lwd = 2, col = cols.perf[i])
    }
    # show best model
    abline(v = ind.best, col = "orange", lty = 2)
    # add legend
    legend("bottomright", colnames(p.mu), col = cols.perf, lwd = 2, bg = "white")

    # plot 2 :  support size vs lambda #
    #----------------------------------#
    n.feats = cv.res$global.model$model.stats$n.feats
    plot(n.feats, type = "l", lwd = 2, xlab = "lambda index", ylab = "number of features", main = paste(plot.id, ": support size"))
    grid()
    # show best model
    abline(v = ind.best, col = "orange", lty = 2)
    mtext(n.feats[ind.best], side = 3, at = ind.best, col = "orange")
   
    # plot 3 : support size vs delta(modsel.criterion) #
    #--------------------------------------------------#
    # define tolerances to consider
    best.eps.list = seq(0, 2, by = 0.25)
    # initialize variables
    n.feats.delta = c()
    # process each tolerance value
    for(b in best.eps.list){
        ind = which(p.mu[,modsel.criterion] >= (max(p.mu[,modsel.criterion])-b) )
        ind = ind[1]
        n.feats.delta = c(n.feats.delta, n.feats[ind])
    }
    # define color
    cols.delta = rep("grey", length(n.feats.delta))
    cols.delta[which(best.eps.list == best.eps)] = "orange"
    # plot
    bargraph = barplot(n.feats.delta, ylim = c(0,1.2*max(n.feats.delta)), col = cols.delta, names = best.eps.list, xlab = paste0("delta(",modsel.criterion,")"), ylab = "number of features", main = paste(plot.id, ": support size\nvs performance tolerance -", modsel.criterion))
    text(bargraph, n.feats.delta, n.feats.delta, pos = 3, col = cols.delta)

}



#' show_cv_best
#'
#' This function generates a figure (made of 4 individual plots) summarizing the best model identified by cross-validation
#' @param cv.res list containing cross-validation results (produced by lasso_cv() or clusterlasso_cv() ) 
#' @param plot.id string containing a label to consider for the plots (typically, the drug corresponding to the cv.res object)  
#' @param modsel.criterion prediction performance criterion to select best model (balanced.accuracy.best or auc) (defaut = balanced.accuracy.best)
#' @param best.eps best-macro-accuracy tolerance considered to flag best model (default = 1)
#' @param method string defining the method considered ("lasso" or "clusterlasso" - default = "lasso")
#' @export
#' @examples
#' show_cv_best()
show_cv_best = function(cv.res, plot.id = "cv-results", modsel.criterion = "balanced.accuracy.best", best.eps = 1, method = "lasso"){

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
    
    # specify list of performance indicators #
    #----------------------------------------#
    perf.crit = c("sensi.best","speci.best","accuracy.best", "balanced.accuracy.best","auc")

    # define colors for performance indicators #
    #------------------------------------------#
    tmp = brewer.pal(10, "Paired")
    cols.perf = c(tmp[c(5,6,1,2,4)])    # sensi, speci, accuracy, accuracy.macro, auc
    names(cols.perf) = perf.crit

    # extract average perf across cv repetitions #
    #--------------------------------------------#
    p.mu = c()
    p.sd = c()
    for(perf in perf.crit){
        tt = sapply(cv.res$cv.res, function(x){x$perf$perf[,perf]}) 
        p.mu = cbind(p.mu, apply(tt, 1, mean))
        p.sd = cbind(p.sd, apply(tt, 1, sd))
    }
    colnames(p.mu) = perf.crit
    p.mu = round(p.mu, digits = 1)

    # extract support size #
    #----------------------#
    n.feats = cv.res$global.model$model.stats$n.feats  

    # flag best model : smallest support up to performance tolerance #
    #----------------------------------------------------------------#
    ind.best = which(p.mu[,modsel.criterion] >= (max(p.mu[,modsel.criterion])-best.eps) )
    ind.best = ind.best[1]

    # extract decision thresholds #
    #-----------------------------#
    thresh.best = sapply(cv.res$cv.res, function(x){x$perf$perf$thresh.best})
    thresh.best = round(thresh.best[ind.best,], digits = 3)

    # plot 1 : overall performance of best model #
    #--------------------------------------------#
    # merge perf and support
    perf.best = c(p.mu[ind.best,], n.feats[ind.best])
    names(perf.best) = c(perf.crit, "support")
    # extract values to plot
    x = perf.best
    y = c(p.sd[ind.best,], 0)
    # format
    names(x) = gsub(".best", "", names(x))
    names(x)[names(x)=="balanced.accuracy"] = "balanced\naccuracy"
    # plot
    bargraph = barplot(x, ylim = c(0,110), xpd = F, col = c(cols.perf, "grey"), space = c(rep(0.2,length(x)-2), 1, 0.2), las = 2, main = paste(plot.id, ": best model performance"))
    arrows(x0 = bargraph, y0 = x-y, y1 = x+y, angle = 90, length = 0.05, code = 3)
    text(bargraph, x, round(x,digits=1), pos = 3)
    if(x[6] > 100){
        text(bargraph[6], 100, x[6], pos = 1)
    }

    # plot 2 : ROC curve #
    #--------------------#
    # define colors
    n.repeat = length(cv.res$cv.res)
    cols.roc = brewer.pal(n.repeat, "Set2")
    # plot roc curves obtained across cv repetitions
    plot(cv.res$cv.res[[1]]$perf$roc.curves[[ind.best]], col = cols.roc[1], lwd = 2, main = paste(plot.id, ": best model ROC curve"))
    abline(0, 1, lty = 2)
    grid()
    for(j in 1:length(cv.res$cv.res)){
        plot(cv.res$cv.res[[j]]$perf$roc.curves[[ind.best]], add = TRUE, col = cols.roc[j], lwd = 2)
    }
    legend("bottomright", paste("cv repeat no", seq(n.repeat)), col = cols.roc, lwd = 2, bg = "white")
    # add (averarge) "best" sensi/speci
    se = perf.best["sensi.best"]/100
    sp = perf.best["speci.best"]/100
    points(1-sp, se, pch = 19)

    ## add best roc curve ?
    #ind.best.auc = which.max(p.mu[,"auc"])
    #plot(cv.res$cv.res[[1]]$perf$roc.curves[[ind.best.auc]], add = TRUE, col = "black", lty = 2)

    # plot 3 : distribution of betas #
    #--------------------------------#
    # extract coefficients
    beta = cv.res$global.model$global.fit$beta[,ind.best]
    # extract active unitigs
    ind.active = which(abs(beta)>0)
    beta = beta[ind.active]
    # order by effect
    ind.sort = order(abs(beta), decreasing = T)
    beta = beta[ind.sort]
    # plot
    bargraph = barplot(abs(beta), names = seq(length(beta)), ylim = c(0, 1.2*max(abs(beta))), main = paste(plot.id, ": absolute values of model coefficients"))

    # plot 4 : number of patterns per cluster #
    #-----------------------------------------#
    if(method == "clusterlasso"){
        # extract members of the clusters
        cluster.ids = cv.res$global.model$cluster.ids$clusters                  # a n.extended vector containing cluster membership
        ind.extended.features = cv.res$global.model$cluster.ids$ind.extended    # a n.extented vector containing indices in the original feature matrix
        patterns.cluster = list()
        for(i in 1:length(beta)){
            ind = which(cluster.ids ==  as.numeric( names(beta)[i]) )
            patterns.cluster[[i]] = ind.extended.features[ind]
        }
        # plot 
        def_par=par()$mar
        par(mar=c(2,4,4,2))
        n = sapply(patterns.cluster, length)
        bargraph = barplot(n, names = seq(length(n)), ylim = c(0, 1.2*max(n)), main = paste(plot.id, ": size of clusters"))
        text(bargraph, n, n, pos = 3)
        par(mar=def_par)
    }

    # extract performance of best model #
    #-----------------------------------#
    # add decision thresholds
    t = paste(thresh.best, collapse = "-")
    t.mean = round(mean(thresh.best), digits = 3)
    # store
    perf.best = c(perf.best, t.mean, t)
    names(perf.best)[length(perf.best)-1] = "thresh-best_mean"
    names(perf.best)[length(perf.best)] = "thresh-best_all"
    # return
    return(perf.best)
}


#' heatmap_correlation_signatures
#'
#' This function generates a heatmap figure comparing the signatures of a lasso and a cluster-lasso model (returned by the "extract_best_model" function)
#' @param model.lasso lasso model (returned by "extract_best_model" function)
#' @param model.cluster lasso model (returned by "extract_best_model" function)
#' @param clust.min minimum size of clusters to consider in the representation (default : 10)
#' @param plot.title title for the plot (default = "" - empty)
#' @export
#' @examples
#' heatmap_correlation_signatures()
heatmap_correlation_signatures = function(X, model.lasso, model.cluster, clust.min = 10, plot.title = ""){

    # extract features involved in lasso signature #
    #----------------------------------------------#
    ind.lasso = model.lasso$indices

    # extract features involved in cluster lasso signatures #
    #-------------------------------------------------------#
    clust.ids = model.cluster$indices
    ind.active = which(model.cluster$cluster.ids$clusters %in% clust.ids)
    # extract indices in original matrix and cluster ids
    ids.clust = factor(model.cluster$cluster.ids$clusters[ind.active]) 
    ind.clust = model.cluster$cluster.ids$ind.extended[ind.active]
    # format beta - add a value per feature
    ind.match = match(ids.clust, model.cluster$indices)
    beta.clust =  model.cluster$beta[ind.match]

    # extract common list of features #
    #---------------------------------#
    ind.all = unique(c(ind.lasso, ind.clust))

    # compute correlation matrix #
    #----------------------------#
    G = cor(as.matrix(X[,ind.all]))

    # build annotation data frame #
    #-----------------------------#
    # flag patterns of both signatures
    lasso.sig = ind.all %in% ind.lasso
    cluster.sig = ind.all %in% ind.clust

    # identify clusters 
    cluster.ids = rep(0, length(ind.all))
    ind.match = match(ind.clust, ind.all)
    cluster.ids[ind.match] = ids.clust

    # add lasso model coefficients
    beta.lasso = rep(0, length(ind.all))
    ind.match = match(ind.lasso, ind.all)
    beta.lasso[ind.match] = abs(model.lasso$beta)

    # add cl-lasso model coefficients
    beta.cluster = rep(0, length(ind.all))
    ind.match = match(ind.clust, ind.all)
    beta.cluster[ind.match] = abs(beta.clust)

    # format cluster ids 
        # 1) dicard "empty" : features from the lasso signature not part of the cluster one
    ind.noclust = which(cluster.ids == 0)
    if(length(ind.noclust) > 0){
     cluster.ids[ind.noclust] = NA
    }
        # 2) focus on largest clusters
    tt = table(cluster.ids)
    clust.maj = as.numeric(names(tt)[which(tt >= clust.min)])

    cluster.ids.fmt = rep(NA, length(cluster.ids))
    ind.maj = which(cluster.ids %in% clust.maj)
    cluster.ids.fmt[ind.maj] = cluster.ids[ind.maj]

    cluster.ids.fmt = factor(cluster.ids.fmt)
    levels(cluster.ids.fmt) = seq(nlevels(cluster.ids.fmt))

    # define annotation data.frame for the aheatmap function
    annot = data.frame("lasso.sig"=factor(lasso.sig, levels = c("FALSE","TRUE")), "clusterLasso.sig"=factor(cluster.sig,levels = c("FALSE","TRUE")), "cluster.id"=factor(cluster.ids.fmt), "beta.lasso"=beta.lasso, "beta.cluster"=beta.cluster)

    # define colors for clusters
    myPal = colorRampPalette(brewer.pal(9, "Set1")) 
    cols.clust = myPal(nlevels(annot$cluster.id))
    # define all colors
    annot.cols = list("lasso.sig" = c("white","orange"), "clusterLasso.sig"=c("white","deepskyblue"), "cluster.id"=cols.clust)
    # plot
    aheatmap(G, distfun = as.dist(1-G), annCol = annot, annColors = annot.cols, main = plot.title)
}



