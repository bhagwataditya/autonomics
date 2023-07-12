#' Plot heatmap
#' @param object SummarizedExperiment
#' @param scale_features TRUE or FALSE: whether to scale (i.e. z-score) features
#' @examples 
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, fit = 'limma', subgroupvar = 'Time', block = 'Subject')
#' object %<>% cluster()
feature_sample_heatmap <- function(object){
   
   # Reverse clutser information to have first item on top 
   features_clustered <- all(c('cluster', 'cluster_order') %in% fvars(object))
   if (features_clustered){
      object %<>% extract(rev(fdata(object)$cluster_order), )
      fdata(object)$cluster %<>% factor(rev(levels(.)))
   }
   
   # Z-score features   
   values(object) %<>% t() %>% scale() %>% t()
   
   # Plot image
   fields::image.plot(
      x = seq(1,ncol(object)), 
      y = seq(1,nrow(object)),
      z = t(values(object)), 
      col = hcl.colors(12, "YlOrRd", rev = TRUE),
      xaxt= "n", yaxt= "n", xlab = "", ylab = "")
   axis( 1, at=seq(1,length.out=ncol(object) ), labels= colnames(object), las= 2, cex.axis = 0.7)
   axis( 2, at=seq(1,length.out=nrow(object) ), labels= rownames(object), las= 2, cex.axis = 0.7)
   
   # Add subgroup lines. Add cluster lines
   abline(v = cumsum(table(object[[subgroupvar]]))+0.5)
   #abline(v = cumsum(table(split_subgroup_values(object)$V1))+0.5, lwd = 3)
   
   if (features_clustered){
      abline(h = cumsum(table(fdata(object)$cluster))+0.5)
   }
   
}


fsplit <- function(object, fvar){
   Map(function(x) object[fdata(object)[[fvar]] == x, ], 
       flevels(object, fvar))
}

plot_cluster_contrastograms <- function(object, subgroupvar){
   n <- length(flevels(object, 'cluster'))
   par(mfrow = c(ceiling(sqrt(n)), floor(sqrt(n))),     # 2x2 layout
       oma = c(0, 0, 0, 0), # bottom, left, top, right space
       mar = c(2, 2, 2, 2), # space for one row of text at ticks and to separate plots
       mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
       xpd = NA)
   lapply(fsplit(object, 'cluster'), plot_contrastogram, subgroupvar=subgroupvar)
}

#' Cluster features
#' @examples 
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, fit = 'limma', subgroupvar = 'Time', block = '')
cluster <- function(object, formula = ~ Subject + Time){

# Fit
    for (col in all.vars(formula)){
        if (is.character(object[[col]])) object[[col]] %<>% factor()
    }
    factors <- vapply(sdata(object)[, all.vars(formula), drop=FALSE], is.factor, logical(1))
    factors <- names(factors)[factors]
    contrast.arg <- rep('contr.sum', length(factors))
    names(contrast.arg) <- factors
    mode(contrast.arg) <- 'list'
    contrastdefs <- create_design(object, formula = formula)
    fit_limma(object, formula=formula, contrastdefs = contrastdefs)

# Filter
    idx <- rowAnys(metadata(object )[[fit]][,,'p'] < filter_p, na.rm=TRUE)
    idx[is.na(idx)] <- FALSE
    object %<>% extract(idx, )
    
# t
    mat <- metadata(object)[[fit]][rownames(object), , quantity]
    #dt <- sumexp_to_longdt(object, fvars='feature_id', svars = c('sample_id', subgroupvar))
    #dt %<>% extract(, .(value = median(value, na.rm = TRUE)), by=c('feature_id', subgroupvar))
    #dt %<>% dcast.data.table(formula=as.formula(sprintf('feature_id ~ %s', subgroupvar)))
    #mat <- dt2mat(dt)
    #mat %<>% extract(rownames(object), )
    
    #tmat[is.na(tmat)] <- 0
    cormat <- propagate::bigcor(t(mat))
    cormat %<>% extract(seq_len(nrow(.)), seq_len(ncol(.)))
    apres <- apcluster::apcluster(s=cormat, details=TRUE, q=0)
    fdata(object)$cluster <- NA_character_
    exemplars <- rownames(tmat)[apres@exemplars]
    cl <- apres@clusters
    fdata(object)$cluster[unlist(cl)] <- rep(exemplars, vapply(cl, length, integer(1)))
    fdata(object)$cluster %<>% factor(exemplars)
    #heatmapres <- apcluster::heatmap(apres, cormat, col = rev(heat.colors(12)))
    
    formula <- ~ Subject + Time
    contrastdefs <- colnames(create_design(object, formula = formula))[-1]
    object %<>% fit_limma(formula = formula, contrastdefs = contrastdefs)
    values(object)[1:3, 1:3]
    
    blockeffects <- limma(object)[, levels(factor(object$Subject))[-1], 'effect']
    blockeffects %<>% data.table(keep.rownames = TRUE)
    setnames(blockeffects, 'rn', 'feature_id')
    blockeffects %<>% melt.data.table(id.vars = 'feature_id', variable.name = block, value.name = 'blockeffect')
    
    valuedt <- sumexp_to_longdt(object, svars = c('Subject', 'Time'), fvars = c('feature_id'))
    valuedt %<>% merge(blockeffects, by = c('feature_id', 'Subject'), all.x = TRUE)
    valuedt[is.na(blockeffect), blockeffect := 0]
    valuedt[, correctedvalue := value - blockeffect]
    ggplot(valuedt[feature_id==feature_id[3]], aes(x = Time, y = value,          color = Time, group = Subject)) + geom_point() + geom_line() + theme_bw()
    ggplot(valuedt[feature_id==feature_id[3]], aes(x = Time, y = correctedvalue, color = Time, group = Subject)) + geom_point() + geom_line() + theme_bw()
    
    object
}


difftype <- function(object, fit='limma', plot=TRUE){
    fitres <- metadata(object)[[fit]][rownames(object), ]
    type <- matrix('.', nrow(fitres), ncol(fitres))
    rownames(type) <- rownames(fitres)
    colnames(type) <- colnames(fitres)
    
    type[fitres[, , 'effect'] < 0 & fitres[, , 'p'     ] < 0.05] <- 'd'
    type[fitres[, , 'effect'] > 0 & fitres[, , 'p'     ] < 0.05] <- 'u'
    fdata(object)$difftype <- apply(type, 1, paste, collapse='')
    
    if (plot){
        subject <- filter_features(object, stri_detect_fixed(difftype, '.', negate = TRUE), verbose=TRUE)
        fdata(subject)$difftype %<>% factor()
        fdata(subject)$difftype
        idx <-c(which(fdata(subject)$difftype==levels(fdata(subject)$difftype)[1])[1:2],
                which(fdata(subject)$difftype==levels(fdata(subject)$difftype)[2])[1:2],
                which(fdata(subject)$difftype==levels(fdata(subject)$difftype)[3])[1:2])
        plot_subgroup_points(subject[idx, ], subgroup = Time, block = Subject, nrow=3)
    }
        
    type %<>% extract(rowAlls(!is.na(type)), )
    type %<>% apply(1, paste, collapse='')
    fdata(object)$difftype <- NA_character_
    fdata(object)$difftype
    table(fdata(object)$difftype)
    
    names(type)
    type
    
    type[data.table(fitres %>% extract(, colnames(.)[1], ))[, effect<0 & p<0.05]] <- 'd'
    type[data.table(fitres %>% extract(, colnames(.)[1], ))[, effect>0 & p<0.05]] <- 'u'
}

