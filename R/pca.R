#=============================================================================
#
#                  flip_sign_if_all_exprs_are_negative.
#                  evenify_upwards
#
#=============================================================================

#' Flip sign if all expr values are negative
#' @param object SummarizedExperiment
#' @param verbose TRUE (default) or FALSE
#' @return updated object
#' @noRd
flip_sign_if_all_exprs_are_negative <- function(object, verbose=TRUE){
    idx <- !is.na(exprs(object))
    if (all(sign(exprs(object)[idx])==-1)){
        if (verbose) cmessage(
            '\t\tAll values negative: flip signs to prevent singularities.')
        exprs(object) %<>% multiply_by(-1)
    }
    object
}



#' Is even/odd?
#' @param x integer
#' @return TRUE or FALSE
#' @examples
#' is_even(13)
#' is_even(12)
#' is_odd(13)
#' is_odd(12)
#' @noRd
is_even <- function(x)   (x %% 2) == 0
is_odd  <- function(x)    !is_even(x)


#' Has even/odd length?
#' @param x vector
#' @return logical
#' @examples
#' has_even_length(1:2)
#' has_odd_length(1:2)
#' has_even_length(1:3)
#' has_odd_length(1:3)
#' @noRd
has_even_length <- function(x)   is_even(length(x))
has_odd_length  <- function(x)   is_odd(length(x))


#' Evenify upwards
#' @param x integer
#' @return integer
#' @examples
#' evenify_upwards(3)
#' @noRd
evenify_upwards <- function(x)   if (is_odd(x)) x+1 else x



#=============================================================================
#
#                               merge_sdata
#                               merge_fdata
#
#=============================================================================

merge_fdata <- function(object, df, by = 'feature_id'){
    df %<>% as.data.frame() # convert matrix
    if (!'feature_id' %in% names(df))  df$feature_id <- rownames(df)
    duplicate_cols <- setdiff(intersect(fvars(object), names(df)), 'feature_id')
    fdata(object)[duplicate_cols] <- NULL
    fdata(object) %<>% merge(df, by = by, all.x = TRUE, sort = FALSE)
    fnames(object) <- fdata(object)$feature_id # merging drops them!
    object
}


#'@examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% merge_sdata( data.frame(sample_id = object$sample_id,
#'                                     number = seq_along(object$sample_id)))
#' sdata(object)
#' @noRd
merge_sdata <- function(object, df, by = 'sample_id'){
    df %<>% as.data.frame() # convert matrix to df
    if (!'sample_id' %in% names(df))  df$sample_id <- rownames(df)
    duplicate_cols <- setdiff(intersect(svars(object), names(df)), 'sample_id')
    sdata(object)[duplicate_cols] <- NULL
    sdata(object) %<>% merge(df, by = by, all.x = TRUE, sort = FALSE)
    snames(object) <- object$sample_id # merging drops them!
    if ('subgroup'  %in% svars(object)) object$subgroup  %<>% as.character()
    if ('replicate' %in% svars(object)) object$replicate %<>% as.character()
    object
}


#============================================================================
#
#                    pca, sma, lda, pls, spls, ropls
#
#============================================================================

#' Add PCA, SMA, LDA, or PLS
#'
#' Perform a dimension reduction.
#' Add sample scores, feature loadings, and dimension variances to object.
#'
#' @param object  SummarizedExperiment
#' @param ndim    number
#' @param verbose TRUE (default) or FALSE
#' @return        SummarizedExperiment
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' pca(object)  # Principal Component Analysis
#' pls(object)  # Partial Least Squares
#' lda(object)  # Linear Discriminant Analysis
#' sma(object)  # Spectral Map Analysis
#' pca(object, ndim=3)
#' @author Aditya Bhagwat, Laure Cougnaud (LDA)
#' @export
pca <- function(object, ndim = 2, verbose = TRUE){
# Assert
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_less_than_or_equal_to(ndim, ncol(object))
    . <- NULL
    if (verbose)  message('\tAdd PCA')
# Prepare
    tmpobj <- object
    tmpobj %<>% inf_to_na(verbose=verbose)
    tmpobj %<>% nan_to_na(verbose=verbose)
    tmpobj %<>% filter_features_available_in_some_sample(verbose = verbose)
# (Double) center and (global) normalize
    row_means <- rowMeans(exprs(tmpobj), na.rm=TRUE)
    col_means <- colWeightedMeans(exprs(tmpobj), abs(row_means), na.rm = TRUE)
    global_mean <- mean(col_means)
    exprs(tmpobj) %<>% apply(1, '-', col_means)   %>%   # Center columns
                        apply(1, '-', row_means)  %>%   # Center rows
                        add(global_mean)          %>%   # Add doubly subtracted
                        divide_by(sd(., na.rm=TRUE))    # Normalize
# Perform PCA
    pca_res  <- pcaMethods::pca(t(exprs(tmpobj)),
        nPcs = ndim, scale = 'none', center = FALSE, method = 'nipals')
    samples   <- pca_res@scores
    features  <- pca_res@loadings
    variances <- round(100*pca_res@R2)
    colnames(samples)  <- sprintf('pca%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('pca%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('pca%d', seq_len(length(variances)))
# Add
    object %<>% merge_sdata(samples)
    object %<>% merge_fdata(features)
    metadata(object)$pca <- variances
# Return
    object
}



#' @rdname pca
#' @export
sma <- function(object, ndim = 2, verbose = TRUE){
# Assert
    if (!requireNamespace('mpm', quietly = TRUE)){
        message("First Biocinstaller::install('mpm'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    . <- NULL
# Preprocess
    tmpobj <- object
    tmpobj %<>% minusinf_to_na(verbose = verbose)   # else SVD singular
    tmpobj %<>% flip_sign_if_all_exprs_are_negative(verbose = verbose)
    tmpobj %<>% filter_features_available_in_all_samples(verbose = verbose)
# Transform
    df <- data.frame(feature = rownames(tmpobj), exprs(tmpobj))
    mpm_tmp <- mpm::mpm(
                df, logtrans = FALSE, closure = 'none', center = 'double',
                normal = 'global', row.weight = 'mean', col.weight = 'constant')
    ncomponents <- length(mpm_tmp$contrib)
    mpm_out <- mpm::plot.mpm(mpm_tmp, do.plot=FALSE, dim = seq_len(ncomponents))
# Extract
    samples   <- mpm_out$Columns
    features  <- mpm_out$Rows
    variances <- round(100*mpm_tmp$contrib[seq_len(ncomponents)])
    names(samples)   <- sprintf('sma%d', seq_len(ncol(samples)))
    names(features)  <- sprintf('sma%d', seq_len(ncol(features)))
    names(variances) <- sprintf('sma%d', seq_len(length(variances)))
# Restrict
    if (is.infinite(ndim)) ndim <- ncol(samples)
    samples   %<>% extract(, seq_len(ndim), drop = FALSE)
    features  %<>% extract(, seq_len(ndim), drop = FALSE)
    variances %<>% extract(  seq_len(ndim))
# Add
    object %<>% merge_sdata(samples)
    object %<>% merge_fdata(features)
    metadata(object)$sma <- variances
# Return
    object
}


#' @rdname pca
#' @export
lda <- function(object, ndim=2, verbose = TRUE){
# Assert
    assert_is_valid_sumexp(object)
    nsubgroup <- length(subgroup_levels(object))
    if (is.infinite(ndim))  ndim <- nsubgroup - 1
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, nsubgroup-1)
    if (ndim > (nsubgroup-1)) stop(
        sprintf('LDA requires ndim (%d) <= nsubgroup-1 (%d)',ndim, nsubgroup-1))
    . <- NULL
# Preprocess
    tmpobj <- object
    tmpobj %<>% minusinf_to_na(verbose = verbose)         # SVD singular
    tmpobj %<>% flip_sign_if_all_exprs_are_negative(verbose = verbose)
    tmpobj %<>% filter_features_available_in_all_samples(verbose = verbose)
# Transform
    exprs_t  <- t(exprs(tmpobj))
    lda_out  <- suppressWarnings(
                    MASS::lda( exprs_t,grouping = sdata(object)$subgroup))
    features <- lda_out$scaling
    if (ncol(features)==1) features %<>% cbind(LD2 = 0)
    exprs_t %<>% scale(center = colMeans(lda_out$means), scale = FALSE)
    samples  <- exprs_t %*% features
    variances <- round((lda_out$svd^2)/sum(lda_out$svd^2)*100)
    if (length(variances)==1) variances <- c(LD1 = variances, LD2 = 0)
# Rename
    colnames(samples)  <- sprintf('lda%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('lda%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('lda%d', seq_len(length(variances)))
# Restrict
    samples   %<>% extract(, seq_len(ndim), drop = FALSE)
    features  %<>% extract(, seq_len(ndim), drop = FALSE)
    variances %<>% extract(  seq_len(ndim))
# Merge
    object %<>% merge_sdata(samples)
    object %<>% merge_fdata(features)
    metadata(object)$sma <- variances
# Return
    object
}


#' @rdname pca
#' @export
pls <- function(object, ndim=2, verbose = FALSE){
# Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        stop("BiocManager::install('mixOmics'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    . <- NULL
# Transform
    x <- t(exprs(object))
    y <- subgroup_values(object)
    pls_out <- mixOmics::plsda( x, y, ncomp = ndim)
    samples   <- pls_out$variates$X
    features  <- pls_out$loadings$X
    variances <- round(100*pls_out$explained_variance$X)
    colnames(samples)  <- sprintf('pls%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('pls%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('pls%d', seq_len(length(variances)))
# Add
    object %<>% merge_sdata(samples)
    object %<>% merge_fdata(features)
    metadata(object)$pls <- variances
# Return
    object
}


#' @param object  SummarizedExperiment
#' @param ndim    number
#' @return        SummarizedExperiment
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' spls(object)
#' @noRd
spls <- function(object, ndim=2){
# Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        stop("BiocManager::install('mixOmics'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    . <- NULL
# Transform
    x <- t(exprs(object))
    y <- subgroup_values(object)
    pls_out <- mixOmics::splsda( x, y, ncomp = ndim)
    samples   <- pls_out$variates$X
    features  <- pls_out$loadings$X
    variances <- round(100*pls_out$explained_variance$X)
    colnames(samples)  <- sprintf('pls%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('pls%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('pls%d', seq_len(length(variances)))
# Add
    object %<>% merge_sdata(samples)
    object %<>% merge_sdata(features)
    metadata(object)$spls <- variances
# Return
    object
}


#' @param object  SummarizedExperiment
#' @param ndim    number
#' @return        SummarizedExperiment
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' opls(object)
#' @noRd
opls <- function(object, ndim=2){
# Assert
    if (!requireNamespace('ropls', quietly = TRUE)){
        message("BiocManager::install('ropls'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    . <- NULL
# Transform
    x <- t(exprs(object))
    y <- subgroup_values(object)
    pls_out <- ropls::opls(x, y, predI = ndim, permI = 0, fig.pdfC = FALSE)
    samples   <- pls_out@scoreMN
    features  <- pls_out@loadingMN
    variances <- round(pls_out@modelDF$R2X*100)
    colnames(samples)  <- sprintf('pls%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('pls%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('pls%d', seq_len(length(variances)))
# Add
    object %<>% merge_sdata(samples)
    object %<>% merge_fdata(features)
    metadata(object)$opls <- variances
# Return
    object
}


#=============================================================================
#
#                       biplot()
#
#==============================================================================

#' Multi-biplot
#' @param object       SummarizedExperiment
#' @param method      'pca', 'pls', 'lda', 'sma'
#' @param ndim         number
#' @param color        variable mapped to color (default subgroup)
#' @param colorscale   vector (name = subgroup, value = colordef)
#' @param fixed        fixed ggplot aesthetics
#' @param verbose      TRUE (default) or FALSE
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' multibiplot(object)
#' @export
multibiplot <- function(
    object, method = 'pca', ndim=9,
    color = subgroup, colorscale = default_colorscale(object, !!enquo(color)),
    ...,
    fixed = list(shape=15, size=3), verbose = TRUE
){

    baredt <- data.table(sdata(object))
    object %<>% get(method)(ndim=ndim, verbose=verbose)
    scoredt <- data.table(sdata(object))
    xvar <- paste0(method, 1)
    yvar <- paste0(method, 2)
    plotdt <- baredt %>% cbind(x = scoredt[[xvar]], y = scoredt[[yvar]])
    plotdt$facet <- sprintf('X1 X2   %d %% %d %%',
                        metadata(object)[[method]][[xvar]],
                        metadata(object)[[method]][[yvar]])
    #  1   2   3   4
    # 1:2 3:4 5:6 7:8
    for (i in seq(2, ndim %/% 2)){
        xdim <- 2*i-1
        ydim <- 2*i
        xvar <- paste0(method, xdim)
        yvar <- paste0(method, ydim)
        tmpdt <- baredt %>% cbind(x = scoredt[[xvar]], y = scoredt[[yvar]])
        tmpdt$facet <- sprintf('X%d X%d   %d %% %d %%',
                               xdim, ydim,
                        metadata(object)[[method]][[xvar]],
                        metadata(object)[[method]][[yvar]])
        plotdt %<>% rbind(tmpdt)
    }
    plotdt$facet %<>% factor(unique(.))

    plot_data(plotdt, x=x, y=y, color=!!enquo(color), ..., fixed = fixed) +
    facet_wrap(~facet, scales = 'free') +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle(method) +
    theme(plot.title = element_text(hjust = 0.5))

}

#' Biplot
#' @param object         SummarizedExperiment
#' @param method         'pca', 'pls', 'lda', or 'sma'
#' @param xdim           number (default 1)
#' @param ydim           number (default 2)
#' @param color          svar mapped to color (symbol)
#' @param colorscale     vector(names = svarlevels, values = colordefs)
#' @param ...            additional svars mapped to aesthetics
#' @param feature_label  fvar mapped to (loadings) label
#' @param fixed          fixed plot aesthetics
#' @param nloadings      number of loadings per half-axis to plot
#' @return ggplot object
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' biplot(object, x=pca1, y=pca2)
#' biplot(object, x=pls1, y=pls2)
#' biplot(object, xdim=3, ydim=4)
#' biplot(object, nloadings = 0)
#' biplot(object, color = TIME_POINT)
#' biplot(object, color = TIME_POINT, xdim=3, ydim=4)
#' biplot(object, color = NULL)
#' @export
biplot <- function(
    object, x, y,
    color = subgroup, colorscale = default_colorscale(object, !!enquo(color)),
    feature_label = feature_name,
    ...,
    fixed = list(shape=15, size=3), nloadings = 1
){
    x     <- enquo(x)
    y     <- enquo(y)
    xstr <- as_name(x)
    ystr <- as_name(y)
    methodx <- substr(xstr, 1, 3)
    methody <- substr(ystr, 1, 3)
    xdim <- xstr %>% substr(4, nchar(.)) %>% as.numeric()
    ydim <- ystr %>% substr(4, nchar(.)) %>% as.numeric()

    object %<>% get(methodx)(ndim=xdim, verbose = FALSE)
    object %<>% get(methody)(ndim=ydim, verbose = FALSE)
    color <- enquo(color)
    feature_label <- enquo(feature_label)
    dots  <- enquos(...)
    fixed %<>% extract(setdiff(names(fixed), names(dots)))

    xlab  <- paste0(xstr, ' : ', metadata(object)[[methodx]][[xstr]],'% ')
    ylab  <- paste0(ystr, ' : ', metadata(object)[[methody]][[ystr]],'% ')

    p <- ggplot() + theme_bw() + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
    p %<>% add_loadings(object, !!x, !!y, label = !!feature_label, nloadings = nloadings)
    p %<>% add_scores(object, !!x, !!y, color = !!color, !!!dots, fixed = fixed)
    p %<>% add_colorscale(!!color, colorscale)

    p
}

add_scores <- function(
    p, object, x = pca1, y = pca2, color = subgroup, ...,
    fixed = list(shape=15, size=3)
){
    x     <- enquo(x)
    y     <- enquo(y)
    color <- enquo(color)

    p + layer(  geom = 'point',
                mapping = aes(x = !!x, y = !!y, color = !!color, ...),
                stat    = "identity",
                data    = sdata(object),
                params  = fixed,
                position= 'identity')

}

utils::globalVariables('feature_name')

add_loadings <- function(
    p, object, x = pca1, y = pca2, label = feature_name, nloadings = 1
){
# Process args
    if (nloadings==0) return(p)
    x     <- enquo(x)
    y     <- enquo(y)
    label <- enquo(label)
    xstr <- rlang::as_name(x)
    ystr <- rlang::as_name(y)
# Loadings
    xloadings <- fdata(object)[[xstr]]
    yloadings <- fdata(object)[[ystr]]
    idx <- unique(c(headtail(order(xloadings, na.last=NA), nloadings),
                    headtail(order(yloadings, na.last=NA), nloadings)))
# Scale loadings to scoreplot
    xscores <- sdata(object)[[xstr]]
    yscores <- sdata(object)[[ystr]]
    maxscore <- min(abs(min(c(xscores, yscores, na.rm=TRUE))),
                    abs(max(c(xscores, yscores, na.rm=TRUE))), na.rm=TRUE)
    scorefactor <- maxscore/max(abs(c(xloadings, yloadings)),  na.rm=TRUE)

    plotdt <- fdata(object)
    plotdt[[xstr]] %<>% multiply_by(scorefactor)
    plotdt[[ystr]] %<>% multiply_by(scorefactor)
    plotdt %<>% extract(idx, )

# Plot
    feature_name <- NULL
    if (!'feature_name' %in% names(plotdt)){
        setnames(plotdt, 'feature_id', 'feature_name')}
    p + layer(  geom     = 'segment',
                mapping  = aes(x=0, y=0, xend=!!x, yend=!!y),
                stat     = "identity",
                data     = plotdt,
                params   = list(alpha = 0.1, size=1, na.rm = TRUE),#params   = list(alpha = 0.05, size=3),
                position = "identity") +
        layer(  geom     = "text",
                mapping  = aes(x = !!x, y = !!y, label = !!label),
                stat     = "identity",
                data     = plotdt,
                params   = list(alpha = 0.5, na.rm = TRUE),
                position ='identity')
}



headtail <- function(x, n){
    c(x[seq(1, n)], x[seq(length(x)+1-n, length(x))])
}

add_colorscale <- function(p, color, colorscale){
    if (!rlang::quo_is_null(enquo(color))){
        p <- p + scale_color_manual(values = colorscale)
    }
    p
}

