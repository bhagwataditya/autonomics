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
#                    .add_pca
#                    .add_sma
#                    .add_lda
#                    .add_pls
#                    .add_spls
#                    .add_ropls
#
#============================================================================

#' @param object  SummarizedExperiment
#' @param ndim    number
#' @param verbose TRUE or FALSE (default)
#' @return        SummarizedExperiment
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' .add_pca(object)
#' @noRd
.add_pca <- function(object, ndim = 2, verbose = TRUE){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_less_than_or_equal_to(ndim, ncol(object))
    . <- NULL
# Prepare
    tmpobj <- object
    tmpobj %<>% inf_to_na(verbose=verbose)
    tmpobj %<>% nan_to_na(verbose=verbose)
    tmpobj %<>% filter_features_available_in_some_sample()
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



#' @param object  SummarizedExperiment
#' @param ndim    number
#' @return        SummarizedExperiment
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' .add_sma(object)
#' @noRd
.add_sma <- function(object, ndim = 2, verbose = TRUE){
# Assert
    if (!requireNamespace('mpm', quietly = TRUE)){
        message("First Biocinstaller::install('mpm'). Then re-run.")
        return(object)
    }
    assert_is_all_of(object, 'SummarizedExperiment')
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


#' @param object  SummarizedExperiment
#' @param ndim    number
#' @return        SummarizedExperiment
#' @author Aditya Bhagwat, Laure Cougnaud
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' .add_lda(object)
#' @noRd
.add_lda <- function(object, ndim=2, verbose = TRUE){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
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



#' @param object  SummarizedExperiment
#' @param ndim    number
#' @return        SummarizedExperiment
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' .add_pls(object)
#' @noRd
.add_pls <- function(object, ndim=2){
# Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        stop("BiocManager::install('mixOmics'). Then re-run.")
        return(object)
    }
    assert_is_all_of(object, 'SummarizedExperiment')
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
#' .add_spls(object)
#' @noRd
.add_spls <- function(object, ndim=2){
# Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        stop("BiocManager::install('mixOmics'). Then re-run.")
        return(object)
    }
    assert_is_all_of(object, 'SummarizedExperiment')
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
#' .add_opls(object)
#' @noRd
.add_opls <- function(object, ndim=2){
# Assert
    if (!requireNamespace('ropls', quietly = TRUE)){
        message("BiocManager::install('ropls'). Then re-run.")
        return(object)
    }
    assert_is_all_of(object, 'SummarizedExperiment')
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
#     plot_data()
#
#==============================================================================

#' Plot data
#' @param data    data.frame
#' @param geom    geom_point, etc.
#' @param ...     mapped aesthetics
#' @param fixed   fixed  aesthetics list
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% add_pca(plot = FALSE)
#' plot_data(sdata(object), x = pca1, y = pca2)
#' plot_data(sdata(object), x = pca1, y = pca2, color = TIME_POINT)
#' plot_data(sdata(object), x = pca1, y = pca2, color = TIME_POINT,
#'             fixed = list(shape=15, size = 3))
#' plot_data(sdata(object), x = pca1, y = pca2, shape = TIME_POINT,
#'             fixed = list(shape=15, size = 3))
#' @author Aditya Bhagwat, Johannes Graumann
#' @export
plot_data <- function(
    data, geom = geom_point, ..., fixed = list()
){
    dots <- enquos(...)
    fixed %<>% extract(setdiff(names(fixed), names(dots)))
    # https://stackoverflow.com/a/55816211
    p <- ggplot(data = data, mapping = eval(expr(aes(!!!dots))))
    p <- p + do.call(geom, fixed)
    p <- p + theme_bw()
    p
}

#' Plot sample scores
#' @param object  SummarizedExperiment
#' @param method  string: 'pca', 'pls', 'lda', 'sma'
#' @param xdim    number (default 1): x axis dimension
#' @param ydim    number (default 2): y axis dimension
#' @param color   sdata variable mapped to color
#' @param ...     additional svars mapped to aesthetics
#' @param fixed   fixed plot aesthetics
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% add_pca(plot = FALSE, ndim = 4)  # Principal Component Analysis
#' plot_sample_scores(object, 'pca')
#' plot_sample_scores(object, 'pca', x=3, y=4)
#' plot_sample_scores(object, 'pca', color = TIME_POINT)
#' @export
plot_sample_scores <- function(object, method, xdim = 1, ydim = 2,
    color = subgroup, ..., fixed = list(shape=15, size=3)
){
    x <- paste0(method, xdim)
    y <- paste0(method, ydim)
    xlab  <- paste0(x, ' : ', metadata(object)[[method]][[x]], '% ')
    ylab  <- paste0(y, ' : ', metadata(object)[[method]][[y]], '% ')

    p <- plot_data(
            sdata(object), x = !!sym(x), y = !!sym(y), color = !!ensym(color),
            ..., fixed = fixed)
    p <- p + ggplot2::xlab(xlab)
    p <- p + ggplot2::ylab(ylab)
    p
}



#============================================================================
#
#
#       add_pca, add_sma, add_pls, add_lda
#
#=============================================================================

#' Add PCA, SMA, LDA, or PLS
#'
#' Perform a dimension reduction.
#' Add sample scores, feature loadings, and dimension variances to object.
#'
#' @param object  SummarizedExperiment
#' @param ndim    number
#' @param plot    TRUE (default) or FALSE
#' @param xdim    number (default 1): x axis dimension
#' @param ydim    number (default 2): y axis dimension
#' @param color   sdata variable mapped to color
#' @param ...     additional svar to aesthetic mappings
#' @param fixed   list with fixed ggplot aesthetics
#' @param verbose TRUE (verbose) or FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' add_pca(object)  # Principal Component Analysis
#' add_pls(object)  # Partial Least Squares
#' add_lda(object)  # Linear Discriminant Analysis
#' add_sma(object)  # Spectral Map Analysis
#' add_pca(object, color = TIME_POINT)
#' add_pca(object, fixed = list(size=3, shape=1))
#' add_pca(object, xdim = 3, ydim = 4)
#' @author Aditya Bhagwat, Laure Cougnaud (LDA)
#' @export
add_pca <- function(
    object, ndim = 2, plot = TRUE, xdim = 1, ydim = 2, color = subgroup,
        ..., fixed = list(shape=15, size=3), verbose = TRUE
){
    if (verbose) message('\tAdd Principal Component Analysis')
    object %<>% .add_pca(ndim = max(ndim, xdim, ydim), verbose = verbose)
    if (plot) print(plot_sample_scores(
            object, 'pca', xdim = xdim, ydim = ydim, color = !!ensym(color),
            ..., fixed = fixed))
    object
}


#' @rdname add_pca
#' @export
add_pls <- function(
    object, ndim=2, plot=TRUE, xdim = 1, ydim = 2, color = subgroup,
    ..., fixed = list(shape = 15, size=3), verbose = TRUE
){
    if (verbose) message('\tAdd Partial Least Squares Analysis')
    object %<>% .add_pls(ndim = max(ndim, xdim, ydim))
    if (plot) print(plot_sample_scores(
            object, 'pls', xdim = xdim, ydim = ydim, color = !!ensym(color),
            ..., fixed = fixed))
    object
}


#' @rdname add_pca
#' @export
add_lda <- function(
    object, ndim = 2, plot = TRUE, xdim = 1, ydim = 2, color = subgroup,
    ..., fixed = list(shape = 15, size=3), verbose = TRUE
){
    if (verbose) message('\tAdd Linear Discriminant Analysis')
    object %<>% .add_lda(ndim = max(ndim, xdim, ydim), verbose = verbose)
    if (plot) print(plot_sample_scores(
            object, 'lda', xdim = xdim, ydim = ydim, color = !!ensym(color),
            ..., fixed = fixed))
    object
}


#' @rdname add_pca
#' @export
add_sma <- function(
    object, ndim=2, plot=TRUE, xdim=1, ydim=2, color = subgroup,
    ..., fixed = list(shape = 15, size=3), verbose = TRUE
){
    if (verbose) message('\tAdd Spectral Map Analysis')
    object %<>% .add_sma(ndim = max(ndim, xdim, ydim), verbose = verbose)
    if (plot) print(plot_sample_scores(
            object, 'sma', xdim = xdim, ydim = ydim, color = !!ensym(color),
            ..., fixed = fixed))
    object
}



