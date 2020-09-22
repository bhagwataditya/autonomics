#=============================================================================
#
#                  flip_sign_if_all_exprs_are_negative.
#                  evenify_upwards
#
#=============================================================================

#' Flip sign if all expr values are negative
#' @param object SummarizedExperiment, eSet, or EList
#' @return updated object
#' @noRd
flip_sign_if_all_exprs_are_negative <- function(object){
   idx <- !is.na(exprs(object))
   if (all(sign(exprs(object)[idx])==-1)){
      cmessage('\t\tAll values negative: flip signs to prevent singularities.')
      exprs(object) %<>% magrittr::multiply_by(-1)
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



#==============================================================================
#
#               pca, sma, lda, plsda, splsda, ropls
#
#==============================================================================


#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' str(pca(object))
#' @noRd
pca <- function(object, ndim = 2){

# Prepare
    tmpobj <- object
    tmpobj %<>% inf_to_na(verbose=TRUE)
    tmpobj %<>% nan_to_na(verbose=TRUE)
    tmpobj %<>% filter_features_nonzero_in_some_sample()
# (Double) center and (global) normalize
    row_means <- rowMeans(exprs(tmpobj), na.rm=TRUE)
    col_means <- colWeightedMeans(exprs(tmpobj), abs(row_means), na.rm = TRUE)
    global_mean <- mean(col_means)
    exprs(tmpobj) %<>% apply(1, '-', col_means)  %>%   # Center columns
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
    names(variances) <- sprintf('pca%d', seq_len(length(variances)))
# Return
    list(samples = samples, features = features, variances = variances)
}



#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' str(sma(object))
#' @noRd
sma <- function(object, ndim = 2){
# Assert
    if (!requireNamespace('mpm', quietly = TRUE)){
        message("First Biocinstaller::install('mpm'). Then re-run.")
        return(object)
    }
# Preprocess
    tmpobj <- object
    tmpobj %<>% minusinf_to_na()
    tmpobj %<>% flip_sign_if_all_exprs_are_negative() # else SVD singular
    tmpobj %<>% filter_features_available_in_all_samples()
# Transform
    df <- data.frame(feature = rownames(tmpobj), exprs(tmpobj))
    mpm_tmp <- mpm::mpm(
                df, logtrans = FALSE, closure = 'none', center = 'double',
                normal = 'global', row.weight = 'mean', col.weight = 'constant')
    ncomponents <- evenify_upwards(sum((100*mpm_tmp$contrib) > 1))
    if(ncomponents < ndim)  stop('\'ndim\' = \'', ndim, '\', but only \'',
                                 ncomponents, 'can be provided.')
    npairs  <- ncomponents/2
    pairs <- split(1:ncomponents, rep(1:npairs, each = 2))
    sma_coords <- function(x){
        y <- suppressPackageStartupMessages(
                mpm::plot.mpm(mpm_tmp, do.plot = FALSE, dim = x))
        list(features = y$Rows[   , c('X', 'Y')],
            samples   = y$Columns[, c('X', 'Y')])
    }
    mpm_out <- lapply(pairs, sma_coords)
# Extract
    samples  <- mpm_out %>% lapply(extract2, 'samples')  %>% do.call(cbind, .)
    features <- mpm_out %>% lapply(extract2, 'features') %>% do.call(cbind, .)
    variances <- round(100*mpm_tmp$contrib[1:ncomponents])
    colnames(samples)  <- sprintf('sma%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('sma%d', seq_len(ncol(features)))
    names(variances) <- sprintf('sma%d', seq_len(length(variances)))
# Return
    list(samples  = samples[ , seq_len(ndim), drop = FALSE],
        features  = features[, seq_len(ndim), drop = FALSE],
        variances = variances[seq_len(ndim)])
}




#' @author Aditya Bhagwat, Laure Cougnaud
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' str(lda(object))
#' @noRd
lda <- function(object, ndim=2){
# Assert
    nsubgroup <- length(subgroup_levels(object))
    assert_all_are_greater_than(nsubgroup, 1)
    if (ndim > (nsubgroup-1)) stop(
        sprintf('LDA requires ndim (%d) <= nsubgroup-1 (%d)',ndim, nsubgroup-1))
# Preprocess
    tmpobj <- object
    tmpobj %<>% minusinf_to_na()
    tmpobj %<>% flip_sign_if_all_exprs_are_negative()      # else SVD singular
    tmpobj %<>% filter_features_available_in_all_samples()
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
    names(variances) <- sprintf('lda%d', seq_len(length(variances)))
# Return
    list(samples  = samples[, seq_len(ndim),  drop = FALSE],
        features  = features[, seq_len(ndim), drop = FALSE],
        variances = variances[seq_len(ndim)])
}



#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' str(pls(object))
#' @noRd
pls <- function(object, ndim=2){
    # Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        stop("BiocManager::install('mixOmics'). Then re-run.")
        return(object)
    }

    # Reduce
    x <- t(exprs(object))
    y <- subgroup_values(object)
    pls_out <- mixOmics::plsda( x, y, ncomp = ndim)
    samples   <- pls_out$variates$X
    features  <- pls_out$loadings$X
    variances <- round(100*pls_out$explained_variance$X)
    colnames(samples)  <- sprintf('pls%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('pls%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('pls%d', seq_len(length(variances)))

    # Return
    list(samples = samples, features = features, var = var)
}


#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' str(spls(object))
#' @noRd
spls <- function(object, ndim=2){
    # Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        stop("BiocManager::install('mixOmics'). Then re-run.")
        return(object)
    }

    # Reduce
    x <- t(exprs(object))
    y <- subgroup_values(object)
    pls_out <- mixOmics::splsda( x, y, ncomp = ndim)
    samples   <- pls_out$variates$X
    features  <- pls_out$loadings$X
    variances <- round(100*pls_out$explained_variance$X)
    colnames(samples)  <- sprintf('pls%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('pls%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('pls%d', seq_len(length(variances)))

    # Return
    list(samples = samples, features = features, var = var)
}


#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' str(opls(object))
#' @noRd
opls <- function(object, ndim=2){
    # Assert
    if (!requireNamespace('ropls', quietly = TRUE)){
        message("BiocManager::install('ropls'). Then re-run.")
        return(object)
    }

    # Reduce
    x <- t(exprs(object))
    y <- subgroup_values(object)
    pls_out <- ropls::opls(x, y, predI = ndim, permI = 0, fig.pdfC = FALSE)
    samples   <- pls_out@scoreMN
    features  <- pls_out@loadingMN
    variances <- round(pls_out@modelDF$R2X*100)
    colnames(samples)  <- sprintf('pls%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('pls%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('pls%d', seq_len(length(variances)))

    # Return
    list(samples = samples, features = features, variances = variances)
}


#=============================================================================
#
#     plot_data()
#
#==============================================================================

#' Plot data
#' @param sdata   sample data.frame
#' @param geom    geom_point, etc.
#' @param x       x variable
#' @param y       y variable
#' @param color   color variable
#' @param ...     additional variable to aesthetic mappings
#' @param fixed   list with fixed aesthetic specifications
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' pcaresults <- pca(object)
#' object %<>% merge_sdata(pcaresults$samples)
#' plot_data(sdata(object), x = pca1, y = pca2)
#' plot_data(sdata(object), x = pca1, y = pca2, color = TIME_POINT)
#' plot_data(sdata(object), x = pca1, y = pca2, color = TIME_POINT,
#'             fixed = list(shape=15, size = 3))
#' @author Aditya Bhagwat, Johannes Graumann
#' @export
plot_data <- function(
    data, geom = geom_point, color = subgroup, ..., fixed = list()
){
    color <- enquo(color)
    p <- ggplot(
            data = data, mapping = aes(color = !!color, ...))
    p <- p + do.call(geom, fixed)
    p <- p + theme_bw()
    print(p)
    p
}



#============================================================================
#
#         add_projection
#
#============================================================================


#' Add PCA, SMA, LDA, or PLS
#'
#' Perform a dimension reduction.
#' Add sample scores, feature loadings, and dimension variances to object.
#'
#' @param object  SummarizedExperiment
#' @param ndim    number
#' @param plot    TRUE (default) or FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' add_projection(object, 'pca')
#' add_projection(object, 'pca', ndim=3, x=2, y=3)
#' add_projection(object, 'pls')
#' add_projection(object, 'lda')
#' add_projection(object, 'sma')
add_projection <- function(
  object, method, ndim = 2, plot = TRUE, x = 1, y = 2,
  fixed = list(shape=15, size=2), ...
){
    # x <- ensym(x)
    # y <- ensym(y)

    results <- get(method)(object, ndim)
    object %<>% merge_sdata(results$samples)
    object %<>% merge_fdata(results$features)
    metadata(object)$lda <- results$variances
    if (plot) plot_data(sdata(object),
                        x     = !!sym(paste0(method, x)),
                        y     = !!sym(paste0(method, y)),
                        ...,
                        fixed = fixed)
    object
}


merge_fdata <- function(object, df, by = 'feature_id'){
    df %<>% as.data.frame() # convert matrix
    if (!'feature_id' %in% names(df))  df$feature_id <- rownames(df)
    fdata(object) %<>% merge(df, by = by, all.x = TRUE, sort = FALSE)
    object
}


#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' pcaresults <- pca(object)
#' object %<>% merge_sdata(pcaresults$samples)
merge_sdata <- function(object, df, by = 'sample_id'){
    df %<>% as.data.frame() # convert matrix to df
    if (!'sample_id' %in% names(df))  df$sample_id <- rownames(df)
    sdata(object) %<>% merge(df, by = by, all.x = TRUE, sort = FALSE)
    object
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
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('glutaminase.metabolon.xlsx')
#' object <- read_metabolon(file)
#' add_pca(object)  # Principal Component Analysis
#' add_pls(object)  # Partial Least Squares
#' add_lda(object)  # Linear Discriminant Analysis
#' add_sma(object)  # Spectral Map Analysis
#' add_pca(object, color = TIME_POINT)  # Principal Component Analysis
#' add_pca(object, fixed = list(shape = 1))  # Principal Component Analysis
#' @export
add_pca <- function(
    object, ndim = 2, plot = TRUE, x = 1, y = 2, ...
){
    object %<>% add_projection(
        'pca', ndim = ndim, plot = plot, x = x, y = y, ...)
}


#' @rdname add_pca
#' @export
add_pls <- function(
  object, ndim = 2, plot = TRUE, x = 1, y = 2, ...
){
    object %<>% add_projection(
        'pls', ndim = ndim, plot = plot, x = x, y = y, ...)
}


#' @rdname add_pca
#' @export
add_lda <- function(
    object, ndim = 2, plot = TRUE, x = 1, y = 2, ...
){
    object %<>% add_projection(
        'lda', ndim = ndim, plot = plot, x = x, y = y, ...)
}


#' @rdname add_pca
#' @export
add_sma <- function(
    object, ndim = 2, plot = TRUE, x = 1, y = 2, ...
){
    object %<>% add_projection(
        'sma', ndim = ndim, plot = plot, x = x, y = y, ...)
}



