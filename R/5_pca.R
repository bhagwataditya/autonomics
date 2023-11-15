#=============================================================================
#
#                  flip_sign_if_all_exprs_are_negative.
#                  evenify_upwards
#
#=============================================================================

#' Flip sign if all expr values are negative
#' @param x matrix
#' @param verbose TRUE (default) or FALSE
#' @return updated matrix
#' @noRd
flip_sign_if_all_exprs_are_negative <- function(x, verbose=TRUE){
    idx <- !is.na(x)
    if (all(sign(x[idx])==-1)){
        if (verbose)  message('\t\tAll values negative: ', 
                                'flip signs to prevent singularities.')
        x %<>% multiply_by(-1)
    }
    x
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



#============================================================================
#
#               pca sma lda pls spls ropls
#                   .filter_minvar
#
#============================================================================

.filter_minvar <- function(object, method, minvar) {
    variances <- metadata(object)[[method]]
    discard_components <- variances[variances < minvar] %>% names()

    sdata(object)[discard_components] <- NULL
    fdata(object)[discard_components] <- NULL
    metadata(object)[[method]] <-
        variances[!names(variances) %in% discard_components]
    object
}


scorenames <- function(
    method = 'pca', by, dims = 1:2
){
    sprintf('effect~%s~%s%d', by, method, dims)
}

loadingnames <- function(
    method = 'pca', by, dims = 1:2
){
    sprintf('effect~%s~%s%d', by, method, dims)
}

methodname <- function(
    method = 'pca', by
){
    sprintf('%s~%s', by, method)
}

variancenames <- function(dims = 1:2){
    sprintf('effect%d', dims)
}

variances <- function(
    object, method = 'pca', by = biplot_by(object, method), dims = 1:2
){
    y <- metadata(object)
    y %<>% extract2(methodname(method, by))
    y %<>% extract(variancenames(dims))
    y
}

#' Extract scores/loadings
#' @param object SummarizedExperiment
#' @param method 'pca', 'pls', etc.
#' @param by      svar (string)
#' @param dims    numeric vector
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% pca()
#'     scores(object)[1:2]
#'   loadings(object)[1:2]
#'   scoremat(object)[1:2, ]
#' loadingmat(object)[1:2, ]
#' @export
scoremat <- function(
    object, method = 'pca', by = biplot_by(object, method), dims = 1:2
){
    vars <- 'sample_id'
    vars %<>% c(scorenames(method, by, dims))
    mat <- sdt(object)[, vars, with = FALSE]
    mat %<>% dt2mat()
    mat
}

#' @rdname scoremat
#' @export
scores <- function(
    object, method = 'pca', by = biplot_by(object, method), dim = 1
){
    sdt(object)[[scorenames(method, by, dim)]]
}

#' @rdname scoremat
#' @export
loadingmat <- function(
    object, method = 'pca', by = biplot_by(object, method), dims = 1:2
){
    vars <- 'feature_id'
    vars %<>% c(loadingnames(method, by, dims))
    mat <- fdt(object)[, vars, with = FALSE]
    mat %<>% dt2mat()
    mat
}

#' @rdname scoremat
#' @export
loadings <- function(
    object, method = 'pca', by = biplot_by(object, method), dim = 1
){
    fdt(object)[[loadingnames(method, by, dim)]]
}


#' PCA, SMA, LDA, PLS, SPLS, OPLS
#'
#' Perform a dimension reduction.
#' Store sample scores, feature loadings, and dimension variances.
#'
#' @param object          SummarizedExperiment
#' @param by              svar or NULL
#' @param assay           string
#' @param ndim            number
#' @param minvar          number
#' @param center_samples  TRUE/FALSE: center samples prior to pca ?
#' @param verbose         TRUE/FALSE: message ?
#' @param plot            TRUE/FALSE: plot ?
#' @param ...             passed to biplot
#' @return                SummarizedExperiment
#' @examples
#'  file <- download_data('atkin.metabolon.xlsx')
#'  object <- read_metabolon(file)
#'  pca(object, plot = TRUE)    # Principal Component Analysis
#'  pls(object, plot = TRUE)    # Partial Least Squares
#'  lda(object, plot = TRUE)    # Linear Discriminant Analysis
#'  sma(object, plot = TRUE)    # Spectral Map Analysis
#' spls(object, plot = TRUE)    # Sparse PLS
#' # opls(object, plot = TRUE)  # OPLS # outcommented because it produces a file named FALSE
#' @author Aditya Bhagwat, Laure Cougnaud (LDA)
#' @export
pca <- function(
    object, by = 'sample_id', assay = assayNames(object)[1], ndim = 2, minvar = 0, 
    center_samples = TRUE, verbose = TRUE, plot = FALSE, ...
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_scalar(assay); assert_is_subset(assay, assayNames(object))
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_less_than_or_equal_to(ndim, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    assert_is_a_bool(center_samples)
    assert_is_a_bool(verbose)
    assert_is_a_bool(plot)
    if (verbose)  message('\t\tAdd PCA')
    . <- NULL
# Prepare
    tmpobj <- object
    assays(tmpobj)[[assay]] %<>% inf_to_na(verbose=verbose)
    assays(tmpobj)[[assay]] %<>% nan_to_na(verbose=verbose)
    tmpobj %<>% rm_missing_in_all_samples(verbose = verbose)
# (Double) center and (global) normalize
    row_means <- rowMeans(        assays(tmpobj)[[assay]], na.rm = TRUE)
    col_means <- colWeightedMeans(assays(tmpobj)[[assay]], abs(row_means), na.rm = TRUE)
    global_mean <- mean(col_means)
    if (center_samples)  assays(tmpobj)[[assay]] %<>% apply(1, '-', col_means)  %>% t()  # Center columns (samples)
                         assays(tmpobj)[[assay]] %<>% apply(2, '-', row_means)           # Center rows (features)
    if (center_samples)  assays(tmpobj)[[assay]] %<>% add(global_mean)                   # Add doubly subtracted
                         assays(tmpobj)[[assay]] %<>% divide_by(sd(., na.rm=TRUE))       # Normalize
# Perform PCA
    pca_res  <- pcaMethods::pca(t(assays(tmpobj)[[assay]]),
        nPcs = ndim, scale = 'none', center = FALSE, method = 'nipals')
    samples   <- pca_res@scores
    features  <- pca_res@loadings
    variances <- round(100*pca_res@R2)
    colnames(samples)  <-    scorenames(method = 'pca', by = by, dims = seq_len(ndim))
    colnames(features) <-  loadingnames(method = 'pca', by = by, dims = seq_len(ndim))
    names(variances)   <- variancenames(seq_len(ndim))
# Add
    object %<>% merge_sdt(mat2dt(samples,   'sample_id'))
    object %<>% merge_fdt(mat2dt(features, 'feature_id'))
    metavar <- methodname(method = 'pca', by = by)
    metadata(object)[[metavar]] <- variances
# Filter for minvar
    object %<>% .filter_minvar('pca', minvar)
# Return
    if (plot)  print(biplot(object, method = 'pca', dims = seq(1,ndim)[1:2], ...))
    object
}

#' @rdname pca
#' @export
pls <- function(
    object,  by = 'subgroup', 
    assay = assayNames(object)[1], ndim = 2, 
    minvar = 0, verbose = FALSE, plot = FALSE, ...
){
# Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        message("BiocManager::install('mixOmics'). Then re-run.")
        return(object) }
    assert_is_valid_sumexp(object)
    assert_is_scalar(assay);  assert_is_subset(assay, assayNames(object))
    assert_is_subset(by, svars(object))
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Transform
    obj <- object[, !is.na(object[[by]]) ]
    x <- t(assays(obj)[[assay]])
    y <- svalues(obj, by)
    pls_out <- mixOmics::plsda( x, y, ncomp = ndim)
    samples   <- pls_out$variates$X
    features  <- pls_out$loadings$X
    variances <- round(100*pls_out$prop_expl_var$X)
    colnames(samples)  <-    scorenames(method = 'pls', by = by, dims = seq_len(ndim))
    colnames(features) <-  loadingnames(method = 'pls', by = by, dims = seq_len(ndim))
    names(variances)   <- variancenames(seq_len(ndim))
# Add
    object %<>% merge_sdt(mat2dt(samples,   'sample_id'))
    object %<>% merge_fdt(mat2dt(features, 'feature_id'))
    metavar <- methodname(method = 'pls', by = by)
    metadata(object)[[metavar]] <- variances
# Filter for minvar
    object %<>% .filter_minvar('pls', minvar)
# Return
    if (plot)  print(biplot(object, method = 'pls', dims = seq(1,ndim)[1:2], ...))
    object
}


#' @rdname pca
#' @export
sma <- function(
    object, by = 'sample_id', assay = assayNames(object)[1], ndim = 2, minvar = 0,
    verbose = TRUE, plot = FALSE, ...
){
# Assert
    if (!requireNamespace('mpm', quietly = TRUE)){
        message("First Biocinstaller::install('mpm'). Then re-run.")
        return(object)}
    assert_is_valid_sumexp(object)
    assert_is_scalar(assay);  assert_is_subset(assay, assayNames(object))
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Preprocess
    tmpobj <- object
    assays(tmpobj)[[assay]] %<>% minusinf_to_na(verbose = verbose)   # else SVD singular
    assays(tmpobj)[[assay]] %<>% flip_sign_if_all_exprs_are_negative(verbose = verbose)
    tmpobj        %<>% rm_missing_in_some_samples(verbose = verbose)
# Transform
    df <- data.frame(feature = rownames(tmpobj), assays(tmpobj)[[assay]])
    mpm_tmp <- mpm::mpm(
                df, logtrans = FALSE, closure = 'none', center = 'double',
                normal = 'global', row.weight = 'mean', col.weight = 'constant')
    ncomponents <- length(mpm_tmp$contrib)
    mpm_out <- mpm::plot.mpm(mpm_tmp, do.plot = FALSE, dim = seq_len(ncomponents))
# Extract
    samples   <- mpm_out$Columns
    features  <- mpm_out$Rows
    variances <- round(100*mpm_tmp$contrib[seq_len(ncomponents)])
    names(samples)   <-    scorenames(method = 'sma', by = by, dims = seq_len(ndim))
    names(features)  <-  loadingnames(method = 'sma', by = by, dims = seq_len(ndim))
    names(variances) <- variancenames(dims = seq_len(ndim))
# Restrict
    if (is.infinite(ndim)) ndim <- ncol(samples)
    samples   %<>% extract(, seq_len(ndim), drop = FALSE)
    features  %<>% extract(, seq_len(ndim), drop = FALSE)
    variances %<>% extract(  seq_len(ndim))
# Add
    samples  %<>% cbind( sample_id = rownames(.), .)
    features %<>% cbind(feature_id = rownames(.), .)
    object %<>% merge_sdt(data.table(samples),  'sample_id')
    object %<>% merge_fdt(data.table(features), 'feature_id')
    metavar <- methodname(method = 'sma', by = by)
    metadata(object)[[metavar]] <- variances
# Filter for minvar
    object %<>% .filter_minvar('sma', minvar)
# Return
    if (plot)  print(biplot(object, method = 'sma', dims = seq(1,ndim)[1:2], ...))
    object
}


#' @rdname pca
#' @export
lda <- function(
    object, assay = assayNames(object)[1], by = 'subgroup', ndim = 2, 
    minvar = 0, verbose = TRUE, plot = FALSE, ...
){
# Assert
    if (!requireNamespace('MASS', quietly = TRUE)){
        message("BiocManager::install('MASS'). Then re-run.")
        return(object)}
    assert_is_valid_sumexp(object)
    assert_is_scalar(assay);  assert_is_subset(assay, assayNames(object))
    assert_is_subset(by, svars(object))
    nsubgroup <- length(slevels(object, by))
    if (is.infinite(ndim))  ndim <- nsubgroup - 1
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, nsubgroup-1)
    if (ndim > (nsubgroup-1)) stop(
        sprintf('LDA requires ndim (%d) <= nsubgroup-1 (%d)',ndim, nsubgroup-1))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Preprocess
    tmpobj <- object
    assays(tmpobj)[[assay]] %<>% minusinf_to_na(verbose = verbose)         # SVD singular
    assays(tmpobj)[[assay]] %<>% flip_sign_if_all_exprs_are_negative(verbose = verbose)
    tmpobj %<>% rm_missing_in_some_samples(verbose = verbose)
# Transform
    exprs_t  <- t(assays(tmpobj)[[assay]])
    lda_out  <- suppressWarnings(MASS::lda( exprs_t,grouping = object[[by]]))
    features <- lda_out$scaling
    if (ncol(features)==1) features %<>% cbind(LD2 = 0)
    exprs_t %<>% scale(center = colMeans(lda_out$means), scale = FALSE)
    samples  <- exprs_t %*% features
    variances <- round((lda_out$svd^2)/sum(lda_out$svd^2)*100)
    features  %<>% extract(, seq_len(ndim))
    samples   %<>% extract(, seq_len(ndim))
    variances %<>% extract(  seq_len(ndim))
    if (length(variances)==1) variances <- c(LD1 = variances, LD2 = 0)
# Rename
    colnames(samples)  <-    scorenames(method = 'lda', by = by, dims = seq_len(ndim))
    colnames(features) <-  loadingnames(method = 'lda', by = by, dims = seq_len(ndim))
    names(variances)   <- variancenames(ndim)
# Restrict
    samples   %<>% extract(, seq_len(ndim), drop = FALSE)
    features  %<>% extract(, seq_len(ndim), drop = FALSE)
    variances %<>% extract(  seq_len(ndim))
# Merge - Filter - Return
    object %<>% merge_sdt(mat2dt(samples,   'sample_id'))
    object %<>% merge_fdt(mat2dt(features, 'feature_id'))
    metavar <- methodname(method = 'lda', by = by)
    metadata(object)[[metavar]] <- variances
    object %<>% .filter_minvar('lda', minvar)
    if (plot)  print(biplot(object, method = 'lda', dims = seq(1,ndim)[1:2], ...))
    object
}


#' @rdname pca
#' @export
spls <- function(
    object, assay = assayNames(object)[1], by = 'subgroup', ndim = 2, 
    minvar = 0, plot = FALSE, ...
){
# Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        message("BiocManager::install('mixOmics'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object);  assert_is_scalar(assay)
    assert_is_subset(assay, assayNames(object))
    assert_is_subset(by, svars(object))
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Transform
    x <- t(assays(object)[[assay]])
    y <- object[[by]]
    pls_out <- mixOmics::splsda( x, y, ncomp = ndim)
    samples   <- pls_out$variates$X
    features  <- pls_out$loadings$X
    variances <- round(100*pls_out$prop_expl_var$X)
    colnames(samples)  <-    scorenames(method = 'spls', by = by, dims = seq_len(ndim))
    colnames(features) <-  loadingnames(method = 'spls', by = by, dims = seq_len(ndim))
    names(variances)   <- variancenames(seq_len(ndim))
# Add
    object %<>% merge_sdt(mat2dt(samples,  'sample_id'))
    object %<>% merge_fdt(mat2dt(features,'feature_id'))
    metavar <- methodname(method = 'spls', by = by)
    metadata(object)[[metavar]] <- variances
# Filter for minvar
    object %<>% .filter_minvar('spls', minvar)
# Return
    if (plot)  print(biplot(object, method = 'spls', dims = seq(1,ndim)[1:2], ...))
    object
}


#' @rdname pca
#' @export
opls <- function(
    object, by = 'subgroup', assay = assayNames(object)[1], ndim = 2, minvar = 0, 
    verbose = FALSE, plot = FALSE, ...
){
# Assert
    if (!requireNamespace('ropls', quietly = TRUE)){
        message("BiocManager::install('ropls'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    assert_is_scalar(assay);  assert_is_subset(assay, assayNames(object))
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Transform
    x <- t(assays(object)[[assay]])
    y <- svalues(object, by)
    pls_out <- ropls::opls(x, y, predI = ndim, permI = 0, fig.pdfC = FALSE)
    samples   <- pls_out@scoreMN
    features  <- pls_out@loadingMN
    variances <- round(pls_out@modelDF$R2X*100)
    colnames(samples)  <-    scorenames(method = 'opls', by = by, dims = seq_len(ndim))
    colnames(features) <-  loadingnames(method = 'opls', by = by, dims = seq_len(ndim))
    names(variances)   <- variancenames(dims = seq_len(ndim))
# Add
    object %<>% merge_sdt(mat2dt(samples,  'sample_id'))
    object %<>% merge_fdt(mat2dt(features, 'feature_id'))
    metavar <- methodname(method = 'opls', by = by)
    metadata(object)[[metavar]] <- variances
# Filter for minvar
    object %<>% .filter_minvar('opls', minvar)
# Return
    if (plot)  print(biplot(object, method = 'opls', dims = seq(1,ndim)[1:2], ...))
    object
}


#=============================================================================
#
#           biplot
#               add_scores
#               add_loadings
#                   headtail
#
#==============================================================================

num2char <- function(x){
    if (is.null(x))     return(x)
    if (is.numeric(x))  return(as.character(x))
    return(x)
}

add_scores <- function(
    p, object, x = 'pca1', y = 'pca2', color = 'subgroup', 
    shape = if ('replicate' %in% svars(object)) 'replicate' else NULL,
    size = NULL, alpha = NULL, group = NULL, linetype = NULL,
    fixed = list(shape = 15, size = 3, na.rm = TRUE)
){
# manual colors require non-numerics
    if (!is.null(color)){
        if (is.numeric(color)){
            levs <- as.character(unique(sort(object[[color]])))
            object[[color]] %<>% num2char() %>% factor(levs)
        }
    }
    xsym <- sym(x)
    ysym <- sym(y)
    colorsym <- if (is.null(color)) quo(NULL) else sym(color)
    shapesym <- if (is.null(shape)) quo(NULL) else sym(shape)
    sizesym  <- if (is.null(size))  quo(NULL) else sym(size )
    alphasym <- if (is.null(alpha)) quo(NULL) else sym(alpha)
    groupsym <- if (is.null(group)) quo(NULL) else sym(group)
    linetypesym <- if (is.null(linetype)) quo(NULL) else sym(linetype)
# Points    
    p <- p + layer(  
                geom     = 'point',
                mapping  = aes(x = !!xsym, 
                               y = !!ysym, 
                           color = !!colorsym, 
                           shape = !!shapesym, 
                           size  = !!sizesym, 
                           alpha = !!alphasym),
                stat     = "identity", 
                data     = sdt(object), 
                params   = fixed,
                position = 'identity')
# Paths
    if (!is.null(group))  p <- p + layer(
                geom     = 'path',
                mapping  = aes(x = !!xsym, 
                               y = !!ysym, 
                           color = !!colorsym, 
                           group = !!groupsym, 
                           linetype = !!linetypesym),
                stat     = "identity",
                data     = sdt(object),
                params   = list(size = 0.1, na.rm = TRUE),
                position = 'identity')
    p
}

headtail <- function(x, n){
    c(x[seq(1, n)], x[seq(length(x)+1-n, length(x))])
}

pca1 <- pca2 <- NULL
add_loadings <- function(
    p, object, x = 'pca1', y = 'pca2', label = 'feature_name', nx = 1, ny = 1
){
# Process args
    if (nx==0 & ny==0) return(p)
    assert_is_subset(c(x, y), fvars(object))
    assert_is_subset(c(x, y), svars(object))
    axis <- angle <- NULL
# Loadings
    xloadings <- fdt(object)[[x]]
    yloadings <- fdt(object)[[y]]
    xscores   <- sdt(object)[[x]]
    yscores   <- sdt(object)[[y]]
    maxscore <- min(abs(min(c(xscores, yscores, na.rm = TRUE))),
                    abs(max(c(xscores, yscores, na.rm = TRUE))), na.rm = TRUE)
    scorefactor <- maxscore/max(abs(c(xloadings, yloadings)),  na.rm = TRUE)
    idx1 <- order(abs(xloadings), decreasing = TRUE)[seq_len(nx)] 
    idx2 <- order(abs(yloadings), decreasing = TRUE)[seq_len(ny)]
    #idx1 <- headtail(order(xloadings, na.last = NA), nx)
    #idx2 <- headtail(order(yloadings, na.last = NA), ny)
    #idx <- unique(c(idx1, idx2))
    idx <- c(idx1, idx2)
    loadingdt1 <- fdt(object)[idx1, c(label, x, y), with = FALSE]
    loadingdt2 <- fdt(object)[idx2, c(label, x, y), with = FALSE]
    loadingdt1[, axis := split_extract_fixed(x, '~', 3)]
    loadingdt2[, axis := split_extract_fixed(y, '~', 3)]
    loadingdt <- rbind(loadingdt1, loadingdt2)
    loadingdt[[x]] %<>% multiply_by(scorefactor) # bring them on same scale
    loadingdt[[y]] %<>% multiply_by(scorefactor)
    loadingdt[[x]] %<>% multiply_by(1.5)         # bring them somewhat outside
    loadingdt[[y]] %<>% multiply_by(1.5)
    loadingdt$angle <- loadingdt[[y]] / loadingdt[[x]]
    loadingdt$angle <- atan(loadingdt$angle)
    loadingdt$angle %<>% multiply_by(180/pi)

# Plot
    p <- p + geom_segment(
               data = loadingdt, 
               aes(x = 0, xend = !!sym(x), y = 0, yend = !!sym(y), linetype = axis), color = 'gray85')
    p <- p + geom_text(
                data = loadingdt, 
                aes(x = !!sym(x), y = !!sym(y), label = !!sym(label), angle = angle), 
                hjust = 'inward', color = 'gray30')
    p
}

pca1 <- pca2 <- feature_name <- NULL

#' Make alpha palette
#' @param object SummarizedExperiment
#' @param alpha string
#' @return character vector
#' @examples 
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' make_alpha_palette(object, 'Time')
#' @export
make_alpha_palette <- function(object, alpha){
# Assert
    assert_is_valid_sumexp(object)
    if (is.null(alpha))  return(NULL)
    assert_scalar_subset(alpha, svars(object))
# Create
    levels <- slevels(object, alpha)
    palette <- seq(1, 0.4, length.out = length(levels))
    names(palette) <- levels
# Return
    palette
}
    
biplot_methods <- function(object){
    y <- grep('(pca|pls)', svars(object), value = TRUE)
    y %<>% split_extract_fixed('~', 3)
    y <- gsub('[0-9]+', '', y)
    y %<>% unique()
    y
}

biplot_by <- function(object, method = 'pca'){
    y <- grep(method, svars(object), value = TRUE, fixed = TRUE)
    y %<>% split_extract_fixed('~', 2)
    y %<>% unique()
    y
}

biplot_dims <- function(
    object, method = 'pca', by = biplot_by(object, method)
){
    x <- sprintf('effect~%s~%s', by, method)
    y <- grep(x, svars(object), value = TRUE, fixed = TRUE)
    y <- gsub(x, '', y)
    y %<>% as.numeric()
    y
}

#' Biplot
#' @param object         SummarizedExperiment
#' @param method         'pca', 'pls', 'lda', 'spls', 'opls', 'sma'
#' @param by             svar
#' @param dims           numeric vector: e.g. 1:2
#' @param alpha          svar
#' @param color          svar
#' @param shape          svar
#' @param size           svar
#' @param label          svar
#' @param group          svar
#' @param linetype       svar
#' @param feature_label  fvar
#' @param fixed          fixed plot aesthetics
#' @param nx     number of x features to plot
#' @param ny     number of y features to plot
#' @param colorpalette   character vector
#' @param alphapalette   character vector
#' @param title          string
#' @param theme          ggplot2::theme output
#' @return ggplot object
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% pca(ndim = 4)
#' object %<>% pls(ndim = 4)
#' biplot(object)
#' biplot(object, nx = 1)
#' biplot(object, dims = 3:4, nx = 1)
#' biplot(object, method = 'pls')
#' biplot(object, method = 'pls', dims = 3:4)
#' biplot(object, method = 'pls', dims = 3:4, group = 'Subject')
#' @export
biplot <- function(
    object, 
    method        = biplot_methods(object)[1],
    by            = biplot_by(object, method)[1], 
    dims          = biplot_dims(object, method, by)[1:2],
    color         = 'subgroup', 
    shape         = NULL, 
    size          = NULL, 
    alpha         = NULL,
    group         = NULL, 
    linetype      = NULL,
    label         = NULL, 
    feature_label = if ('gene' %in% fvars(object)) 'gene' else 'feature_id', 
    fixed         = list(shape = 15, size = 3), 
    nx    = 0,
    ny    = 0,
    colorpalette  =  make_svar_palette(object, color),
    alphapalette  = make_alpha_palette(object, alpha), 
    title         = sprintf('%s~%s', method, by), 
    theme         = ggplot2::theme(plot.title = element_text(hjust = 0.5), 
                                   panel.grid = element_blank())
){
# Assert / Process
    assert_is_all_of(object, 'SummarizedExperiment')
    if (!is.null(color)){ assert_is_a_string(color)
                          assert_is_subset(color, svars(object)) }
    if (!is.null(group)){ assert_is_a_string(group)
                          assert_is_subset(group, svars(object)) }
    if (!is.null(shape)){ assert_is_a_string(shape)
                          assert_is_subset(shape, svars(object)) 
                          fixed %<>% extract(names(.) %>% setdiff('shape'))}
    if (!is.null(size)){  assert_is_a_string(size)
                          assert_is_subset(size,  svars(object)) 
                          fixed %<>% extract(names(.) %>% setdiff('size'))}
    
    ndim <- max(dims)
    x <- scorenames(method, by = by, dims = dims[[1]])
    y <- scorenames(method, by = by, dims = dims[[2]])
    vars <- round(variances(object, method = method, by = by, dims = dims))
    xlab <- sprintf('X%d : %d%%', dims[[1]], vars[[1]])
    ylab <- sprintf('X%d : %d%%', dims[[2]], vars[[2]])

# Loadings
# Plot
    p <- ggplot() + theme_bw() + theme
    p <- p + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) 
    p <- p + ggtitle(title)
    p %<>% add_loadings(object, x = x, y = y, label = feature_label, nx = nx, ny = ny)
    p %<>% add_scores(object, x = x, y = y, color = color, shape = shape, 
                      size = size, alpha = alpha, group = group, linetype = linetype, fixed = fixed)
    if (!is.null(colorpalette))  p <- p + scale_color_manual(values = colorpalette, na.value = 'gray80')
    if (!is.null(alphapalette))  p <- p + scale_alpha_manual(values = alphapalette)
    if (!is.null(label  ))  p <- p + geom_text_repel(
                    aes(x = !!sym(x), y = !!sym(y), label = !!sym(label)), 
                    data = sdt(object), na.rm = TRUE)
    if (!is.null(shape)){
        n <- length(slevels(object, shape))
        if (n > 6)  p <- p + scale_shape_manual(values = seq(15, 15+n-1))
            # Warning messages: The shape palette can deal with a maximum 
            # of 6 discrete values
            # https://stackoverflow.com/questions/16813278
    }
# Return
    p
}




#=============================================================================
#
#                       biplot_corrections()
#                       biplot_covariates()
#
#==============================================================================

#' @export
#' @rdname biplot_corrections
plot_corrections <- function(...){
    .Deprecated("biplot_corrections")
    biplot_corrections(...)
}

subgroup <- NULL
#' Biplot batch corrections
#'
#' @param object      SummarizedExperiment
#' @param method      'pca', 'pls', 'lda', or 'sma'
#' @param color       variable mapped to color (symbol)
#' @param covariates  covariates to be batch-corrected
#' @param varcols     number of covariate columns
#' @param plot        TRUE/FALSE: plot?
#' @param ...         used to maintain deprecated functions
#' @return  grid object
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, pca=TRUE, plot = FALSE)
#' biplot_corrections(
#'     object,  color = Group, covariates = c('SEX', 'T2D', 'SUB', 'SET'))
#' @seealso biplot_covariates
#' @export
biplot_corrections <- function(
    object, method = 'pca', color = subgroup, covariates = character(0),
    varcols = ceiling(sqrt(1+length(covariates))), plot = TRUE
){
    x <- paste0(method, "1")
    y <- paste0(method, "2")
    p <- biplot(object, !!sym(x), !!sym(y), color = !!enquo(color),
                nloadings=0)
    p <- p + ggtitle('INPUT')
    p <- p + guides(color=FALSE, fill=FALSE)
    plotlist <- list(p)
    for (ibatch in covariates){
        tmp_object <- object
        values(tmp_object) %<>%
            removeBatchEffect(batch=sdata(tmp_object)[[ibatch]])
        tmp_object <- get(method)(tmp_object, ndim=2, verbose=FALSE)
        p <- biplot(tmp_object, !!sym(x), !!sym(y), color = !!enquo(color),
                    nloadings=0)
        p <- p + ggtitle(paste0(' - ', ibatch))
        p <- p + guides(color=FALSE, fill=FALSE)
        plotlist %<>% c(list(p))
    }
    pp <- arrangeGrob(grobs = plotlist, ncol = varcols)
    if (plot) grid::grid.draw(pp)
    invisible(pp)
}


#' @rdname biplot_covariates
#' @export
plot_covariates <- function(...){
    .Deprecated('biplot_covariates')
    biplot_covariates(...)
}


#' Biplot covariates
#'
#' @param object     SummarizedExperiment
#' @param method     'pca', 'pls', 'lda', or 'sma'
#' @param covariates  covariates: mapped to color or batch-corrected
#' @param ndim        number of dimensions to plot
#' @param dimcols     number of dimension columns
#' @param varcols     number of covariate columns
#' @param plot        TRUE or FALSE: whether to plot
#' @param ...         used to maintain deprecated functions
#' @return  ggplot object
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, pca = TRUE, plot = FALSE)
#' biplot_covariates(object, covariates = 'Group', ndim = 12, dimcols = 3)
#' biplot_covariates(object, covariates = c('SEX', 'T2D', 'SUB', 'SET'))
#' biplot_covariates(object, covariates = c('SEX', 'T2D', 'SUB', 'SET'), ndim=2)
#' biplot_covariates(object, covariates = c('Group'), dimcols = 3)
#' @seealso biplot_corrections
#' @export
biplot_covariates <- function(
    object, method = 'pca', covariates = 'subgroup', ndim = 6,
    dimcols = 1, varcols = length(covariates), plot = TRUE
){
    x <- y <- NULL
    plotdt <- prep_covariates(object, method = method, ndim=ndim)
    plotlist <- list()
    for (covar in covariates){
        p <- plot_data(plotdt, geom = geom_point, x=x, y=y, color=!!sym(covar),
                        fixed = list(shape=15, size=3))
        p <- p + facet_wrap(~dims, ncol = dimcols, scales = 'free')
        p <- p + xlab(NULL) + ylab(NULL) + ggtitle(covar)
        p <- p + theme(legend.position = 'bottom', legend.title=element_blank())
        plotlist %<>% c(list(p))
    }
    pp <- gridExtra::arrangeGrob(grobs = plotlist, ncol = varcols)
    if (plot) grid::grid.draw(pp)
    invisible(pp)
}

prep_covariates <- function(object, method='pca', ndim=6){
    . <- NULL
    plotdt <- cbind(sdata(object)[FALSE,], x= character(0), y = character(0))
    projdt <- data.table(sdata(get(method)(object, ndim=ndim, verbose=FALSE)))
    alldims <- names(projdt) %>% extract(stri_detect_fixed(., method)) %>%
                stri_replace_first_fixed(method, '') %>% as.numeric()
    ndim <- min(c(max(alldims), ndim))
    npairs <- ndim %/% 2
    for (idim in seq_len(npairs)){
        dim1 <- idim*2-1
        dim2 <- idim*2
        xvar <- paste0(method, dim1)
        yvar <- paste0(method, dim2)
        tmpdt <- data.table::copy(projdt)
        setnames(tmpdt, c(xvar, yvar), c('x', 'y'))
        tmpdt %<>% extract(, stri_detect_fixed(
                                names(.), method, negate = TRUE), with = FALSE)
        tmpdt$dims <- paste0(dim1, ':', dim2)
        plotdt %<>% rbind(tmpdt)
    }
    plotdt$dims %<>% factor(unique(.))
    plotdt
}


