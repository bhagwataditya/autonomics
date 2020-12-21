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
        if (verbose) cmessage(
            '\t\tAll values negative: flip signs to prevent singularities.')
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
#' file <- download_data('halama18.metabolon.xlsx')
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
    sdata(object) %<>%  merge(df, by = by, all.x = TRUE, sort = FALSE) %>%
                        set_rownames(rownames(sdata(object))) # merging drops!
    if ('subgroup'  %in% svars(object)) object$subgroup  #%<>% as.character()
    if ('replicate' %in% svars(object)) object$replicate #%<>% as.character()
    object
}


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


#' Add PCA, SMA, LDA, PLS
#'
#' Perform a dimension reduction.
#' Add sample scores, feature loadings, and dimension variances to object.
#'
#' @param object  SummarizedExperiment
#' @param ndim    number
#' @param minvar  number
#' @param verbose TRUE (default) or FALSE
#' @param plot    TRUE/FALSE
#' @param ...     passed to biplot
#' @return        SummarizedExperiment
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' pca(object, plot=TRUE)  # Principal Component Analysis
#' pls(object)  # Partial Least Squares
#' lda(object)  # Linear Discriminant Analysis
#' sma(object)  # Spectral Map Analysis
#' pca(object, ndim=3)
#' pca(object, ndim=Inf, minvar=5)
#' @seealso mofa
#' @author Aditya Bhagwat, Laure Cougnaud (LDA)
#' @export
pca <- function(
    object, ndim = 2, minvar = 0, verbose = TRUE, plot = FALSE, ...
){
# Assert
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_less_than_or_equal_to(ndim, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
    if (verbose)  message('\tAdd PCA')
# Prepare
    tmpobj <- object
    exprs(tmpobj) %<>% inf_to_na(verbose=verbose)
    exprs(tmpobj) %<>% nan_to_na(verbose=verbose)
    tmpobj %<>% rm_missing_in_all_samples(verbose = verbose)
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
# Filter for minvar
    object %<>% .filter_minvar('pca', minvar)
# Return
    pca1 <- pca2 <- NULL
    if (plot)  print(biplot(object, sym(pca1), pca2, ...))
    object
}



#' @rdname pca
#' @export
sma <- function(
    object, ndim = 2, minvar = 0, verbose = TRUE, plot = TRUE, ...
){
# Assert
    if (!requireNamespace('mpm', quietly = TRUE)){
        message("First Biocinstaller::install('mpm'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Preprocess
    tmpobj <- object
    exprs(tmpobj) %<>% minusinf_to_na(verbose = verbose)   # else SVD singular
    exprs(tmpobj) %<>% flip_sign_if_all_exprs_are_negative(verbose = verbose)
    tmpobj        %<>% rm_missing_in_some_samples(verbose = verbose)
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
# Filter for minvar
    object %<>% .filter_minvar('sma', minvar)
# Return
    sma1 <- sma2 <- NULL
    if (plot)  print(biplot(object, sma1, sma2, ...))
    object
}


#' @rdname pca
#' @export
lda <- function(
    object, ndim = 2, minvar = 0, verbose = TRUE, plot = FALSE, ...
){
# Assert
    if (!requireNamespace('MASS', quietly = TRUE)){
        message("BiocManager::install('MASS'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    nsubgroup <- length(subgroup_levels(object))
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
    exprs(tmpobj) %<>% minusinf_to_na(verbose = verbose)         # SVD singular
    exprs(tmpobj) %<>% flip_sign_if_all_exprs_are_negative(verbose = verbose)
    tmpobj %<>% rm_missing_in_some_samples(verbose = verbose)
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
# Filter for minvar
    object %<>% .filter_minvar('lda', minvar)
# Return
    lda1 <- lda2 <- NULL
    if (plot)  print(biplot(object, lda1, lda2, ...))
    object
}


#' @rdname pca
#' @export
pls <- function(
    object, ndim = 2, minvar = 0, verbose = FALSE, plot = FALSE, ...
){
# Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        stop("BiocManager::install('mixOmics'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
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
# Filter for minvar
    object %<>% .filter_minvar('pls', minvar)
# Return
    pls1 <- pls2 <- NULL
    if (plot)  print(biplot(object, pls1, pls2, ...))
    object
}


#' @param object  SummarizedExperiment
#' @param ndim    number
#' @return        SummarizedExperiment
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' spls(object)
#' @noRd
spls <- function(object, ndim = 2, minvar = 0, plot = FALSE, ...){
# Assert
    if (!requireNamespace('mixOmics', quietly = TRUE)){
        stop("BiocManager::install('mixOmics'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Transform
    x <- t(exprs(object))
    y <- subgroup_values(object)
    pls_out <- mixOmics::splsda( x, y, ncomp = ndim)
    samples   <- pls_out$variates$X
    features  <- pls_out$loadings$X
    variances <- round(100*pls_out$explained_variance$X)
    colnames(samples)  <- sprintf('spls%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('spls%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('spls%d', seq_len(length(variances)))
# Add
    object %<>% merge_sdata(samples)
    object %<>% merge_sdata(features)
    metadata(object)$spls <- variances
# Filter for minvar
    object %<>% .filter_minvar('spls', minvar)
# Return
    spls1 <- spls2 <- NULL
    if (plot)  print(biplot(object, spls1, spls2, ...))
    object
}


#' @param object  SummarizedExperiment
#' @param ndim    number
#' @return        SummarizedExperiment
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' opls(object)
#' @noRd
opls <- function(
    object, ndim = 2, minvar = 0, verbose = FALSE, plot = FALSE, ...
){
# Assert
    if (!requireNamespace('ropls', quietly = TRUE)){
        message("BiocManager::install('ropls'). Then re-run.")
        return(object)
    }
    assert_is_valid_sumexp(object)
    if (is.infinite(ndim)) ndim <- ncol(object)
    assert_is_a_number(ndim)
    assert_all_are_in_range(ndim, 1, ncol(object))
    assert_is_a_number(minvar)
    assert_all_are_in_range(minvar, 0, 100)
    . <- NULL
# Transform
    x <- t(exprs(object))
    y <- subgroup_values(object)
    pls_out <- ropls::opls(x, y, predI = ndim, permI = 0, fig.pdfC = FALSE)
    samples   <- pls_out@scoreMN
    features  <- pls_out@loadingMN
    variances <- round(pls_out@modelDF$R2X*100)
    colnames(samples)  <- sprintf('opls%d', seq_len(ncol(samples)))
    colnames(features) <- sprintf('opls%d', seq_len(ncol(features)))
    names(variances)   <- sprintf('opls%d', seq_len(length(variances)))
# Add
    object %<>% merge_sdata(samples)
    object %<>% merge_fdata(features)
    metadata(object)$opls <- variances
# Filter for minvar
    object %<>% .filter_minvar('opls', minvar)
# Return
    opls1 <- opls2 <- NULL
    if (plot)  print(biplot(object, opls1, opls2, ...))
    object
}


#' mofa/snf
#'
#' Perform a platform-aggregating dimension reduction.
#'
#' Add sample scores, feature loadings, and dimension variances to object.
#'
#' @param object    MultiAssayExperiment
#' @param mofafile  MOFA results file
#' @param verbose   TRUE/FALSE
#' @param plot      TRUE/FALSE
#' @param ...       passed to biplot
#' @return          MultiAssayExperiment
#' @examples
#' somascan <-read_somascan( download_data('atkin18.somascan.adat'),plot=FALSE)
#' metabolon<-read_metabolon(download_data('atkin18.metabolon.xlsx'),plot=FALSE)
#' object <- sumexp2mae(list(somascan=somascan, metabolon=metabolon))
#' object %<>% mofa(plot=TRUE)
#' object %<>%  snf(plot=TRUE)
#' object$type <- substr(object$replicate,1,1)
#' biplot(object, snf1, snf2, color=type)
#' biplot(object, mofa1, mofa2, color=type)
#' @export
mofa <- function(
    object, mofafile = file.path(tempdir(), 'mofaout.hdf'),
    verbose = FALSE, plot = FALSE, ...
){
# Assert
    if (!requireNamespace('MOFA2', quietly = TRUE)){
        message("BiocManager::install('MOFA2'). Then re-run.")
        return(object)}
    assert_is_all_of(object, "MultiAssayExperiment")
    expnames <- names(experiments(object))
    for (exp in expnames)   assert_all_are_not_matching_fixed(
                                rownames(experiments(object)[[exp]]),'*')
# Run/Load mofa
    if (file.exists(mofafile)){
        if (verbose)  message('\t\tLoad ', mofafile)
        mofaobj <- MOFA2::load_model(mofafile)
    } else {
        assert_is_all_of(object, "MultiAssayExperiment")
        mofaobj <- MOFA2::prepare_mofa(MOFA2::create_mofa(object))
        mofaobj %<>% MOFA2::run_mofa(outfile = mofafile)
    }
# Add scores, loadings, variances
    scores <- mofaobj@expectations$Z$group1
    colnames(scores) %<>% gsub('Factor', 'mofa', .)
    object %<>% merge_sdata(scores)
    for (exp in expnames){
        loadings <- mofaobj@expectations$W[[exp]]
        colnames(loadings) %<>% gsub('Factor', 'mofa', .)
        experiments(object)[[exp]] %<>% merge_fdata(loadings)
    }
    variances <- rowSums(mofaobj@cache$variance_explained$r2_per_factor$group1)
    names(variances) %<>% gsub('Factor', 'mofa', .)
    metadata(object)$mofa <- variances
# Return
    mofa1 <- mofa2 <- NULL
    if (plot)  print(biplot(object, mofa1, mofa2, ...))
    object
}


#' Create MultiAssayExperiment from SummarizedExperiment list
#' @param experiments named list of SummarizedExperiments
#' @return MultiAssayExperiment
#' @examples
#' somascanfile  <- download_data('atkin18.somascan.adat')
#' metabolonfile <- download_data('atkin18.metabolon.xlsx')
#' somascan      <- read_somascan(somascanfile, plot=FALSE)
#' metabolon     <- read_metabolon(metabolonfile,plot=FALSE)
#' object        <- sumexp2mae(list(somascan=somascan, metabolon=metabolon))
#' @export
sumexp2mae <- function(experiments){
    assert_is_list(experiments)
    assert_has_names(experiments)
    for (experiment in experiments){
        assert_is_all_of(experiment, 'SummarizedExperiment')
        assert_is_subset(c('sample_id', 'subgroup'), svars(experiment))
    }
    for (i in seq_along(experiments))  experiments[[i]] %<>%
                                            extract(, order(colnames(.)))
    extract_sdata <- function(sumexp){
        extractvars <- c('sample_id', 'subgroup', 'replicate')
        extractvars %<>% intersect(svars(sumexp))
        sdata(sumexp)[, extractvars, drop=FALSE]
    }
    sdata1 <- unique(Reduce(rbind, lapply(experiments, extract_sdata)))
    sdata1 %<>% extract(order(.$sample_id), )
    assert_all_are_true(table(sdata1$sample_id)==1)
    MultiAssayExperiment(experiments = experiments, colData = sdata1)
}


snf <- function(object, ndim = 2, plot = FALSE, ...){
    # Assert
    if (!requireNamespace('SNFtool', quietly = TRUE)){
        message("BiocManager::install('SNFtool'). Then re-run.")
        return(object)}
    assert_is_all_of(object, "MultiAssayExperiment")
    assert_all_are_true(Reduce('==', colnames(object)))
    affinitymatrices <- Map(
        function(exp){  exprsmat <- exprs(experiments(object)[[exp]])
                        k <- which(is.na(exprsmat), arr.ind=TRUE)
                        exprsmat[k] <- rowMeans(exprsmat, na.rm=TRUE)[k[, 1]]
                            # https://stackoverflow.com/questions/6918086
                        exprsmat %<>% t()
                        exprsmat %<>% SNFtool::standardNormalization()
                        distmat <- SNFtool::dist2(exprsmat, exprsmat)
                        SNFtool::affinityMatrix(distmat)},
        names(experiments(object)))
    fusedaffinities <- SNFtool::SNF(affinitymatrices)
    mdsscores <- cmdscale(max(fusedaffinities) - fusedaffinities, k = ndim)
    colnames(mdsscores) <- sprintf('snf%d', seq_len(ndim))
    object %<>% merge_sdata(mdsscores)
    metadata(object)$snf <-
        structure(rep(NA, ndim), names = sprintf('snf%d', seq_len(ndim)))
    snf1 <- snf2 <- NULL
    if (plot)  print(biplot(object, x=snf1, y=snf2, ...))
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


headtail <- function(x, n){
    c(x[seq(1, n)], x[seq(length(x)+1-n, length(x))])
}


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
    assert_is_subset(xstr, names(rowData(object)))
    assert_is_subset(ystr, names(rowData(object)))
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
                data     = plotdt, # list(alpha = 0.05, size=3),
                params   = list(alpha = 0.1, size=1, na.rm = TRUE),
                position = "identity") +
        layer(  geom     = "text",
                mapping  = aes(x = !!x, y = !!y, label = !!label),
                stat     = "identity",
                data     = plotdt,
                params   = list(alpha = 0.5, na.rm = TRUE),
                position ='identity')
}


#' Biplot
#' @param object         SummarizedExperiment
#' @param x              pca1, etc.
#' @param y              pca2, etc.
#' @param color          svar mapped to color (symbol)
#' @param label          svar mapped to label (symbol)
#' @param ...            additional svars mapped to aesthetics
#' @param feature_label  fvar mapped to (loadings) label
#' @param fixed          fixed plot aesthetics
#' @param nloadings      number of loadings per half-axis to plot
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% pca(plot=FALSE, ndim=4)
#' object %<>% pls(plot=FALSE)
#' biplot(object)
#' biplot(object, pca3, pca4)
#' biplot(object, pls1, pls2)
#' biplot(object, nloadings=1)
#' biplot(object, color = TIME_POINT)
#' biplot(object, color = NULL)
#' @export
biplot <- function(object, x=pca1, y=pca2, color = subgroup, label = NULL,
    feature_label = feature_name, ...,
    fixed = list(shape=15, size=3), nloadings = 0
){
    x     <- enquo(x)
    y     <- enquo(y)
    label <- enquo(label)
    xstr <- as_name(x)
    ystr <- as_name(y)
    methodx <- gsub('[0-9]+', '', xstr)
    methody <- gsub('[0-9]+', '', ystr)
    assert_is_subset(xstr, names(colData(object)))
    assert_is_subset(ystr, names(colData(object)))

    #object %<>% get(methodx)(ndim=xdim, verbose = FALSE)
    #object %<>% get(methody)(ndim=ydim, verbose = FALSE)
    color <- enquo(color)
    feature_label <- enquo(feature_label)
    dots  <- enquos(...)
    fixed %<>% extract(setdiff(names(fixed), names(dots)))

    xvar <- round(metadata(object)[[methodx]][[xstr]], 1)
    yvar <- round(metadata(object)[[methody]][[ystr]], 1)
    xlab  <- paste0(xstr, ' : ', xvar,'% ')
    ylab  <- paste0(ystr, ' : ', yvar,'% ')

    p <- ggplot() + theme_bw() + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
    p <- p + ggtitle(paste0(xstr, ':', ystr))
    p %<>% add_loadings(
            object, !!x, !!y, label = !!feature_label, nloadings = nloadings)
    p %<>% add_scores(object, !!x, !!y, color = !!color, !!!dots, fixed = fixed)
    p %<>% add_color_scale(!!color, data = sdata(object))

    if (!quo_is_null(label)){
        p <- p + geom_text_repel(
                    aes(x=!!x, y=!!y, label=!!label), data=sdata(object))}

    p
}

#' @rdname biplot
#' @export
plot_biplot <- function(...){
    .Deprecated('biplot')
    biplot(...)
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
#' object <- read_metabolon(file, plot = FALSE)
#' biplot_corrections(object,  covariates = c('SEX', 'T2D', 'SUB', 'SET'))
#' @seealso biplot_covariates
#' @export
biplot_corrections <- function(
    object, method = 'pca', color = subgroup, covariates = character(0),
    varcols = ceiling(sqrt(1+length(covariates))), plot = TRUE
){
    p <- biplot(object, pca1, pca2, color = !!enquo(color), nloadings=0)
    p <- p + ggtitle('INPUT')
    p <- p + guides(color=FALSE, fill=FALSE)
    plotlist <- list(p)
    for (ibatch in covariates){
        exprs(object) %<>% removeBatchEffect(batch=sdata(object)[[ibatch]])
        p <- biplot(object, pca1, pca2, color = !!enquo(color), nloadings=0)
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
#' object <- read_metabolon(file, plot = FALSE)
#' biplot_covariates(object, ndim = 12, dimcols = 3)
#' biplot_covariates(object, covariates = c('SEX', 'T2D', 'SUB', 'SET'))
#' biplot_covariates(object, covariates = c('SEX', 'T2D', 'SUB', 'SET'), ndim=2)
#' biplot_covariates(object, covariates = c('subgroup'), dimcols = 3)
#' @seealso biplot_corrections
#' @export
biplot_covariates <- function(
    object, method = 'pca', covariates = 'subgroup', ndim = 6,
    dimcols = 1, varcols = length(covariates), plot = TRUE
){
    x <- y <- NULL
    plotdt <- prep_covariates(object, method = 'pca', ndim=ndim)
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


