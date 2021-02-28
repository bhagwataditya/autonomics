#==============================================================================
#
#                   Switch between nondetect representations
#
#==============================================================================


#' Switch between nondetect representations
#' @param x    matrix
#' @param verbose   logical(1)
#' @return Updated matrix
#' @examples
#' x <- matrix(c(0, 2, 3, 4, 5, 6), nrow=2, byrow=TRUE)
#' x
#' zero_to_na(x)
#' @export
zero_to_na <- function(x, verbose = FALSE){
    selector <- x == 0
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0('\t\tReplace 0 -> NA for %d/%d values ',
                            '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector), na.rm=TRUE), nrow(x),
                    sum(colAnys(selector), na.rm=TRUE), ncol(x))
        x[selector] <- NA_real_
    }
    x
}


#' Convert NaN to NA
#' @param x matrix
#' @param verbose TRUE/FALSE
#' @return matrix
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' invert_subgroups <- c('EM_E', 'BM_E', 'BM_EM')
#' x <- exprs(read_proteingroups(
#'               file, invert_subgroups=invert_subgroups, plot=FALSE))
#' nan_to_na(x, verbose=TRUE)
#' @export
nan_to_na <- function(x, verbose = FALSE){
    selector <- is.nan(x)
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0( '\t\tReplace NaN -> NA for %d/%d values ',
                            '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector)), nrow(x),
                    sum(colAnys(selector)), ncol(x))
        x[selector] <- NA_real_
    }
    x
}


#' Convert NA to zero
#' @param x matrix
#' @param verbose TRUE/FALSE
#' @return matrix
#' @examples
#' x <- matrix(c(NA, 2, 3, 4, 5, 6), nrow=2, byrow=TRUE)
#' x
#' na_to_zero(x)
#' @export
na_to_zero <- function(x, verbose = FALSE){
    selector <- is.na(x)
    if (any(selector)){
        if (verbose) cmessage(
                        paste0( '\t\tReplace NA -> 0 for %d/%d values ',
                                '(in %d/%d features and %d/%d samples)'),
                        sum(selector), nrow(selector)*ncol(selector),
                        sum(rowAnys(selector)), nrow(x),
                        sum(colAnys(selector)), ncol(x))
        x[selector] <- 0
    }
    x
}


#' Convert Inf to NA
#' @param x matrix
#' @param verbose TRUE/FALSE
#' @return matrix
#' @examples
#' x <- matrix(c(-Inf, 2, 3, 4, 5, 6), nrow=2, byrow=TRUE)
#' x
#' inf_to_na(x)
#' @export
inf_to_na <- function(x, verbose = FALSE){
    selector <- is.infinite(x)
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0(
                        '\t\tReplace -Inf -> NA for %d/%d values ',
                        '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector)), nrow(x),
                    sum(colAnys(selector)), ncol(x))
        x[selector] <- NA_real_
    }
    x
}


#' Convert -Inf to NA
#' @param x matrix
#' @param verbose TRUE/FALSE
#' @return matrix
#' @examples
#' x <- matrix(c(-Inf, 2, 3, 4, 5, 6), nrow=2, byrow=TRUE)
#' x
#' minusinf_to_na(x)
#' @export
minusinf_to_na <- function(x, verbose = FALSE){
    selector <- x==-Inf
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0( '\t\tReplace -Inf -> NA for %d/%d values ',
                            '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector)), nrow(x),
                    sum(colAnys(selector)), ncol(x))
        x[selector] <- NA_real_
    }
    x
}


na_to_string <- function(x){
    x[is.na(x)] <- ''
    x
}


#=============================================================================
#
#                   normimpute
#                   halfnormimpute
#                   zeroimpute
#
#=============================================================================

#' @rdname halfnormimpute
#' @export
normimpute <- function(x, selector = is.na(x)){
    x[selector] <- rnorm(length(x[selector]), sd = sd(x[!is.na(x)]))
    x
}


#' Impute from half-normal distribution around 0
#' @param x          NA-containing numeric vector
#' @param selector   which values to impute
#' @return numeric vector of same length
#' @examples
#' require(data.table)
#' x <- rnorm(1e5)
#' idx <- runif(length(x))>0.9
#' x[idx] <- NA
#' dt1 <- data.table(value = normimpute(x), distr = 'norm')
#'
#' x <- abs(rnorm(1e5)); x[idx] <- NA
#' dt2 <- data.table(value = halfnormimpute(x), distr = 'halfnorm')
#'
#' x <- abs(rnorm(1e5)); x[idx] <- NA
#' dt3 <- data.table(value = zeroimpute(x), distr = 'zero')
#'
#' require(ggplot2)
#' ggplot(rbind(dt1,dt2,dt3), aes(x=value, fill=distr)) +
#' geom_density(alpha=0.5)
#' @export
halfnormimpute <- function(x, selector = is.na(x)){
    x[selector] <- abs(rnorm(length(x[selector]), sd = 2*sd(x[!is.na(x)])))
    x
}


#' @rdname halfnormimpute
#' @export
zeroimpute <- function(x, selector = is.na(x)){
    x[selector] <- 0
    x
}

#=============================================================================
#
#                        split_by_svar
#
#==============================================================================

#' Split by svar
#' @param object SummarizedExperiment
#' @param svar   svar to split on
#' @return list of SummarizedExperiment
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute=FALSE, plot = FALSE)
#' split_by_svar(object)
#' @export
split_by_svar <- function(object, svar = subgroup){
    svar <- enquo(svar)
    svarstr <- as_name(svar)

    if (is.null(svar)) return(list(object))
    extract_samples  <- function(sg){
                            idx <- sdata(object)[[svarstr]] == sg
                            object[, idx]
                        }
    Map(extract_samples, slevels(object, svarstr))
}



#=============================================================================
#
#                     is_systematic_detect
#                     is_random_detect
#                     is_full_detect
#
#=============================================================================


is_systematic_detect <- function(object, subgroup = subgroup){
    subgroup <- enquo(subgroup)
    split_by_svar(object, !!subgroup) %>%
    lapply(function(x) rowAlls(is.na(exprs(x)))) %>%
    Reduce("|", .)
}

is_random_detect <- function(object, subgroup = subgroup){
    subgroup <- enquo(subgroup)
    rowAnys(is.na(exprs(object))) & !is_systematic_detect(object, !!subgroup)
}


is_full_detect <- function(object){
    rowAlls(!is.na(exprs(object)))
}


#=============================================================================
#
#                     venn_detects
#
#=============================================================================


#' Venn detects
#'
#' Venn diagram full/systematic/random detects
#'
#' @param object SummarizedExperiment
#' @return  \code{NULL}
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute=FALSE, plot = FALSE)
#' venn_detects(object)
#' @export
venn_detects <- function(
    object, 
    subgroup = if ('subgroup' %in% svars(object)) subgroup else NULL
){
    subgroup <- enquo(subgroup)
    limma::vennDiagram(as.matrix(cbind(
        systematic = is_systematic_detect(object, !!subgroup),
        random     = is_random_detect(    object, !!subgroup),
        full       = is_full_detect(      object))))
}


#=============================================================================
#
#                     impute_systematic_nondetects
#
#==============================================================================


#' Impute systematic nondetects
#' @param object    SummarizedExperiment
#' @param subgroup  subgroup svar
#' @param fun       imputation function
#' @param plot      TRUE or FALSE
#' @param verbose   TRUE or FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute = FALSE, plot = FALSE)
#' impute_systematic_nondetects(object)
#' @export
impute_systematic_nondetects <- function(object, subgroup = subgroup, 
    fun = halfnormimpute, plot = TRUE, verbose = TRUE
){
# Process
    absent <- replicated <- systematic <- NULL
    subgroup <- enquo(subgroup)
    subgroupvar <- as_name(subgroup)
# Filter
    object %<>% filter_exprs_replicated_in_some_subgroup(
                    subgroupvar = subgroupvar, verbose = verbose)
# Impute
    dt  <-  sumexp_to_long_dt(object, svars = subgroupvar)
    dt[, absent     := all(is.na(value)),    by = c('feature_id', subgroupvar)]
    dt[, replicated := sum(!is.na(value))>1, by = c('feature_id', subgroupvar)]
    dt[, systematic := any(absent) & any(replicated), by = 'feature_id']
    dt[, is_imputed := systematic & absent]
    dt[, value := fun(value, is_imputed), by='feature_id']
# Update object
    ff <- fnames(object)
    ss <- snames(object)
    exprs(object) <- dt2exprs(dt)[ff, ss]
    assays(object)$is_imputed <- dt2mat(data.table::dcast(
                    dt, feature_id ~ sample_id, value.var = 'is_imputed'))
    fdata(object)$imputed <- rowAnys(assays(object)$is_imputed)
# Plot
    nrowimputed <- sum(rowAnys(is_imputed(object)))
    ncolimputed <- sum(colAnys(is_imputed(object)))
    if (verbose & nrowimputed>0)  cmessage(
        "\t\tImpute systematic nondetects for %d/%d features in %d/%d samples",
        nrowimputed, nrow(object),
        ncolimputed, ncol(object))
    if (plot)    print(plot_detections(object, subgroup = !!subgroup))
# Return
    object
}


#==============================================================================
#
#                      cluster_order_features
#                      detect_order_features
#
#==============================================================================

cluster_order_features <- function(object){
    if (nrow(object) < 3) return(object)
    idx <- is.na(exprs(object))
    exprs(object)[idx] <- 0
    order <- hclust(dist(exprs(object)))$order
    exprs(object)[idx] <- NA
    object %<>% extract(order, )
    object
}

detect_order_features <- function(object, subgroup){
    x <- object
    exprs(x)[is_imputed(x)] <- NA
    idx1 <- fnames(cluster_order_features(
                        x[is_systematic_detect(x, !!subgroup),]))
    idx2 <- fnames(cluster_order_features(
                        x[is_random_detect(x, !!subgroup), ]))
    idx3 <- fnames(cluster_order_features(
                        x[is_full_detect(x),]))
    SummarizedExperiment::rbind(object[idx1,], object[idx2,], object[idx3,])
}


#==============================================================================
#
#                           plot_detections
#
#==============================================================================


#' @rdname plot_detections
#' @export
plot_detects <- function(...){
    .Deprecated('plot_detections')
    plot_detections(...)
}


#' Plot detections
#'
#' Plot detections
#'
#' \code{plot_detections} plots feature x sample detections. It shows per
#' feature/sample nondetects (white), imputes (light colored), and detects
#' (full color).
#'
#' \code{plot_summarized_detections} gives an summarized view, plotting
#' featuretype x subgroup detections. It visualizes the subgroup-wise nondetect
#' structure often seen in mass spectrometry proteomics data (across e.g.
#' different cell types)
#' @param object     SummarizedExperiment
#' @param subgroup   subgroup svar sym
#' @param fill       fill svar
#' @param na_imputes whether to NA imputes prior to plottin (TRUE/FALSE)g
#' @param ...      for backward compatibilty
#' @return ggplot object
#' @examples
#' require(magrittr)
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute=FALSE, plot = FALSE)
#' plot_summarized_detections(object)
#' plot_detections(object)
#' plot_detections(impute_systematic_nondetects(object, plot=FALSE))
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, impute = FALSE, plot = FALSE)
#' plot_summarized_detections(object)
#' plot_detections(object)
#' @export
plot_detections <- function(object, subgroup = subgroup, fill = subgroup){
# Process
    detection <- feature_id <- NULL
    subgroup <- enquo(subgroup);         fill     <- enquo(fill)
    subgroupvar <- as_name(subgroup);    fillstr  <- as_name(fill)
# Reorder samples
    sdata(object)[[subgroupvar]] %<>% factor()
    object %<>% extract(, order(sdata(.)[[subgroupvar]]))
# Reorder/block features
    object %<>% detect_order_features(!!subgroup)
    y <- object; exprs(y)[is_imputed(y)] <- NA
    nfull       <- sum(is_full_detect(y))
    nsystematic <- sum(is_systematic_detect(y, subgroup=!!subgroup))
    nrandom     <- sum(is_random_detect(y, subgroup=!!subgroup))
# Melt
    plotdt  <-  sumexp_to_long_dt(object)
    alpha <- NULL
    plotdt[,             detection := 'detect']
    plotdt[is.na(value), detection := 'nondetect']
    if ('is_imputed' %in% SummarizedExperiment::assayNames(object)){
        plotdt[is_imputed==TRUE, detection := 'impute']}
    plotdt[, detection := factor(detection, c('nondetect', 'impute', 'detect'))]
    plotdt[, sample_id  := factor( sample_id, unique(snames(object)))]
    plotdt[, feature_id := factor(feature_id, rev(unique(fnames(object))))]
    colors <- make_colors(unique(sdata(object)[[fillstr]]))
    colors %<>% c(nondetect = "#FFFFFF")
# Plot
    ggplot(plotdt) +
    geom_tile(aes(x=sample_id, y=feature_id, fill=!!fill, alpha=detection)) +
    scale_fill_manual(values = colors) +
    scale_alpha_manual(values = c(nondetect=0, impute=0.3, detect = 1)) +
    ylab('Features') +
    xlab('Samples') +
    ggtitle(sprintf('detects: %d full, %d random, %d systematic',
                    nfull, nrandom, nsystematic)) +
    theme_bw() +
    theme(  axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
            axis.text.x  = element_text(angle = 90, vjust = 0.5),
            legend.title = element_blank(),
            panel.grid   = element_blank()) +
    geom_hline(yintercept = cumsum(c(nfull, nrandom, nsystematic))) +
    guides(alpha=FALSE)

}

#==============================================================================
#
#                           plot_summarized_detections
#
#==============================================================================


get_subgroup_combinations <- function(object, subgroupvar){
    type <- NULL
    subgroups <- slevels(object, subgroupvar)
    subgroups  %>%
        lapply(function(x) c(0,1) %>% set_names(rep(x,2))) %>%
        set_names(subgroups) %>%
        expand.grid() %>%
        data.table() %>%
        extract(rev(order(rowSums(.)))) %>%
        extract(, type := 0:(.N-1)) %>%
        extract()
}


#' @rdname plot_detections
#' @export
plot_quantifications <- function(...){
    .Deprecated('plot_summarized_detections')
    plot_summarized_detections(...)
}


#' @rdname plot_detections
#' @export
plot_summarized_detections <- function(object, subgroup = subgroup,
                                        fill = subgroup, na_imputes = TRUE){
# Assert
    assert_is_all_of(object, "SummarizedExperiment")
    subgroup <- enquo(subgroup)
    if (quo_is_null(subgroup))  return(ggplot() + geom_blank())
    subgroupvar <- as_name(subgroup)
    fill <- enquo(fill);     fillstr <- as_name(fill)
    assert_is_subset(subgroupvar, svars(object))
    assert_is_subset(fillstr,  svars(object))
    xmin <- xmax <- ymin <- ymax <- nfeature <- quantified <- NULL
# Prepare
    object %<>% filter_samples(!is.na(!!subgroup), verbose=TRUE)
    exprs(object) %<>% zero_to_na()  #### TODO fine-tune
    featuretypes <- get_subgroup_combinations(object, subgroupvar)
    dt <- sumexp_to_long_dt(object, svars = subgroupvar)
    if (na_imputes) if ('is_imputed' %in% names(dt))  dt[is_imputed==TRUE,
                                                        value := NA]
    dt %<>% extract(, .(quantified   = as.numeric(any(!is.na(value)))),
                    by = c(subgroupvar, 'feature_id'))
    dt %<>% data.table::dcast.data.table(
        as.formula(paste0('feature_id ~ ', subgroupvar)),value.var='quantified')
    dt %<>% merge(featuretypes, by = setdiff(names(featuretypes), 'type'))
    dt %<>% extract(,.(nfeature=.N),by='type')
    dt %<>% merge(featuretypes,by='type')
    dt[, ymax := cumsum(nfeature)]
    dt[, ymin := c(0,ymax[-.N])]
    dt %<>% data.table::melt.data.table(
        id.vars = c('type', 'nfeature', 'ymin', 'ymax'),
        variable.name = subgroupvar, value.name='quantified')
    dt$quantified %<>% as.factor()
    nsampledt <- data.table(sdata(object))[, .N, by=subgroupvar] %>% # preserves
                set_names(c(subgroupvar, 'xmax'))                # factor order!
    setorderv(nsampledt, subgroupvar)
    nsampledt[, xmax := cumsum(xmax)]
    nsampledt[, xmin := c(0, xmax[-.N])]
    dt %<>% merge(nsampledt, by = subgroupvar)
# Plot
    colors <- make_colors(slevels(object, fillstr))
    ggplot(dt) + geom_rect(aes( xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                fill=!!fill, alpha=quantified)) +
                geom_segment(aes(x=xmin, xend=xmax, y = ymax, yend=ymax)) +
                geom_segment(aes(x=xmin, xend=xmax, y = ymin, yend=ymin)) +
                geom_segment(aes(x=xmax, xend=xmax, y = ymin, yend=ymax)) +
                geom_segment(aes(x=xmin, xend=xmin, y = ymin, yend=ymax)) +
                theme_minimal() + xlab('Samples') + ylab('Features') +
                theme(panel.grid = element_blank()) + guides(alpha=FALSE) +
                scale_fill_manual(values = colors) +
                scale_alpha_manual(values=c(`0`=0, `1`=1)) +
                ggtitle('detections')
}


#==============================================================================
#
#                           explore_imputations
#
#==============================================================================


#' Explore imputations
#' @param object SummarizedExperiment
#' @param xbiplot biplot x axis. Default pca1 (symbol)
#' @param ybiplot biplot y axis. Default pca2 (symbol)
#' @param ... aesthetic mappings
#' @return ggplot object
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute = FALSE, pca = TRUE, plot = FALSE)
#' explore_imputations(object)
#' explore_transformations(object)
#' @export
explore_imputations <- function(object, xbiplot = pca1, ybiplot = pca2, ...){
    imputed     <- impute_systematic_nondetects(object, plot=FALSE)
    zeroed <- impute_systematic_nondetects(
                object, fun = zeroimpute, plot = FALSE)
    legend <- gglegend(biplot(object))

    do_plot_sample_detections <- function(obj, ...){
        plot_detections(obj, ...) + guides(color=FALSE, fill=FALSE)
    }

    do_biplot <- function(obj, ...){
        biplot(obj, x=!!enquo(xbiplot), y=!!enquo(ybiplot), nloadings=0,...) +
        guides(color=FALSE, fill=FALSE) +
        ggtitle(NULL)
    }

    do_plot_sample_densities <- function(obj, ...){
        plot_sample_densities(obj, ...) +
        guides(color=FALSE, fill=FALSE) +
        ggtitle(NULL)
    }

    p1 <- do_plot_sample_detections(object, ...)  + ggtitle('Original')
    p2 <- do_plot_sample_detections(imputed, ...) + ggtitle('Halfnorm imputed')
    p3 <- do_plot_sample_detections(zeroed, ...)  + ggtitle('Zero imputed')
    p4 <- do_plot_sample_densities(object,  ...)
    p5 <- do_plot_sample_densities(imputed, ...)
    p6 <- do_plot_sample_densities(zeroed,  ...)
    p7 <- do_biplot(object, ...)
    p8 <- do_biplot(imputed, ...)
    p9 <- do_biplot(zeroed, ...)
    grid.arrange(arrangeGrob(p1,p2,p3,p4,p5,p6,p7,p8,p9, nrow=3),
                            legend, ncol=2, widths = c(8, 1))
}


#==============================================================================

#'@title Get/set is_imputed
#'@description Get/Set is_imputed
#'@param object SummarizedExperiment
#'@param value matrix
#'@return matrix (get) or updated object (set)
#'@examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' sum(is_imputed(object))
#' @rdname is_imputed
#' @export
setGeneric("is_imputed",  function(object) standardGeneric("is_imputed") )

#' @rdname is_imputed
setMethod("is_imputed", signature("SummarizedExperiment"),  function(object){
    if ('is_imputed' %in% names(assays(object))){
        assays(object)$is_imputed
    } else {
        matrix(FALSE, nrow = nrow(object), ncol = ncol(object),
                dimnames = dimnames(object))
    }
})

#' @rdname is_imputed
#' @export
setGeneric(
    "is_imputed<-",
    function(object, value)  standardGeneric("is_imputed<-") )

#' @rdname is_imputed
setReplaceMethod(
    "is_imputed",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$is_imputed <- value; object})

#' @rdname is_imputed
setReplaceMethod(
    "is_imputed",
    signature("SummarizedExperiment", "NULL"),
    function(object, value){object})

