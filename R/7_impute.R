#==============================================================================
#
#                   Change nondetect representation
#
#==============================================================================


#' Change nondetect representation
#' @param x    matrix
#' @param verbose   logical(1)
#' @return Updated matrix
#' @examples
#' require(magrittr)
#' matrix(c(0, 7), nrow=1)
#' matrix(c(0, 7), nrow=1)    %>% zero_to_na(verbose=TRUE)
#' 
#' matrix(c(NA, 7), nrow=1)
#' matrix(c(NA, 7), nrow=1)   %>% na_to_zero(verbose=TRUE)
#' 
#' matrix(c(NaN, 7), nrow=1)
#' matrix(c(NaN, 7), nrow=1)  %>% nan_to_na(verbose=TRUE)
#' 
#' matrix(c(Inf, 7), nrow=1)
#' matrix(c(Inf, 7), nrow=1)  %>% inf_to_na(verbose=TRUE)
#' 
#' matrix(c(-Inf, 7), nrow=1)
#' matrix(c(-Inf, 7), nrow=1) %>% minusinf_to_na(verbose=TRUE)
#' @export
zero_to_na <- function(x, verbose = FALSE){
    selector <- x == 0
    if (any(c(selector), na.rm = TRUE)){
        if (verbose)  message('\t\tReplace 0->NA for ', 
            sum(selector, na.rm=TRUE), '/', nrow(selector)*ncol(selector), 
            ' values (in ',  sum(rowAnys(selector), na.rm=TRUE), '/', nrow(x), 
            ' features of ', sum(colAnys(selector), na.rm=TRUE), '/', ncol(x), 
            ' samples)')
        x[selector] <- NA_real_
    }
    x
}


#' @rdname zero_to_na
#' @export
nan_to_na <- function(x, verbose = FALSE){
    selector <- is.nan(x)
    if (any(c(selector), na.rm = TRUE)){
        if (verbose)  message('\t\tReplace NaN->NA for ', 
            sum(selector, na.rm=TRUE), '/', nrow(selector)*ncol(selector), 
            ' values (in ',  sum(rowAnys(selector)), '/', nrow(x), 
            ' features of ', sum(colAnys(selector)), '/', ncol(x), ' samples)')
        x[selector] <- NA_real_
    }
    x
}


#' @rdname zero_to_na
#' @export
na_to_zero <- function(x, verbose = FALSE){
    selector <- is.na(x)
    if (any(selector)){
        if (verbose)  message('\t\tReplace NA->0 for ', 
            sum(selector), '/', nrow(selector)*ncol(selector), 
            ' values (in ',  sum(rowAnys(selector)), '/', nrow(x), 
            ' features of ', sum(colAnys(selector)), '/', ncol(x), ' samples)')
        x[selector] <- 0
    }
    x
}


#' @rdname zero_to_na
#' @export
inf_to_na <- function(x, verbose = FALSE){
    selector <- is.infinite(x)
    if (any(c(selector), na.rm = TRUE)){
        if (verbose)  message('\t\tReplace -Inf->NA for ', 
            sum(selector, na.rm=TRUE), '/', nrow(selector)*ncol(selector), 
            ' values (in ',  sum(rowAnys(selector)), '/', nrow(x), 
            ' features of ', sum(colAnys(selector)), '/', ncol(x), ' samples)')
        x[selector] <- NA_real_
    }
    x
}


#' @rdname zero_to_na
#' @export
minusinf_to_na <- function(x, verbose = FALSE){
    selector <- x==-Inf
    if (any(c(selector), na.rm = TRUE)){
        if (verbose)  message('\t\tReplace -Inf->NA for ', 
            sum(selector, na.rm=TRUE), '/', nrow(selector)*ncol(selector), 
            ' values (in ',  sum(rowAnys(selector)), '/', nrow(x), 
            ' features of ', sum(colAnys(selector)), '/', ncol(x), ' samples)')
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
normimpute <- function(x, selector = is.na(x), mean = 0){
    x[selector] <- rnorm(
        length(x[selector]), mean = mean, sd = sd(x[!selector]))
    x
}


#' Impute from half-normal distribution around 0
#' 
#' @param x          NA-containing numeric vector
#' @param selector   which values to impute
#' @param mean       which mean to impute around
#' @param ref        reference (\code{translate})
#' @param pos        position (\code{translate})
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
#' x <- abs(rnorm(1e5)); x[idx] <- NA
#' dt4 <- data.table(value = translate(x), distr = 'translate')
#'
#' require(ggplot2)
#' ggplot(rbind(dt1,dt2,dt3, dt4), aes(x=value, fill=distr)) +
#' geom_density(alpha=0.5)
#' @export
halfnormimpute <- function(x, selector = is.na(x)){
    x[selector] <- abs(
        rnorm(length(x[selector]), sd = 2*sd(x[!selector], na.rm = TRUE)))
    x
}


#' @rdname halfnormimpute
#' @export
zeroimpute <- function(x, selector = is.na(x)){
    x[selector] <- 0
    x
}

#' @rdname halfnormimpute
#' @export
translate <- function(
    x, ref = c(min, mean, median, max)[[1]], pos = 3*sd(x, na.rm = TRUE)
){
    assert_any_are_true(sapply(c(min, mean, median, max), identical, ref))
    shift <- ref(x, na.rm = TRUE) - pos
    x - shift
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

has_complete_nondetects <- function(object) any(rowAlls(is.na(values(object))))

has_consistent_nondetects <- function(object, ssym){
    any( vapply(split_by_svar(object, !!enquo(ssym)), 
                has_complete_nondetects, 
                logical(1)))
}


is_systematic_detect <- function(object, subgroup = subgroup){
    . <- NULL
    subgroup <- enquo(subgroup)
    split_by_svar(object, !!subgroup) %>%
    lapply(function(x) rowAlls(is.na(values(x)))) %>%
    Reduce("|", .)
}

is_random_detect <- function(object, subgroup = subgroup){
    subgroup <- enquo(subgroup)
    rowAnys(is.na(values(object))) & !is_systematic_detect(object, !!subgroup)
}


is_full_detect <- function(object){
    rowAlls(!is.na(values(object)))
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
#' @param subgroup subgroup symbol
#' @return  \code{NULL}
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute=FALSE, plot = FALSE)
#' venn_detects(object, subgroup)
#' @export
venn_detects <- function(object, subgroup){
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
#' @param subgroup     subgroup svar
#' @param fun       imputation function
#' @param plot      TRUE or FALSE
#' @param verbose   TRUE or FALSE
#' @param ...       passed to `fun`
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute = FALSE, plot = FALSE)
#' impute_systematic_nondetects(object)
#' @export
impute_systematic_nondetects <- function(object, subgroup = subgroup,
    fun = halfnormimpute, plot = TRUE, verbose = TRUE, ...
){
# Process
    absent <- replicated <- systematic <- value <- NULL
    subgroup <- enquo(subgroup)
    groupvar <- as_name(subgroup)
# Filter
    object %<>% filter_exprs_replicated_in_some_subgroup(
                    subgroupvar = groupvar, verbose = verbose)
# Impute
    dt  <-  sumexp_to_long_dt(object, svars = groupvar)
    dt[, absent     := all(is.na(value)),    by = c('feature_id', groupvar)]
    dt[, replicated := sum(!is.na(value))>1, by = c('feature_id', groupvar)]
    dt[, systematic := any(absent) & any(replicated), by = 'feature_id']
    dt[, is_imputed := systematic & absent]
    dt[, value := fun(value, is_imputed, ...), by='feature_id']
# Update object
    ff <- fnames(object)
    ss <- snames(object)
    values(object) <- dt2exprs(dt)[ff, ss]
    assays(object)$is_imputed <- dt2mat(data.table::dcast(
                    dt, feature_id ~ sample_id, value.var = 'is_imputed'))
    fdata(object)$imputed <- rowAnys(assays(object)$is_imputed)
# Plot
    nrowimputed <- sum(rowAnys(is_imputed(object)))
    ncolimputed <- sum(colAnys(is_imputed(object)))
    if (verbose & nrowimputed>0)  message('\t\tImpute systematic nondetects ', 
        'for ', nrowimputed, '/', nrow(object), ' features ', 
        'in ',  ncolimputed, '/', ncol(object), ' samples')
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
    idx <- is.na(values(object))
    values(object)[idx] <- 0
    order <- hclust(dist(values(object)))$order
    values(object)[idx] <- NA
    object %<>% extract(order, )
    object
}

detect_order_features <- function(object, subgroup){
    subgroup <- enquo(subgroup)
    x <- object
    values(x)[is_imputed(x)] <- NA
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
#' @param subgroup      subgroup var (sym)
#' @param fill       fill var (sym)
#' @param na_imputes whether to NA imputes prior to plottin (TRUE/FALSE)g
#' @param ...        for backward compatibilty
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
#' plot_summarized_detections(object, Group)
#' plot_detections(object, Group)
#' @export
plot_detections <- function(
    object, subgroup = subgroup, fill = !!enquo(subgroup)
){
# Process
    . <- detection <- feature_id <- sample_id <- value <- NULL
    subgroup <- enquo(subgroup);         fill     <- enquo(fill)
    groupvar <- as_name(subgroup);    fillstr  <- as_name(fill)
# Reorder samples
    sdata(object)[[groupvar]] %<>% factor()
    object %<>% extract(, order(sdata(.)[[groupvar]]))
# Reorder/block features
    object %<>% detect_order_features(!!subgroup)
    y <- object; values(y)[is_imputed(y)] <- NA
    nfull       <- sum(is_full_detect(y))
    nsystematic <- sum(is_systematic_detect(y, subgroup=!!subgroup))
    nrandom     <- sum(is_random_detect(y, subgroup=!!subgroup))
# Melt
    plotdt  <-  sumexp_to_long_dt(object, svars = c(groupvar, fillstr))
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
    guides(alpha = 'none')

}

#==============================================================================
#
#                           plot_summarized_detections
#
#==============================================================================


get_subgroup_combinations <- function(object, subgroupvar){
    . <- type <- NULL
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
plot_summarized_detections <- function(
    object, subgroup = subgroup, fill = !!enquo(subgroup), na_imputes = TRUE){
# Assert
    . <- value <- NULL
    assert_is_all_of(object, "SummarizedExperiment")
    subgroup <- enquo(subgroup)
    if (quo_is_null(subgroup))  return(ggplot() + geom_blank())
    groupvar <- as_name(subgroup)
    fill <- enquo(fill);     fillstr <- as_name(fill)
    assert_is_subset(groupvar, svars(object))
    assert_is_subset(fillstr,  svars(object))
    xmin <- xmax <- ymin <- ymax <- nfeature <- quantified <- NULL
# Prepare
    object %<>% filter_samples(!is.na(!!subgroup), verbose=TRUE)
    values(object) %<>% zero_to_na()  #### TODO fine-tune
    featuretypes <- get_subgroup_combinations(object, groupvar)
    dt <- sumexp_to_long_dt(object, svars = groupvar)
    if (na_imputes) if ('is_imputed' %in% names(dt))  dt[is_imputed==TRUE,
                                                        value := NA]
    dt %<>% extract(, .(quantified   = as.numeric(any(!is.na(value)))),
                    by = c(groupvar, 'feature_id'))
    dt %<>% dcast.data.table(
        as.formula(paste0('feature_id ~ ', groupvar)),value.var='quantified')
    dt %<>% merge(featuretypes, by = setdiff(names(featuretypes), 'type'))
    dt %<>% extract(,.(nfeature=.N),by='type')
    dt %<>% merge(featuretypes,by='type')
    dt[, ymax := cumsum(nfeature)]
    dt[, ymin := c(0,ymax[-.N])]
    dt %<>% melt.data.table(id.vars = c('type', 'nfeature', 'ymin', 'ymax'),
                            variable.name = groupvar, value.name='quantified')
    dt$quantified %<>% as.factor()
    nsampledt <- sdt(object)[, .N, by=groupvar] %>% # preserves
                set_names(c(groupvar, 'xmax'))                # factor order!
    setorderv(nsampledt, groupvar)
    nsampledt[, xmax := cumsum(xmax)]; nsampledt[, xmin := c(0, xmax[-.N])]
    dt %<>% merge(nsampledt, by = groupvar)
# Plot
    colors <- make_colors(slevels(object, fillstr))
    ggplot(dt) + geom_rect(aes( xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                fill=!!fill, alpha=quantified)) +
                geom_segment(aes(x=xmin, xend=xmax, y = ymax, yend=ymax)) +
                geom_segment(aes(x=xmin, xend=xmax, y = ymin, yend=ymin)) +
                geom_segment(aes(x=xmax, xend=xmax, y = ymin, yend=ymax)) +
                geom_segment(aes(x=xmin, xend=xmin, y = ymin, yend=ymax)) +
                theme_minimal() + xlab('Samples') + ylab('Features') +
                theme(panel.grid = element_blank()) + guides(alpha = 'none') +
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
#' @param subgroup subgroup (sym)
#' @param xbiplot biplot x axis. Default pca1 (symbol)
#' @param ybiplot biplot y axis. Default pca2 (symbol)
#' @param ... aesthetic mappings
#' @return ggplot object
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, impute = FALSE, pca = TRUE, plot = FALSE)
#' explore_imputations(object, subgroup=subgroup)
#' explore_transformations(object, subgroup=subgroup)
#' @export
explore_imputations <- function(
    object, subgroup, xbiplot = pca1, ybiplot = pca2, ...
){
    subgroup <- enquo(subgroup)
    imputed  <- impute_systematic_nondetects(object, plot=FALSE)
    zeroed <- impute_systematic_nondetects(
                object, subgroup = !!subgroup, fun = zeroimpute, plot = FALSE)
    legend <- gglegend(biplot(object))

    do_plot_sample_detections <- function(obj, ...){
        plot_detections(obj, ...) + guides(color = 'none', fill = 'none')
    }

    do_biplot <- function(obj, ...){
        biplot(obj, x = !!enquo(xbiplot), y = !!enquo(ybiplot), 
                color = !!subgroup,  nloadings = 0,...) +
        guides(color = 'none', fill = 'none') +
        ggtitle(NULL)
    }

    do_plot_sample_densities <- function(obj, ...){
        plot_sample_densities(obj, ...) +
        guides(color = 'none', fill = 'none') +
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

