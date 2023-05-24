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
#' matrix(c(0, 7), nrow = 1)
#' matrix(c(0, 7), nrow = 1)  %>% zero_to_na(verbose = TRUE)
#' 
#' matrix(c(NA, 7), nrow = 1)
#' matrix(c(NA, 7), nrow = 1) %>% na_to_zero(verbose = TRUE)
#' 
#' matrix(c(NaN, 7), nrow = 1)
#' matrix(c(NaN, 7), nrow = 1)  %>% nan_to_na(verbose = TRUE)
#' 
#' matrix(c(Inf, 7), nrow = 1)
#' matrix(c(Inf, 7), nrow = 1)  %>% inf_to_na(verbose = TRUE)
#' 
#' matrix(c(-Inf, 7), nrow = 1)
#' matrix(c(-Inf, 7), nrow = 1) %>% minusinf_to_na(verbose = TRUE)
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


#' @rdname zero_to_na
#' @export
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

#' Impute
#' 
#' Impute NA values
#'
#' Imputes NA values from N(mean - 2.5 sd, 0.3 sd)
#' @param object   numeric vector, SumExp
#' @param assay    string
#' @param by       svar
#' @param shift    number: sd units
#' @param width    number: sd units
#' @param frac     fraction: fraction of available samples should be greater 
#'                           than this value for a subgroup to be called available
#' @param verbose  TRUE or FALSE
#' @param plot     TRUE or FALSE
#' @param n        number of samples to plot
#' @param palette  color vector
#' @param ...      required for s3 dispatch
#' @return numeric vector, matrix or SumExp
#' @examples
#' # Simple Design
#'    file <- download_data('fukuda20.proteingroups.txt')
#'    object <- read_maxquant_proteingroups(file)
#'    impute(values(object)[, 1], plot = TRUE)[1:3]              # vector
#'    impute(values(object),      plot = TRUE)[1:3, 1:3]         # matrix
#'    impute(object, plot = TRUE)                                # sumexp
#' # Complex Design
#'    file <- download_data('atkin18.metabolon.xlsx')
#'    object <- read_metabolon(file)
#'    invisible(impute(values(object)[1:3, 1   ]))               # vector
#'    invisible(impute(values(object)[1:3, 1:5 ]))               # matrix
#'    object %>%  filter_samples(T2D == 'Control') %>% impute()  # sumexp
#' @export
impute <- function(object, ...) UseMethod('impute')

#' @rdname impute
#' @export
impute.numeric <- function(
    object, shift = 2.5, width = 0.3, verbose = TRUE, plot = FALSE, ...
){
# Original
    count <- imputed <- NULL
    sd1    <- sd(object, na.rm = TRUE)
    mean1  <- mean(object, na.rm = TRUE)
# Imputed
    mean0 <- mean1 - shift*sd1
    sd0 <- width*sd1
    idx    <- is.na(object)
    if (verbose)  message('\tImpute ', sum(idx), ' / ', length(idx), ' values')
    n <- length(object[idx])
    object[idx] <- rnorm(n, mean = mean0, sd = sd0)
# Plot and Return
    if (plot){
        dt <- data.table(x = object, imputed = idx)
        p <- ggplot(dt) + 
             geom_density(aes(x = x, y = after_stat(count), fill = imputed))
        print(p)
    }
    object
}

#' @rdname impute
#' @export
impute.matrix <- function(
    object, 
    shift   = 2.5, 
    width   = 0.3, 
    verbose = TRUE, 
    plot    = FALSE, 
    n       = min(9, ncol(object)),  
    palette = make_colors(colnames(object)), 
    ...
){
    count <- imputed <- sample_id <- value <- NULL
    idx <- is.na(object)
    if (verbose){
        message(sprintf('\tImpute (out of %d) features per sample: ', nrow(object)))
        message_df('\t\t%s', colSums(idx[, 1:n]))
    }
    object %<>% apply(2, impute.numeric,
                 shift = shift, width = width, verbose = FALSE, plot = FALSE)
    if (plot){
        dt1 <- mat2dt(object[,  1:n], 'feature_id')
        dt2 <- mat2dt(idx[,1:n], 'feature_id')
        dt1 %<>% melt.data.table(id.vars = 'feature_id', variable.name = 'sample_id', value.name = 'value')
        dt2 %<>% melt.data.table(id.vars = 'feature_id', variable.name = 'sample_id', value.name = 'imputed')
        dt <- merge(dt1, dt2, by = c('feature_id', 'sample_id'))
        p <- ggplot(dt) + 
             geom_density(aes(x = value, y = after_stat(count), fill = sample_id, 
                              group = interaction(sample_id, imputed))) +
            scale_fill_manual(values = palette)
        print(p)
    }
    object
}

#' @rdname impute
#' @export 
impute.SummarizedExperiment <- function(
    object,
    assay    = assayNames(object)[1],
    by       = 'subgroup',
    shift    = 2.5, 
    width    = 0.3, 
    frac     = 0.5,
    verbose  = TRUE, 
    plot     = FALSE, 
    palette  = make_colors(colnames(object)), 
    n        = min(9, ncol(object)), 
    ...
){
# Assert
    assert_is_scalar(assay); assert_is_subset(assay, assayNames(object))
    assert_is_a_number(shift)
    assert_is_a_number(width)
    assert_is_a_bool(verbose)
    assert_is_a_bool(plot)
    assert_is_character(palette)
    assert_has_names(palette)
    assert_is_a_number(n)
    consistent.na <- imputed <- isNa <- isValue <- na.group <- value <- NULL
    value.group <- NULL
# Impute systematic NAs
    dt <- sumexp_to_longdt(object, assay = assay, svars = by)
    dt[, imputed := impute(value, shift = shift, width = width, 
                           verbose = FALSE, plot = FALSE), by = 'sample_id']
    dt[, isNa    :=  is.na(value)]
    dt[, isValue := !is.na(value)]
    dt[, na.group    := sum(isNa)   == .N,   by = c('feature_id', by)]
    dt[, value.group := sum(isValue) > frac*.N, by = c('feature_id', by)]
    dt[, consistent.na := na.group & any(value.group), by = 'feature_id']
    dt[consistent.na==TRUE, value := imputed]
# Update object
    mat <- dcast(dt, feature_id ~ sample_id, value.var = 'value')
    mat %<>% dt2mat()
    mat %<>% extract(rownames(object), )
    mat %<>% extract(, colnames(object))
    is_imputed(object) <- is.na(values(object))  &  !is.na(mat)
    fdt(object)$imputed <- rowAnys(is_imputed(object))
    values(object) <- mat
    if (verbose & any(is_imputed(object))){
        message(sprintf('\tImputed %d/%d features in the following groups:', 
                        sum(is_imputed_feature(object)), nrow(object)))
        message_df('\t\t%s', n_imputed_features_per_subgroup(object))
    }
# Plot/Return
    if (plot){
        #p1 <- plot_sample_densities(object, fill = 'subgroup') + guides(fill = 'none')
        #dt <- sumexp_to_longdt(object)
        # p1 <- ggplot(dt) + ggridges::geom_density_ridges(aes(
        #         x = value, y = sample_id, fill = subgroup, group = sample_id)) + 
        #         theme_bw() + ggtitle('Sample densities')
        p1 <- if (ncol(object)<=9){ plot_sample_nas(object) + guides(fill = 'none')
              } else {   plot_subgroup_nas(object) }
        p2 <- plot_sample_violins(object, fill = by) + guides(fill = 'none') + ggtitle(NULL)
        gridExtra::grid.arrange(p1, p2, nrow = 2)
    }
    object
}

# # @rdname impute
# # @export
# impute.list <- function(
#     object,
#     assay = assayNames(object[[1]])[1],
#     by = 'subgroup',
#     shift    = 2.5, 
#     width    = 0.3, 
#     frac     = 0.5,
#     verbose  = TRUE, 
#     plot     = FALSE, 
#     palette  = make_colors(colnames(object)), 
#     n        = min(9, ncol(object))
# ){
# # Assert
#     assert_is_list(object)                                            # is it a list
#     assert_all_are_true(vapply(object, is_valid_sumexp, logical(1)))  # of valid sumexps
#     assert_all_are_true(Reduce(identical, lapply(object, fdt)))       # with identical fdt
# # Impute then Align
#     object %<>% Map(impute.SummarizedExperiment, .) 
#     imputed <- object %>% lapply(fdt)
#     imputed %<>% lapply(extract2, 'imputed') 
#     imputed %<>% Reduce(`|`, .)
#     object %<>% Map(function(obj){ fdt(obj)$imputed <- imputed; obj }, .)
# # Return
#     object
# }

n_imputed_features_per_subgroup <- function(object){ 
    objlist <- split_samples(object, by = 'subgroup')
    n <- vapply(objlist, n_imputed_features, integer(1))
    n[n!=0]
}
n_imputed_features  <- function(object) sum(is_imputed_feature(object))
n_imputed_samples   <- function(object) sum(is_imputed_sample( object))
is_imputed_feature  <- function(object)     rowAnys(is_imputed(object))
is_imputed_sample   <- function(object)     colAnys(is_imputed(object))

# Plot imputation densities
# @param object SummarizedExperiment
# @return 
# @examples 
# @export
#plot_imputation_densities <- function(object, n = min(9, ncol(object))){
# Prepare
#    dt1 <- mat2dt(    values(object)[, 1:n], 'feature_id')
#    dt2 <- mat2dt(is_imputed(object)[, 1:n], 'feature_id')
#    dt1 %<>% melt.data.table(id.vars = 'feature_id', variable.name = 'sample_id', value.name = 'value')
#    dt2 %<>% melt.data.table(id.vars = 'feature_id', variable.name = 'sample_id', value.name = 'imputed')
#    dt <- merge(dt1, dt2, by = c('feature_id', 'sample_id'))
# Measured
#    p <- ggplot(dt[imputed==FALSE])
#    p <- p + geom_density(aes(x = value, y = stat(count), fill = sample_id), na.rm = TRUE)
# Imputed
#    dt %<>% extract(imputed==TRUE)
#    dt %<>% extract(, .SD[.N>2], by = c('sample_id'))
#    if (nrow(dt)>0){
#        p <- p + geom_density(aes(x = value, y = stat(count), fill = sample_id), data = dt, na.rm = TRUE) 
#    }
# Return
#    p <- p + scale_fill_manual(values = palette)
#    p
#}


# @rdname halfnormimpute
# @export
# normimpute <- function(x, selector = is.na(x), mean = 0){
#    x[selector] <- rnorm(
#        length(x[selector]), mean = mean, sd = sd(x[!selector]))
#    x
#}


# Impute from half-normal distribution around 0
# 
# @param x          NA-containing numeric vector
# @param selector   which values to impute
# @param mean       which mean to impute around
# @param ref        reference (\code{translate})
# @param pos        position (\code{translate})
# @return numeric vector of same length
# @examples
# x <- rnorm(1e5, mean = 5)
# idx <- runif(length(x))>0.9
# x[idx] <- NA
# dt0 <- data.table(x = x, method = '0.original')
# dt1 <- data.table(x =     zeroimpute(x)[idx], method = '1.zeroimpute')
# dt2 <- data.table(x =     normimpute(x)[idx], method = '2.normimpute')
# dt3 <- data.table(x = halfnormimpute(x)[idx], method = '3.halfnormimpute')
# dt <- rbindlist(list(dt0, dt1, dt2, dt3, dt4))
# ggplot(dt) + geom_density(aes(x = x, y = stat(count), group = method, fill = method), alpha = 0.5)
# @export
#halfnormimpute <- function(x, selector = is.na(x)){
#    x[selector] <- abs(
#        rnorm(length(x[selector]), sd = 2*sd(x[!selector], na.rm = TRUE)))
#    x
#}


# @rdname halfnormimpute
# @export
#zeroimpute <- function(x, selector = is.na(x)){
#    x[selector] <- 0
#    x
#}

# @rdname halfnormimpute
# @export
#translate <- function(
#    x, ref = c(min, mean, median, max)[[1]], pos = 3*sd(x, na.rm = TRUE)
#){
#    assert_any_are_true(sapply(c(min, mean, median, max), identical, ref))
#    shift <- ref(x, na.rm = TRUE) - pos
#    x - shift
#}

#=============================================================================
#
#                        split_samples
#
#==============================================================================

#' Split samples
#'
#' Split samples by svar
#' @param object  SummarizedExperiment
#' @param objlist SummarizedExperiment list
#' @param by     svar to split by (string)
#' @return  SummarizedExperiment list
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' objlist <- split_features(object, by = 'PLATFORM')
#' objlist <- split_samples(object, 'T2D')
#' objlist %<>% Map(impute, .)
#' object <- cbind_imputed(objlist)
#' @export
split_samples <- function(object, by = 'subgroup'){
    if (!by %in% svars(object))  return(list(object))
    extract_samples  <- function(level){  object %>% extract(, .[[by]] == level)  }
    Map(extract_samples, slevels(object, by))
}


#' @rdname split_samples
#' @export
cbind_imputed <- function(objlist){
    imputed <- objlist %>% lapply(fdt)
    imputed %<>% lapply(extract2, 'imputed') 
    imputed %<>% Reduce(`|`, .)
    objlist %<>% Map(function(obj){ fdt(obj)$imputed <- imputed; obj }, .)
    object <- Reduce(SummarizedExperiment::cbind, objlist)
    object
}


#' @rdname split_samples
#' @export
split_features <- function(object, by){
    if (!by %in% fvars(object))  return(list(object))
    extract_features  <- function(level)  object %>% extract(fdt(.)[[by]] == level, )
    Map(extract_features, flevels(object, by))
}

#=============================================================================
#
#                     systematic_nas
#                     random_nas
#                     no_nas
#
#=============================================================================


# all
# any
# fraction
fraction <- function(x, frac)  sum(x) >= frac*length(x)


#' Is systematic/random/full NA
#' @param object SummarizedExperiment
#' @param by     svar (string)
#' @param frac   fraction
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' head(systematic_nas(object))   # missing in some subgroups, present in others
#' head(random_nas(object))       # missing in some samples, independent of subgroup
#' head(no_nas(object))           # missing in no samples
#' @export
systematic_nas <- function(object, by = 'subgroup', frac = 0.5){
# Assert
    assert_is_valid_sumexp(object)
    assert_scalar_subset(by, svars(object))
    assert_is_fraction(frac)
    value <- nalevel <- valuelevel <- NULL
# Call
    nlevels <- length(unique(object[[by]]))
    dt <- sumexp_to_longdt(object, svars = by)
    dt %<>% extract(, .(nalevel =      all( is.na(value)),
                     valuelevel = fraction(!is.na(value), frac = frac)), by = c('feature_id', by) )
    dt %<>% extract(, .SD[any(nalevel) & any(valuelevel)] , by = by)
    y <- fnames(object) %in% as.character(dt[, feature_id])
# Return
    y
}

#' @rdname systematic_nas
#' @export
random_nas <- function(object, by = 'subgroup'){
    rowAnys(is.na(values(object))) & 
    !systematic_nas(object, by)
}

#' @rdname systematic_nas
#' @export
no_nas <- function(object){
    rowAlls(!is.na(values(object)))
}


#=============================================================================
#
#                     venn_detects
#
#=============================================================================


#' Venn detects
#'
#' Venn diagram full/consistent/random detects
#'
#' @param object  SummarizedExperiment
#' @param by      svar (string)
#' @return  \code{NULL}
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file, impute = FALSE, plot = FALSE)
#' venn_detects(object, 'subgroup')
#' @export
venn_detects <- function(object, by = 'subgroup'){
    limma::vennDiagram(as.matrix(cbind(
        consistent = systematic_nas(object, by),
        random     = random_nas(    object, by),
        full       = no_nas(        object))))
}


#==============================================================================
#
#                      cluster_order_features
#                      detect_order_features
#
#==============================================================================

cluster_order_features <- function(object){
    
    assert_is_identical_to_false(any(is.infinite(values(object))))
        # Otherwise Error in hclustfun(distfun(x)) : NA/NaN/Inf
        # https://stackoverflow.com/questions/23738900
    if (nrow(object) < 3) return(object)
    idx <- is.na(values(object))
    values(object)[idx] <- 0
    order <- hclust(dist(values(object)))$order
    values(object)[idx] <- NA
    object %<>% extract(order, )
    object
}

detect_order_features <- function(object, by){
    objlist <- object %>% split_samples(by)
    block_detected <- function(obj){
        as.character(as.numeric(!rowAlls(is.na(values(obj)))))
    }
    object <- Reduce(cbind, objlist)
    fdt(object)$natype <- objlist %>% Map(block_detected, .) %>% Reduce(paste0, .)
    ngroup <- length(objlist)
    onestring <- paste0(rep('1', ngroup), collapse = '')
    twostring <- paste0(rep('2', ngroup), collapse = '')
    fdt(object)$natype[rowAlls(!is.na(values(object)))] <- twostring

    natypes <- table(fdt(object)$natype)
    natypes %<>% sort(decreasing = TRUE)
    natypes %<>% names()
    natypes %<>% setdiff(c(onestring, twostring))
    natypes %<>% c(onestring, twostring)
    fdt(object)$natype %<>% factor(natypes)
    object %<>% extract(order(fdt(.)$natype), )
    object
}


#==============================================================================
#
#                           plot_sample_nas
#
#==============================================================================

#' @rdname plot_sample_nas
#' @export
plot_detections <- function(...){
    .Deprecated('plot_sample_nas')
    plot_sample_nas(...)
}

#' @rdname plot_sample_nas
#' @export
plot_summarized_detections <- function(...){
    .Deprecated('plot_subgroup_nas')
    plot_subgroup_nas(...)
}


#' Plot missingness per sample / subgroup
#'
#' \code{plot_sample_nas} shows systematic and random missingness 
#' (white), and full detection (bright color) at sample resolution.
#' Imputations are also shown (light color).
#'
#' \code{plot_subgroup_nas} shows systematic missingness at subgroup resolution.
#' Random missingness and full detection are shown together (bright color).
#' Imputations are also shown (light color).
#' @param object       SummarizedExperiment
#' @param by           svar (string)
#' @param fill         svar (string)
#' @param palette      color vector (names = levels, values = colors)
#' @param axis.text.y  passed to ggplot2::theme
#' @param na_imputes   TRUE or FALSE
#' @param ...          used to maintain deprecated functions
#' @return ggplot object
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file, impute = FALSE)
#' plot_subgroup_nas(object)
#' plot_sample_nas(object)
#' plot_sample_nas(impute(object))
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' plot_subgroup_nas(object, 'subgroup')
#' plot_subgroup_nas(object, 'TIME_POINT')
#' plot_sample_nas(object, 'subgroup')
#' plot_sample_nas(object, 'TIME_POINT')
#' @export
plot_sample_nas <- function(
    object, 
    by = 'subgroup', 
    fill = by, 
    palette = make_svar_palette(object, fill), 
    axis.text.y = element_blank()
){
# Process
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_subset(by, svars(object))
    assert_is_subset(fill, svars(object))
    . <- detection <- feature_id <- sample_id <- value <- NULL
# Reorder samples
    object[[by]] %<>% factor()
    object %<>% extract(, order(.[[by]]))
# Reorder/block features
    # object %<>% detect_order_features_old(by = by)
    object %<>% detect_order_features(by = by)
    y <- object; values(y)[is_imputed(y)] <- NA
    nfull       <- sum(no_nas(y))
    nconsistent <- sum(systematic_nas(y, by = by))
    nrandom     <- sum(random_nas(    y, by = by))
# Melt
    plotdt  <-  sumexp_to_longdt(object, svars = unique(c(by, fill)))
    alpha <- NULL
    plotdt[,             detection := 'detect']
    plotdt[is.na(value), detection := 'nondetect']
    if ('is_imputed' %in% SummarizedExperiment::assayNames(object)){
        plotdt[is_imputed==TRUE, detection := 'impute']}
    plotdt[, detection  := factor(detection, c('nondetect', 'impute', 'detect'))]
    plotdt[, sample_id  := factor( sample_id, unique(snames(object)))]
    plotdt[, feature_id := factor(feature_id, rev(unique(fnames(object))))]
    palette %<>% c(nondetect = "#FFFFFF")
# Plot
    p <- ggplot(plotdt) +
        geom_tile(aes(x = sample_id, y = feature_id, fill = !!sym(fill), alpha = detection)) +
        scale_fill_manual( values = palette, breaks = names(palette) %>% setdiff('nondetect')) +
        scale_alpha_manual(values = c(nondetect = 0, impute = 0.3, detect = 1)) +
        scale_x_discrete(position = 'top') +
        ylab(paste0(nrow(object), ' Features')) +
        xlab(paste0(ncol(object), ' Samples')) +
        #ggtitle(sprintf('detects: %d full, %d random, %d consistent',
        #                nfull, nrandom, nconsistent)) +
        theme_bw() +
        theme(   axis.text.y  = axis.text.y, 
                 axis.ticks.y = element_blank(),
                 axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 0),
                 legend.title = element_blank(),
                 panel.grid   = element_blank(), 
                 panel.border = element_rect(linewidth = 1.5))
# Add lines
    hlines <- cumsum(c(nfull, nrandom )) + 0.5  # cumsum(c(nfull, nrandom, nconsistent))
    fainthlines <- cumsum(rev(table(fdt(object)$natype))) + 0.5
    vlines <- table(object[[by]])
    vlines %<>% extract(-length(.))
    vlines %<>% cumsum()
    vlines %<>% add(0.5)
    p <- p + geom_hline(yintercept = hlines, size = 1)
    p <- p + geom_hline(yintercept = fainthlines, size = 0.5)
    p <- p + geom_vline(xintercept = vlines, size = 1)
    p +  guides(alpha = 'none')
}


#==============================================================================
#
#                           plot_subgroup_nas
#
#==============================================================================


get_subgroup_combinations <- function(object, by = 'subgroup'){
    . <- type <- NULL
    subgroups <- slevels(object, by)
    subgroups  %>%
        lapply(function(x) c(0,1) %>% set_names(rep(x,2))) %>%
        set_names(subgroups) %>%
        expand.grid() %>%
        data.table() %>%
        extract(rev(order(rowSums(.)))) %>%
        extract(, type := 0:(.N-1)) %>%
        extract()
}


#' @rdname plot_sample_nas
#' @export
plot_subgroup_nas <- function(
    object, by = 'subgroup', fill = by, palette = NULL, na_imputes = TRUE
){
# Assert
    . <- value <- NULL
    assert_is_all_of(object, "SummarizedExperiment")
    assert_is_subset(c(by, fill), svars(object))
    if (!is.null(palette)){
        assert_is_character(palette)
        assert_is_subset(slevels(object, fill), names(palette))
    }
    assert_is_a_bool(na_imputes)
    xmin <- xmax <- ymin <- ymax <- nfeature <- quantified <- NULL
# Prepare
    object[[by]] %<>% num2char()
    object %<>% extract(, !is.na(.[[by]]))
    values(object) %<>% zero_to_na()  #### TODO fine-tune
    featuretypes <- get_subgroup_combinations(object, by)
    dt <- sumexp_to_longdt(object, svars = by)
    if (na_imputes) if ('is_imputed' %in% names(dt))  dt[is_imputed==TRUE,
                                                        value := NA]
    dt %<>% extract(, .(quantified   = as.numeric(any(!is.na(value)))),
                    by = c(by, 'feature_id'))
    dt %<>% dcast.data.table(
        as.formula(paste0('feature_id ~ ', by)), value.var = 'quantified')
    dt %<>% merge(featuretypes, by = setdiff(names(featuretypes), 'type'))
    dt %<>% extract(,.(nfeature=.N), by = 'type')
    dt %<>% merge(featuretypes, by = 'type')
    dt[, ymax := cumsum(nfeature)]
    dt[, ymin := c(0,ymax[-.N])]
    dt %<>% melt.data.table(id.vars = c('type', 'nfeature', 'ymin', 'ymax'),
                            variable.name = by, value.name = 'quantified')
    dt$quantified %<>% as.factor()
    nsampledt <- sdt(object)[, .N, by = by] %>% # preserves
                set_names(c(by, 'xmax'))      # factor order!
    setorderv(nsampledt, by)
    nsampledt[, xmax := cumsum(xmax)]; nsampledt[, xmin := c(0, xmax[-.N])]
    dt %<>% merge(nsampledt, by = by)
# Plot
    npersubgroup <- table(object[[by]])
    xbreaks <- c(cumsum(npersubgroup)- npersubgroup/2)
    ybreaks <- c(0.01*nrow(object), 0.99*nrow(object))
    if (is.null(palette))  palette <- make_colors(slevels(object, by))
    p <- 
        ggplot(dt) + 
        geom_rect(aes( xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                        fill = !!sym(fill), alpha = quantified)) +
        scale_y_continuous(expand = c(0, 0)) + #, limits = c(0, nrow(object)))  +
        scale_x_continuous(breaks = xbreaks, position = 'top', expand = c(0,0)) + #, limits = c(0, ncol(object)))  + 
        geom_segment(aes(x = xmin, xend = xmax, y = ymax, yend = ymax)) +
        geom_segment(aes(x = xmin, xend = xmax, y = ymin, yend = ymin)) +
        geom_segment(aes(x = xmax, xend = xmax, y = ymin, yend = ymax)) +
        geom_segment(aes(x = xmin, xend = xmin, y = ymin, yend = ymax)) +
        #theme_minimal() + 
        xlab(paste0(formatC(ncol(object), format = 'd', big.mark = ' '), ' Samples')) + 
        ylab(paste0(formatC(nrow(object), format = 'd', big.mark = ' '), ' Features')) +
        theme(panel.grid      = element_blank(), 
              legend.position = 'top', 
              legend.title    = element_blank(), 
              axis.text.x     = element_text(angle = 90), 
              axis.text.y     = element_blank()) + 
        guides(alpha = 'none', fill = 'none') +
        #guides(alpha = 'none', fill = guide_legend(label.position = 'top', nrow = 1, title.hjust = 0.5)) +
        scale_fill_manual(values = palette) +
        scale_alpha_manual(values=c(`0`=0, `1`=1))
    p
}


#==============================================================================

#'@title Get/set is_imputed
#'@description Get/Set is_imputed
#'@param object SummarizedExperiment
#'@param value matrix
#'@return matrix (get) or updated object (set)
#'@examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file, plot = FALSE)
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

