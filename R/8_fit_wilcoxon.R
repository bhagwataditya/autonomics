#==============================================================================
#
#                           fit_wilcoxon
#
#==============================================================================

.wilcoxon_onesample <- function(
    dt, subgroupvar = NULL, subgrouplevels = NULL, block = NULL, verbose = TRUE
){
    . <- value <- NULL
    if (verbose)  message('\t\t\twilcox.test(x = value)')
    suppressWarnings(dt[, 
        .(  p      = wilcox.test(x = value, y = NULL, paired = FALSE)$p.value, 
            t      = wilcox.test(x = value, y = NULL, paired = FALSE)$statistic,
            effect = mean(value, na.rm=TRUE)),
        by = 'feature_id'])
}


.wilcoxon_unpaired <- function(
    dt, subgroupvar, subgrouplevels, block = NULL, verbose = TRUE
){
    . <- value <- NULL
    xx <- subgrouplevels[[1]]
    yy <- subgrouplevels[[2]]
    if (verbose)  message('\t\t\twilcox.test(x = ', xx, ', y = ', yy, ')')
    suppressWarnings(dt[, 
        .(  p = wilcox.test(x      = value[get(subgroupvar)==xx],
                            y      = value[get(subgroupvar)==yy], 
                            paired = FALSE)$p.value, 
            t = wilcox.test(x      = value[get(subgroupvar)==xx],
                            y      = value[get(subgroupvar)==yy], 
                            paired = FALSE)$statistic, 
            effect = if (is.null(yy)){  
                            mean(value[get(subgroupvar)==xx], na.rm=TRUE) 
                } else {    mean(value[get(subgroupvar)==yy], na.rm=TRUE) - 
                            mean(value[get(subgroupvar)==xx], na.rm=TRUE) }),
        by = 'feature_id'])
}


.wilcoxon_paired <- function(
    dt, subgroupvar, subgrouplevels, block, verbose = TRUE
){
    . <- NULL
    dt <- data.table::dcast(dt, as.formula(sprintf(
            'feature_id + feature_name + %s ~ %s', block, subgroupvar)),
            value.var = 'value')
    xx <- subgrouplevels[[1]]
    yy <- subgrouplevels[[2]]
    if (verbose)  message(
        "\t\t\twilcox.test(", "x = ", xx, ", y = ", yy, ", paired = TRUE) - ", 
        "pair on '", block, "'")
    suppressWarnings(dt[!is.na(get(xx)) & !is.na(get(yy)), 
        .(  p = wilcox.test(x = get(xx), y = get(yy), paired = TRUE)$p.value, 
            t = wilcox.test(x = get(xx), y = get(yy), paired = TRUE)$statistic, 
            effect = mean(get(yy) - get(xx), na.rm=TRUE)),
        by = 'feature_id'])
}


.wilcoxon <- function(contrastdef, dt, subgroupvar, block, verbose){
    subgrouplevels <- stri_split_regex(contrastdef, pattern = '[ ]*[-][ ]*')
    subgrouplevels %<>% unlist()
    subgrouplevels %<>% rev()
    assert_is_subset(length(subgrouplevels), c(1,2))
    fun <- if (length(subgrouplevels)==1){ .wilcoxon_onesample
        } else if (is.null(block)){        .wilcoxon_unpaired
        } else {                           .wilcoxon_paired   }
    resdt <- fun(dt, 
                subgroupvar = subgroupvar, 
                subgrouplevels = subgrouplevels, 
                block = block, verbose = verbose)
    data.table::setnames(resdt, c('p', 't', 'effect'), 
                        sprintf('%s.%s', c('p', 't', 'effect'), contrastdef))
    resdt
}

#' @export
#' @rdname fit_limma
fit_wilcoxon <- function(
    object,
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL,
    formula = default_formula(object, subgroupvar, fit = 'wilcoxon'),
    contrastdefs = contrast_coefs(object, formula=formula),
    block = NULL, weightvar = NULL, verbose = TRUE, plot = FALSE
){
# fit
    . <- NULL
    dt <- sumexp_to_long_dt(object, svars = c(subgroupvar, block))
    if (verbose)  message('\t\tWilcoxon')
    fitres <- lapply(vectorize_contrastdefs(contrastdefs), .wilcoxon, 
                    dt, subgroupvar, block, verbose)
    fitres %<>% Reduce(function(x,y) merge(x,y, by='feature_id', all=TRUE), .)
    fitres %<>% merge(data.table(fdata(object))[, 'feature_id', drop = FALSE], 
                        ., by = 'feature_id', all.x = TRUE)
    fitres %<>% add_fdr()
    object %<>% merge_fitres(fitres, fit='wilcoxon')
# extract
    extract_quantity <- function(quantity, fitres){
        quantitydot <- paste0(quantity, '.')
        quantitymat <- fitres[, stri_startswith_fixed(
                        names(fitres), quantitydot), with=FALSE]
        quantitymat %<>% as.matrix()
        rownames(quantitymat) <- fitres$feature_id
        colnames(quantitymat) %<>% stri_replace_first_fixed(quantitydot, '')
        quantitymat }
    
    fitres <- mapply(extract_quantity, c('effect', 'fdr', 'p', 't'), 
                    MoreArgs=list(fitres=fitres), SIMPLIFY=FALSE)
# wrap
    metadata(object)$wilcoxon <- do.call(abind::abind, c(fitres, along = 3))
    names(dimnames(metadata(object)$wilcoxon)) <-
                                        c('feature', 'contrast', 'quantity')
                                      #  c('feature', subgroupvar, 'quantity')
# Return
    if (plot)  print(plot_volcano(object, fit='wilcoxon')) 
    if (verbose)  message_df('\t\t\t%s', summarize_fit(object,'wilcoxon'))
    object
}

