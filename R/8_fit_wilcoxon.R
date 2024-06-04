#==============================================================================
#
#                           fit_wilcoxon
#
#==============================================================================

.wilcoxon_onesample <- function(
    dt, subgroupvar = NULL, subgrouplevels = NULL, block = NULL, verbose = TRUE
){
    . <- value <- NULL
    if (verbose)  cmessage('%swilcox.test(x = value)', spaces(22))
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
    if (verbose)  cmessage('\t\t\twilcox.test(x = %s, y = %s)', xx, yy )
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
            'feature_id + %s ~ %s', block, subgroupvar)),
            value.var = 'value')
    xx <- subgrouplevels[[1]]
    yy <- subgrouplevels[[2]]
    if (verbose)  cmessage("\t\t\twilcox.test(x = %s, y = %s, paired = TRUE) - pair on %s", xx, yy, block)
    suppressWarnings(dt[!is.na(get(xx)) & !is.na(get(yy)), 
        .(  p = wilcox.test(x = get(xx), y = get(yy), paired = TRUE)$p.value, 
            t = wilcox.test(x = get(xx), y = get(yy), paired = TRUE)$statistic, 
            effect = mean(get(yy) - get(xx), na.rm=TRUE)),
        by = 'feature_id'])
}


.wilcoxon <- function(contrastdef, dt, subgroupvar, block, sep, verbose){
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
                        paste(c('p', 't', 'effect'), contrastdef, sep = sep ))
    resdt
}

all_vars <- function(x){
    y <- all.vars(x)
    if (length(y)==0)  y <- NULL
    y
}
    
#' @export
#' @rdname fit
fit_wilcoxon <- function(
    object,
      formula = default_formula(object), 
         drop = NULL,
    codingfun = contr.treatment.explicit, # wilcox is the only one where `contr.treatment` doesnt work
       design = NULL, # only so that fit(.) works
    contrasts = NULL,
        coefs = NULL, 
        block = NULL, 
    weightvar = NULL, 
     statvars = c('effect', 'p'),
          sep = FITSEP,
       suffix = paste0(sep, 'wilcoxon'),
      verbose = TRUE, 
         plot = FALSE
){
# assert
    assert_is_valid_sumexp(object)
    assert_valid_formula(formula, object)
    subgroupvar <- all_vars(formula)[1]
    if (is.null(contrasts)){
        contrasts <- colnames(create_design(object, formula = formula, drop = TRUE, codingfun = codingfun))[-1]
    }
    assert_is_character(contrasts)
    if (!is.null(block))      assert_is_subset(block, svars(object))
    if (verbose)  cmessage('%sFeatures', spaces(14))
    object %<>% reset_fit('wilcoxon')
    obj <- object
    obj %<>% keep_replicated_features(formula, n = 1, verbose = verbose)
    # connected block filtering not required, .wilcoxon doesnt break there
# fit
    . <- NULL
    dt <- sumexp_to_longdt(obj, svars = c(subgroupvar, block))
    if (verbose)  cmessage('%sWilcoxon', spaces(14))
    fitres <- lapply(vectorize_contrasts(contrasts), .wilcoxon, 
                     dt, subgroupvar = subgroupvar, block = block, sep = sep, verbose = verbose)
    fitres %<>% Reduce(function(x, y)  merge(x, y, by = 'feature_id', all = TRUE), .)
    pattern <- sprintf('^(feature_id|%s)',  paste0(statvars, collapse = '|'))   # select statvars
    fitres <- fitres[, .SD, .SDcols = patterns(pattern) ]
    names(fitres)[-1] %<>% paste0(suffix)
    if (verbose)  message_df('\t\t\t%s', summarize_fit(fitres, fit = 'wilcoxon'))
    object %<>% merge_fit(fitres)
# extract
    extract_quantity <- function(quantity, fitres){
        quantitydot <- paste0(quantity, FITSEP)
        quantitymat <- fitres[, stri_startswith_fixed(
                        names(fitres), quantitydot), with = FALSE]
        quantitymat %<>% as.matrix()
        rownames(quantitymat) <- fitres$feature_id
        colnames(quantitymat) %<>% stri_replace_first_fixed(quantitydot, '')
        quantitymat }
# Return
    if (plot)  print(plot_volcano(object, fit = 'wilcoxon')) 
    object
}

