#==============================================================================
#
#                           fit_wilcoxon
#
#==============================================================================

#' @export
#' @rdname fit_limma
fit_wilcoxon <- function(
    object, formula = if ('subgroup' %in% svars(object)) ~ subgroup else ~ 1, 
    subgroupvar = all.vars(formula), 
    contrastdefs = contrast_subgroups(object, create_design(object, formula)), 
    block = NULL, verbose = TRUE, plot = TRUE
){
# fit
    dt <- sumexp_to_long_dt(object, svars = c(subgroupvar, block))
    results <- lapply(vectorize_contrastdefs(contrastdefs), 
                    .fit_wilcoxon, dt, subgroupvar, block)
    results %<>% Reduce(merge, .)
# extract
    extract_quantity <- function(quantity, results){
        quantitydot <- paste0(quantity, '.')
        quantitymat <- results[, stri_startswith_fixed(
                        names(results), quantitydot), with=FALSE]
        quantitymat %<>% as.matrix()
        rownames(quantitymat) <- results$feature_id
        colnames(quantitymat) %<>% stri_replace_first_fixed(quantitydot, '')
        quantitymat }
    
    results <- mapply(extract_quantity, c('p', 'w', 'effect'), 
                    MoreArgs=list(results=results), SIMPLIFY=FALSE)
    results$fdr  <- apply(results$p, 2, p.adjust, method = 'fdr')
    results$bonf <- apply(results$p, 2, p.adjust, method = 'bonf')
# wrap
    metadata(object)$wilcoxon <- do.call(abind::abind, c(results, along = 3))
# Return
    if (plot)  print(plot_volcano(object, fit='wilcoxon')) 
                    # plot_contrastogram(object)
    if (verbose) cmessage_df('\t\t\t%s', summarize_fit(object,'wilcoxon'))
    object
}


.fit_wilcoxon <- function(contrastdef, dt, subgroupvar, block){
    if (!is.null(block))  dt <- data.table::dcast(dt, as.formula(sprintf(
                    'feature_id + feature_name + %s ~ %s', block, subgroupvar)),
                    value.var = 'value')
    terms <- unlist(stri_split_regex(contrastdef,pattern='[ ]*[-][ ]*'))
    x <- terms[[2]]
    y <- terms[[1]]
    resdt <- suppressWarnings(dt[, .(
        p = wilcox.test(  
                x = if (is.null(block)) value[get(subgroupvar)==x] else get(x),
                y = if (is.null(block)) value[get(subgroupvar)==y] else get(y), 
                paired = !is.null(block))$p.value, 
        w = wilcox.test(
                x = if (is.null(block)) value[get(subgroupvar)==x] else get(x),
                y = if (is.null(block)) value[get(subgroupvar)==y] else get(y), 
                paired = !is.null(block))$statistic, 
        effect = if (is.null(block)){
                mean(value[get(subgroupvar)==y])-
                mean(value[get(subgroupvar)==x])
            } else {
                mean(get(y) - get(x), na.rm=TRUE)
            }),
        by = 'feature_id'])
    data.table::setnames(resdt, c('p', 'w', 'effect'), 
                        sprintf('%s.%s', c('p', 'w', 'effect'), contrastdef))
    resdt
}
    

