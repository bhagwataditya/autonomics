

find_medoid <- function(subdt){
    mat <- dcast.data.table(subdt, feature_id ~ sample_id, value.var = 'value')
    mat %<>% dt2mat()
    mat %<>% extract(!matrixStats::rowAnyNAs(.), , drop = FALSE)
    if (nrow(mat) <= 3)  return(as.character(subdt$feature_id[1]))
    names(which.medoid(t(mat)))
}


#' Feature type
#' @param object   SummarizedExperiment
#' @param formula  model formula
#' @param design   design matrix
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object %<>% fit_limma(block = 'Subject')
#' ftype(object)
#' @export
ftype <- function(
       object, 
      formula = default_formula(object), 
         drop = varlevels_dont_clash(object, all.vars(formula)),
          fit = fits(object)[1],
    codingfun = if (fit == 'wilcoxon')  contr.treatment.explicit  else  contr.treatment
){
# Assert
    assert_is_valid_sumexp(object)
    assert_valid_formula(formula, object)
    assert_is_a_bool(drop)
    assert_is_function(codingfun)
    assert_scalar_subset(fit, fits(object))
# Predict
       xmat <- X(object, formula = formula, drop = drop, codingfun = codingfun)
    betamat <- beta(object, fit = fit)
    y <- xmat %*% betamat
# Type
    typedt <- mat2sdt(y)
    typedt %<>% melt.data.table(id.vars = 'sample_id', variable.name = 'feature_id', value.name = 'y')
    typedt[ , x    := seq(0, .N-1), by = 'feature_id']
    typedt[ , xlab := rownames(xmat),  by = 'feature_id' ]
    typedt[ , y := data.table::frank(y, ties.method = 'dense'), by = 'feature_id']
    typedt[ , y := y - rep(y[1], length(y)), by = 'feature_id']
    typedt[ , ystr := sprintf('%3d', y)]
    typedt[ , type := paste0(ystr, collapse = ''), by = 'feature_id']
# Medoid
    valuedt <- sumexp_to_longdt(object, fvars = NULL, svars = NULL)[, c('feature_id', 'sample_id', 'value'), with = FALSE]
    valuedt %<>% merge(typedt[, .SD[1, c('type')], by = 'feature_id'], by = 'feature_id')
    valuedt[ , typemedoid := find_medoid(.SD), by = 'type']
    valuedt %<>% extract(, c('type', 'typemedoid'), with = FALSE)
    valuedt %<>% unique()
    typedt %<>% merge(valuedt, by = 'type', sort = FALSE)
    typedt
# Count
    ndt <- typedt[ , .SD[1], by = 'feature_id'][, .( n = .N ), by = 'type']
    typedt %<>% merge(ndt, by = 'type', sort = FALSE)
    typedt %<>% extract(, c('feature_id', 'type', 'typemedoid', 'x', 'y', 'n', 'xlab'), with = FALSE)
    typedt %<>% unique()
    typedt %<>% extract(order(type, feature_id, x, y))
# Plot & Return
    print(plot_contrast_types(typedt))
    print(plot_contrast_trajectories(typedt))
    object %<>% merge_fdt(unique(typedt[, c('feature_id', 'type', 'typemedoid'), with = FALSE]))
    object
}


plot_contrast_types <- function(typedt){
    plotdt <- unique(typedt[, .SD, .SDcols = setdiff(names(typedt), 'feature_id')])
    plotdt %<>% extract(rev(order(n)))
    plotdt[, facet := paste0(typemedoid, '\n', n) ]
    plotdt[, facet := factor(facet, unique(facet))]
    
    ggplot(plotdt[n>1]) + theme_bw() + theme(panel.grid = element_blank()) + facet_wrap(vars(facet)) + 
    geom_line(aes(x = xlab, y = y, group = type, color = as.factor(n), linewidth = n)) + 
    guides(color = 'none', linewidth = 'none')
}


plot_contrast_trajectories <- function(typedt){
    plotdt <- unique(typedt[, .SD, .SDcols = setdiff(names(typedt), 'feature_id')])
    plotdt <- plotdt[ , .( x0 = x[-.N], 
                              x1 = x[-1], 
                              y0 = y[-.N], 
                              y1 = y[-1], 
                               n = n[-1], 
                            xlab = xlab[-1]), by = 'type' ]
    plotdt <- plotdt[ , .(n = sum(n)) , by = c('x0', 'x1', 'y0', 'y1') ]
    plotdt %<>% extract(rev(order(n)))
    plotdt[                              , alpha := 0.15 ]
    plotdt[ n >= length(unique(n)) * 2/3 , alpha := 1   ]
    plotdt[,     dy  := y1-y0]
    plotdt[,     dx  := x1-x0]
    plotdt[,      dz := sqrt(dx^2+dy^2) ]
    plotdt[,     rad := atan(dy/dx)     ]
    plotdt[, degrees := 360/(2*pi)*rad  ]
    plotdt[,   xtext := x0 + cos(rad)*dz*1/3 ]
    plotdt[,   ytext := y0 + sin(rad)*dz*1/3 ]
    plotdt[,   xarrowstart := x0 + cos(rad)*dz*9/10 ]
    plotdt[,   yarrowstart := y0 + sin(rad)*dz*9/10 ]
    plotdt %<>% extract(order(n))

    ggplot() + theme_bw() + theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    geom_text(    data = plotdt, aes(x = xtext, y = ytext,       label = n,           color = as.factor(x1), alpha = alpha), vjust = -1) + 
    geom_segment( data = plotdt[alpha <1], aes(x = x0, xend = x1,          y = y0, yend = y1,          color = as.factor(x1), alpha = alpha, linewidth = n)) + 
    geom_segment( data = plotdt[alpha==1], aes(x = x0, xend = xarrowstart, y = y0, yend = yarrowstart, color = as.factor(x1), alpha = alpha, linewidth = n)) + 
    geom_segment( data = plotdt[alpha==1], aes(x = x0, xend = x1,          y = y0, yend = y1,          color = as.factor(x1), alpha = alpha),  arrow = arrow(angle = 20, length = unit(0.35, 'inches'), type = 'closed')) + 
    scale_alpha_identity() + 
    guides(linewidth = 'none') + 
    scale_x_continuous(labels = unique(typedt$xlab)) + 
    xlab(NULL) + 
    ylab(assayNames(object)[1]) + 
    guides(color = 'none')
    
}




