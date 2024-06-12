
prep_model_summary <- function(
    object, 
       fit = rep( fits(object), times = length(coefs(object))), 
      coef = rep(coefs(object),  each = length( fits(object))),
      nmax = NULL
){
# Abstract
    color <- coef1 <- coef2 <- y <- n <- xlev <- ylev <- NULL
    object %<>% abstract_fit(fit = fit, coef = coef)
    
     dnames <- paste0( fit, '.', coef )
      upmat <- downmat <- matrix(0, nrow = length(dnames), 
                                    ncol = length(dnames), 
                                dimnames = list( dnames, dnames))
    for (row in dnames){
    for (col in dnames){
        x <- abstractvec(object, fit = split_extract_fixed(row, '.', 1), coef = split_extract_fixed(row, '.', 2))
        y <- abstractvec(object, fit = split_extract_fixed(col, '.', 1), coef = split_extract_fixed(col, '.', 2))
        downmat[ row, col ] <- length(intersect( fnames(object[x=='down', ]),  fnames(object[y=='down', ]) ))
          upmat[ row, col ] <- length(intersect( fnames(object[x=='up',   ]),  fnames(object[y=='up',   ]) ))
    } }

# Prepare
    updt <- mat2dt(  upmat, idvar = 'xlev')
    dndt <- mat2dt(downmat, idvar = 'xlev')
    updt %<>% melt.data.table(id.vars = 'xlev', variable.name = 'ylev', value.name = 'n')
    dndt %<>% melt.data.table(id.vars = 'xlev', variable.name = 'ylev', value.name = 'n')
    updt[, xlev := factor(xlev, unique(xlev))]
    dndt[, xlev := factor(xlev, unique(xlev))]
    updt[, ylev := factor(ylev, unique(ylev))]
    dndt[, ylev := factor(ylev, unique(ylev))]
    dndt[, x := as.numeric(xlev)];  dndt[, y := as.numeric(ylev)]
    updt[, x := as.numeric(xlev)];  updt[, y := as.numeric(ylev)]
    dndt <- dndt[y>=x]  #dndt[, pct := n / nrow(object)*100 ]
    updt <- updt[x>=y]  #updt[, pct := n / nrow(object)*100 ]
    dndt[, dir := 'down']
    updt[, dir := 'up']
    plotdt <- rbind(dndt, updt)
# Color
    if (is.null(nmax))  nmax <- plotdt[, max(abs(n))]
    downcolors <- colorRampPalette(c('white', '#dc3e00'))(nmax+1)
      upcolors <- colorRampPalette(c('white', '#00BBDC'))(nmax+1)
    colordt <- rbind( data.table( dir = 'down', n = 0:nmax,  color = downcolors),
                      data.table( dir = 'up',   n = 0:nmax,  color = upcolors) )
    plotdt %<>% merge(colordt, by = c('dir', 'n'), sort = FALSE)
    plotdt
}

#' Plot model summary
#' @param object  SummarizedExperiment
#' @param fit     character vector (values : 'limma', 'lm', 'Lme', 'lmer', 'wilcoxon' )
#' @param coef    character vector 
#' @param title   character(1)
#' @param sub     character(1)
#' @param xlab    character(1)
#' @param ylab    character(1)
#' @param nmax    integer(1)
#' @param cex.sub number
#' @examples
#' # Read
#'     file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'     object <- read_metabolon(file)
#' # Block
#'     plot_model_summary(fit_wilcoxon(object),                    sub =           '~Time (wilcoxon)', nmax = 10)
#'     plot_model_summary(fit_wilcoxon(object, block = 'Subject'), sub = '~Time | Subject (wilcoxon)', nmax = 10)
#'     plot_model_summary(      fit_lm(object),                    sub =           '~Time (lm)',       nmax = 10)
#'     plot_model_summary(      fit_lm(object, block = 'Subject'), sub = '~Time | Subject (lm)',       nmax = 10)
#'     plot_model_summary(   fit_limma(object),                    sub =           '~Time (limma)',    nmax = 10)
#'     plot_model_summary(   fit_limma(object, block = 'Subject'), sub = '~Time | Subject (limma)',    nmax = 10)
#' # Engines
#'     plot_model_summary(fit_wilcoxon(object, block = 'Subject'), sub = 'Wilcoxon : ~Time | Subject', nmax = 10)
#'     plot_model_summary(      fit_lm(object, block = 'Subject'), sub =       'Lm : ~Time | Subject', nmax = 10)
#'     plot_model_summary(   fit_limma(object, block = 'Subject'), sub =    'Limma : ~Time | Subject', nmax = 10)
#'     plot_model_summary(     fit_lme(object, block = 'Subject'), sub =      'Lme : ~Time | Subject', nmax = 10)
#' @export
plot_model_summary <- function(
     object, 
        fit = rep( fits(object), times = length(coefs(object))), 
       coef = rep(coefs(object),  each = length( fits(object))),
      title = '', 
        sub = '', 
       xlab = '', 
       ylab = '', 
       nmax = NULL, 
    cex.sub = 1.4
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(fit, fits(object))
    assert_is_a_string(title)
    assert_is_a_string(sub)
    assert_is_a_string(xlab)
    assert_is_a_string(ylab)
    if (!is.null(nmax))  assert_is_a_number(nmax)
    x <- y <- color <- n <- xlev <- ylev <- NULL
# Prepare
    plotdt <- prep_model_summary(object, fit = fit, coef = coef, nmax = nmax)
    dd <- 0.30
    k <- length(levels(plotdt$xlev))
# by fails
    # https://github.com/Rdatatable/data.table/issues/6177
        # plotdt <- data.table( x = c(1, 2), 
        #                       y = c(1, 2), 
        #                   color = c('#DAF5F9', '#00BBDC') )
        # plot(1, type = 'n', xlab = '', ylab = '', xlim = c(0,3), ylim = c(-3,0), axes = 0)  # axes = 1
        # for (j in 1:nrow(plotdt)) plotdt[j][ , polygon(x = c(x-1,x,x    ), y = -c(y-1,y-1,y  ), col = color )          ]
        #                           plotdt[    , polygon(x = c(x-1,x,x    ), y = -c(y-1,y-1,y  ), col = color ) , by = x ]
    # By-based code
        # plot(1, type = 'n', xlab = '', ylab = '', xlim = c(-0.3,k), ylim = c(-k,0.3), axes = 0, cex.sub = cex.sub)
        # plotdt[, i := 1:.N]
        # plotdt[ x==y & dir=='up'     ][ , polygon(x = c(x-1,x,x    ), y = -c(y-1,y-1,y  ), col = color                ) , by = i ]
        # plotdt[ x==y & dir=='down'   ][ , polygon(x = c(x-1,x,x-1  ), y = -c(y-1,y,  y  ), col = color                ) , by = i ]
        # plotdt[ x!=y                 ][ , polygon(x = c(x-1,x,x,x-1), y = -c(y-1,y-1,y,y), col = color                ) , by = i ]
        # plotdt[ dir=='up'            ][ ,    text(x = x-dd,           y = -y+1-dd,         col = 'black', label = n   ) , by = i ]
        # plotdt[ dir=='down'          ][ ,    text(x = x-1+dd,         y = -y+dd,           col = 'black', label = n   ) , by = i ]
        # plotdt[ y==1 & dir == 'up'   ][ ,    text(x = x-0.5,          y = 0.2,                            label = xlev) , by = i ]
        # plotdt[ x==1 & dir == 'down' ][ ,    text(x = -0.2,           y = -y+0.5,                         label = ylev) , by = i ]
# for works
    theta <- 0.015
    plot(1, type = 'n', xlab = '', ylab = '', xlim = c(-0.3,k), ylim = c(-k,0.3), axes = 0)  # axes = 1
    subdt <- plotdt[x==y & dir=='up'    ]; for (j in 1:nrow(subdt)) subdt[j][ , polygon(x = c(x-1,x,x,  x-1), y = -c(y-1,y-1,y, y-1), col = color, lwd = 1  )]
    subdt <- plotdt[x==y & dir=='down'  ]; for (j in 1:nrow(subdt)) subdt[j][ , polygon(x = c(x-1,x,x-1,x-1), y = -c(y-1,y,  y, y-1), col = color, lwd = 1  )]
    subdt <- plotdt[x!=y];                 for (j in 1:nrow(subdt)) subdt[j][ , polygon(x = c(x-1,x,x,x-1),   y = -c(y-1,y-1,y, y  ), col = 'white', lwd = 0.5)]
    subdt <- plotdt[dir=='up'           ]; for (j in 1:nrow(subdt)) subdt[j][ , text(x = x-dd, y = -y+1-dd, col = 'black', label = n)]
    subdt <- plotdt[dir=='down'         ]; for (j in 1:nrow(subdt)) subdt[j][ , text(x = x-1+dd, y = -y+dd, col = 'black', label = n)]
    subdt <- plotdt[y==1 & dir == 'up'  ]; for (j in 1:nrow(subdt)) subdt[j][ , text(x = x-0.5,  y = 0.2,    label = xlev) ]
    subdt <- plotdt[x==1 & dir == 'down']; for (j in 1:nrow(subdt)) subdt[j][ , text(x = -0.05,  y = -y+0.5, label = ylev, adj = 1)]
    plotdt[, graphics::segments(x0 = theta,   x1 = max(x),       y0 =  0,       y1 = 0,               col = '#00bbdc', lwd = 4  )]
    plotdt[, graphics::segments(x0 = max(x),  x1 = max(x),       y0 =  0,       y1 = -(max(y)-theta), col = '#00bbdc', lwd = 4  )]
    plotdt[, graphics::segments(x0 = 0,       x1 = 0,            y0 = -theta,   y1 = - max(y),        col = '#dc3e00', lwd = 4  )]
    plotdt[, graphics::segments(x0 = 0,       x1 = max(x)-theta, y0 = -max(y),  y1 = - max(y),        col = '#dc3e00', lwd = 4  )]
    plotdt[, graphics::segments(x0 = 0,       x1 = max(x),       y0 =  0,       y1 = -(max(y)),       col = 'white',   lwd = 0.1)]
    plotdt[, graphics::segments(x0 = 0+theta, x1 = max(x),       y0 =  0,       y1 = -(max(y)-theta), col = '#00bbdc', lwd = 4  )]
    plotdt[, graphics::segments(x0 = 0,       x1 = max(x)-theta, y0 = -theta,   y1 = -(max(y)),       col = '#dc3e00', lwd = 4  )]
  # plotdt[, polygon(x = c(0, max(x), max(x)), y = c(0,         0, -(max(y))), border = '#00bbdc', lwd = 1)]
  # plotdt[, polygon(x = c(0, max(x), 0     ), y = c(0, -(max(y)), -(max(y))), border = '#dc3e00', lwd = 1)]
    title(main = title, sub = sub, line = 1, xlab = xlab, ylab = ylab, cex.sub = cex.sub)
}


aa <- function(){
    file <- download_data('atkin.metabolon.xlsx')
    object <- read_metabolon(file)
    object %<>% fit_lm(      block = 'Subject', codingfun = contr.treatment.explicit, sep = '.')
    object %<>% fit_limma(   block = 'Subject', codingfun = contr.treatment.explicit, sep = '.')
    object %<>% fit_lme(     block = 'Subject', codingfun = contr.treatment.explicit, sep = '.')
    object %<>% fit_wilcoxon(block = 'Subject', codingfun = contr.treatment.explicit, sep = '.')
    summarize_fit(object)
    plot_model_summary(object, fit = 'limma', title = 'limma : ~ Time | Subject')
    
}