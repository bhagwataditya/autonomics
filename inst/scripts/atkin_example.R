# Load
devtools::document()

# Read
mfile <- download_data('atkin.metabolon.xlsx')
pfile <- download_data('atkin.somascan.adat')
mobj <- read_metabolon(mfile)
pobj <- read_somascan( pfile)

# Glucose
table(mobj$Subject, mobj$Time)
mobj %<>% filter_samples(Subject != 'C01')
plot_exprs(mobj['glucose', ], block = 'Subject', geom = 'point', coefs = NULL, shape = 'Diabetes', size = 'Diabetes') + 
    scale_size_manual(values = c(Control = 3, T2DM = 3)) + 
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          panel.grid.major.x = element_line(color = 'gray80'))

# t1 - t0
contrasts(mobj$subgroup) <- MASS::contr.sdif(levels(mobj$subgroup))
contrasts(pobj$subgroup) <- MASS::contr.sdif(levels(pobj$subgroup))

object <- mobj
subgroupvar <- c('subgroup', 'Diabetes')
block <- 'Subject'
object$Diabetes %<>% factor()

svar_coefs <- function(object, svar)  colnames(contrasts(object[[svar]]))


#' Fit across/within/between
#' @param object     SummarizedExperiment
#' @param groupvars  c(svar1, svar2)
#' @param block      svar
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file)
#' object$Diabetes %<>% factor()
#' fdt1 <- .fit_xxx_across( object, groupvars = c('subgroup', 'Diabetes'), block = 'Subject')
#' fdt2 <- .fit_xxx_within( object, groupvars = c('subgroup', 'Diabetes'), block = 'Subject')
#' fdt3 <- .fit_xxx_within( object, groupvars = c('Diabetes', 'subgroup'), block = 'Subject')
#' fdt4 <- .fit_xxx_between(object, groupvars = c('Diabetes', 'subgroup'), block = 'Subject')
#' object %<>% fit_across_within_between(groupvars = c('subgroup', 'Diabetes'), block = 'Subject')
#' @export
fit_xxx_across <- function(object, groupvars, block){
# Formula
    formula <- paste0(groupvars, collapse = ' + ')
    formula %<>% paste0('~ ', .)
    formula %<>% as.formula()
# Difference Contrasts
    coefs <- character(0)
    for (groupvar in groupvars){
        levels <- base::levels(object[[groupvar]])
        contrastmat <- MASS::contr.sdif(levels)
        contrasts(object[[groupvar]]) <- contrastmat
        coefs %<>% c(colnames(contrastmat))
    }
    object %<>% fit_limma(formula = formula, block = block,  coefs = coefs)
# Treatment Contrasts
    coefs <- character(0)
    for (groupvar in groupvars){
        levels <- base::levels(object[[groupvar]])
        contrastmat <- contr.treatment(levels)
        colnames(contrastmat) %<>% paste0('-', levels[1])
        contrasts(object[[groupvar]]) <- contrastmat
        coefs %<>% c(colnames(contrastmat)[-1]) # first contrast already present
    }
    fdt1 <- .fit_limma(object, formula = formula, block = block,  coefs = coefs)
}

#' @rdname dot-fit_xxx_across
#' @export
.fit_xxx_within <- function(object, groupvars, block){
    assert_is_of_length(groupvars, 2)
    Avar <- groupvars[1]
    Bvar <- groupvars[2]
    Alevels <- levels(object[[Avar]])
    Blevels <- levels(object[[Bvar]])
    Acoefs <- svar_coefs(object, Avar)
    Bcoefs <- svar_coefs(object, Bvar)
    
    formula <- as.formula(sprintf('~ %s / %s', Bvar, Avar))
    coefs <- sprintf('%s:%s', rep(Blevels, times = length(Acoefs)), 
                              rep(Acoefs,  each  = length(Blevels)))
    .fit_limma(object, formula = formula, block = block,  coefs = coefs)
}

#' @rdname dot-fit_xxx_across
#' @export
.fit_xxx_between <- function(object, groupvars, block){
    assert_is_of_length(groupvars, 2)
    Avar <- groupvars[[1]]
    Bvar <- groupvars[[2]]
    Alevels <- levels(object[[Avar]])
    Blevels <- levels(object[[Bvar]])
    Acoefs <- svar_coefs(object, Avar)
    Bcoefs <- svar_coefs(object, Bvar)
    
    formula <- as.formula(sprintf('~ %s * %s', Avar, Bvar))
    coefs <- sprintf('%s:%s', rep(Acoefs, times = length(Bcoefs)), 
                              rep(Bcoefs, each  = length(Acoefs)))
    .fit_limma(object, formula = formula, block = block, coefs = coefs)
}


fit_across_within_between <- function(object, groupvars, block){
# Subsequent Differences    
    object %<>% reset_fit('limma')
    Avar <- groupvars[[1]]
    Bvar <- groupvars[[2]]
    Alevels <- levels(object[[Avar]])
    Blevels <- levels(object[[Bvar]])
    contrasts(object[[Avar]]) <- MASS::contr.sdif(Alevels)
    contrasts(object[[Bvar]]) <- MASS::contr.sdif(Blevels)
    fdt1a <- .fit_xxx_across( object,     groupvars,  block = block)
    fdt2a <- .fit_xxx_within( object,     groupvars,  block = block)
    fdt3a <- .fit_xxx_within( object, rev(groupvars), block = block)
    fdt4a <- .fit_xxx_between(object,     groupvars,  block = block)
    sdifvars <- union(names(fdt1), names(fdt2), names(fdt3), names(fdt4))
# Baseline Differences
    contrasts(object[[Avar]]) <- contr.treatment(Alevels)
    contrasts(object[[Bvar]]) <- contr.treatment(Blevels)
    colnames(contrasts(object[[Avar]])) %<>% paste0('-', Alevels[1])
    colnames(contrasts(object[[Bvar]])) %<>% paste0('-', Blevels[1])
    fdt1b <- .fit_xxx_across( object,     groupvars,  block = block)
    fdt2b <- .fit_xxx_within( object,     groupvars,  block = block)
    fdt3b <- .fit_xxx_within( object, rev(groupvars), block = block)
    fdt4b <- .fit_xxx_between(object,     groupvars,  block = block)
    cols5 <- c('feature_id', setdiff(names(fdt5), sdifvars))
    cols6 <- c('feature_id', setdiff(names(fdt6), sdifvars))
    cols7 <- c('feature_id', setdiff(names(fdt7), sdifvars))
    cols8 <- c('feature_id', setdiff(names(fdt8), sdifvars))
    fdt5 %<>% extract(, cols5, with = FALSE)
    fdt6 %<>% extract(, cols6, with = FALSE)
    fdt7 %<>% extract(, cols7, with = FALSE)
    fdt8 %<>% extract(, cols8, with = FALSE)
# Return
    object %<>% merge_fdt(fdt1) # across
    object %<>% merge_fdt(fdt5)
    object %<>% merge_fdt(fdt2)
    object %<>% merge_fdt(fdt6)
    object %<>% merge_fdt(fdt3)
    object %<>% merge_fdt(fdt4)
    object %<>% merge_fdt(fdt7)
    object %<>% merge_fdt(fdt8)
    object
}






pdt2 <- .fit_limma(pobj, formula = ~ Time + Diabetes, block = 'Subject', coefs = c('t1-t0', 't2-t1', 't3-t2', 'T2DM'))
mobj %<>% fit_limma(     formula = ~ Time + Diabetes, block = 'Subject', coefs = c('t1-t0', 't2-t1', 't3-t2', 'T2DM'))
pobj %<>% fit_limma(     formula = ~ Time + Diabetes, block = 'Subject', coefs = c('t1-t0', 't2-t1', 't3-t2', 'T2DM'))
mobj %<>% order_on_p(coef = 't1-t0')
pobj %<>% order_on_p(coef = 't1-t0')
plot_exprs(mobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_exprs(pobj[2,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes', facet = 'EntrezGeneSymbol') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# t2 - t1
mobj %<>% order_on_p(coef = 't2-t1')
pobj %<>% order_on_p(coef = 't2-t1')
fnames(mobj)[1] <- fdt(mobj)$feature_id[1] <- 'cAMP'
plot_exprs(mobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_exprs(pobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes', facet = 'EntrezGeneSymbol') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# t3 - t2
mobj %<>% order_on_p(coef = 't3-t2')
pobj %<>% order_on_p(coef = 't3-t2')
plot_exprs(mobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_exprs(pobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes', facet = 'EntrezGeneSymbol') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# T2DM
mobj %<>% order_on_p(coef = 'T2DM')
pobj %<>% order_on_p(coef = 'T2DM')
plot_exprs(mobj[2,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_exprs(pobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes', facet = 'EntrezGeneSymbol') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# D : t1 - t0
mobj %<>% fit_limma(formula = ~ Diabetes/Time, block = 'Subject')
pobj %<>% fit_limma(formula = ~ Diabetes/Time, block = 'Subject')
mobj %<>% order_on_p(coef = 'T2DM:t1-t0')
pobj %<>% order_on_p(coef = 'T2DM:t1-t0')
plot_exprs(mobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_exprs(pobj[2,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# C : t1 - t0
mobj %<>% order_on_p(coef = 'Control:t1-t0')
pobj %<>% order_on_p(coef = 'Control:t1-t0')
plot_exprs(mobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_exprs(pobj[2,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# D : t2 - t1
mobj %<>% order_on_p(coef = 'T2DM:t2-t1')
pobj %<>% order_on_p(coef = 'T2DM:t2-t1')
plot_exprs(mobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_exprs(pobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes', facet = 'TargetFullName') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# D : t2 - t1
mobj %<>% order_on_p(coef = 'Control:t2-t1')
pobj %<>% order_on_p(coef = 'Control:t2-t1')
plot_exprs(mobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_exprs(pobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes', facet = 'TargetFullName') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# D : t3 - t2
mobj %<>% order_on_p(coef = 'T2DM:t3-t2')
pobj %<>% order_on_p(coef = 'T2DM:t3-t2')
plot_exprs(mobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_exprs(pobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes', facet = 'TargetFullName') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# C : t3 - t2
mobj %<>% order_on_p(coef = 'Control:t3-t2')
pobj %<>% order_on_p(coef = 'Control:t3-t2')
plot_exprs(mobj[1,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
plot_exprs(pobj[2,],
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes', facet = 'TargetFullName') + theme(legend.position = "none") + ylab(NULL) + theme(axis.text.x = element_blank(), axis.text.y = element_blank())

# Explore
mobj %<>% order_on_p(coef = 'T2DM:t3-t2')
pobj %<>% order_on_p(coef = 'T2DM:t3-t2')
plot_exprs(mobj[1:6, ], n = 6, 
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes', 
           size  = 'Diabetes', nrow = 2, ncol = 3) + scale_size_manual(values = c(Control = 3, T2DM = 3))
plot_exprs(pobj[c(1,2,4,5,6,7)], n = 6,
           coef  = NULL,
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes', 
           size  = 'Diabetes',
           ncol  = 3,
           facet = c('TargetFullName', 'EntrezGeneSymbol')) + scale_size_manual(values = c(C = 3, D = 3))


pobj %<>% order_on_p()
mobj %<>% order_on_p()
plot_exprs(mobj[c(29,1,2,3), ], 
           coefs = NULL, 
           block = 'Subject', 
           geom  = 'point', 
           shape = 'Diabetes', 
           size  = 'Diabetes', 
           nrow  = 1, ncol = 4) + scale_size_manual(values = c(Control = 3, T2DM = 3))


mobj %>% fit_limma()
mobj$Time
mobj$Diabetes

# 
mobj

# t0: strong dip and T2D/control difference
contrasts(mobj$subgroup) <- MASS::contr.sdif(levels(mobj$subgroup))
mobj %<>% fit_limma(formula = ~ subgroup/Diabetes, block = 'Subject')
mobj %<>% order_on_p(    coefs = c('t1-t0', 't0:T2DM'), combiner = '&')
mobj[1, ] %>% plot_exprs(coefs = c('t1-t0', 't0:T2DM'), geom = 'point', block = 'Subject', shape = 'Diabetes', size = 'Diabetes') +  scale_size_manual(values = c(Control = 2, T2DM = 3))

# t3: strong recovery AND difference between diabetics and controls
contrasts(mobj$subgroup) <- MASS::contr.sdif(levels(mobj$subgroup))
mobj %<>% fit_limma(formula = ~ subgroup/Diabetes, block = 'Subject')
mobj %<>% order_on_p(    coefs = c('t3-t2', 't3:T2DM'), combiner = '&')
mobj[1, ] %>% plot_exprs(coefs = c('t3-t2', 't3:T2DM'), geom = 'point', block = 'Subject', shape = 'Diabetes', size = 'Diabetes') +  scale_size_manual(values = c(Control = 2, T2DM = 3))

# 
contrasts(mobj$subgroup) <- MASS::contr.sdif(levels(mobj$subgroup))
mobj %<>% fit_limma(formula = ~ subgroup*Diabetes, block = 'Subject')
mobj %<>% order_on_p(coefs = c('t3-t2', 't3-t2:T2DM'), combiner = '&')
mobj[1, ] %>% plot_exprs(coef = NULL, geom = 'point', block = 'Subject', shape = 'Diabetes', facet = 'Diabetes', ncol = 2)


mobj |> fit_limma(formula = ~ Diabetes/subgroup, block = 'Subject') |>
        plot_exprs(coef = 'Diabetes', geom = 'point', block = 'Subject', shape = 'Diabetes', size = 'Diabetes') + 
        scale_size_manual(values = c(Control = 3, T2DM = 3))
    