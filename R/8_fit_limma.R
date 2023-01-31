
#==============================================================================
#
#                       create_design
#                           single_subgroup
#                           are_factor
#                           singlelevel
#                           multilevel
#                               nlevels
#
#==============================================================================

single_subgroup <- function(object){
    assert_is_subset('subgroup', svars(object))
    length(unique(object$subgroup))==1
}


are_factor <- function(df) vapply(df, is.factor, logical(1))

nlevels <- function(object, svar){
    if (!svar %in% svars(object))  return(0)
    length(unique(object[[svar]]))
}

singlelevel <- function(object, svar)   nlevels(object, svar) ==1
multilevel  <- function(object, svar)   nlevels(object, svar) > 1

#' @rdname default_formula
#' @export
default_subgroupvar <- function(object){
    if ('subgroup' %in% svars(object))  'subgroup' else NULL
}

#' Create default formula
#' @param object SummarizedExperiment
#' @param subgroupvar string
#' @param fit 'limma', 'lm', 'lme', 'lmer'
#' @return formula
#' @examples 
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <-.read_metabolon(file)
#' default_subgroupvar(object)
#' default_formula(object)
#' @export 
default_formula <- function(
    object, subgroupvar = default_subgroupvar(object), contrasts = NULL
){
    formula <- if (is.null(subgroupvar))               '~1'
          else if (!subgroupvar %in% svars(object))    '~1'
          else if (singlelevel(object, subgroupvar))   '~1'
          else if (!is.null(contrasts))                sprintf('~0 + %s', subgroupvar)
          else                                         sprintf(   '~ %s', subgroupvar)
    formula %<>% as.formula()
    formula
}

character2factor <- function(x)  if (is.character(x)) factor(x) else x


#' Create design
#'
#'  Create design matrix  for statistical analysis
#' @param object       SummarizedExperiment or data.frame
#' @param subgroupvar  subgroup svar
#' @param formula      formula with svars
#' @param drop         whether to drop predictor names
#' @param verbose      whether to message
#' @param ...          required to s3ify
#' @return design matrix
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot = FALSE)
#' unique(create_design(object))
#'
#' object$subgroup <- 'billing19'
#' unique(create_design(object, drop = TRUE))
#'
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' unique(create_design(object))
#' create_design(object, formula= ~ 0 + SampleGroup + Sex + T2D + age + bmi)
#' object$subgroup <- 'atkin18'
#' unique(create_design(object))
#' @export
create_design <- function(object, ...) UseMethod('create_design')


#' @rdname create_design
#' @export
create_design.SummarizedExperiment <- function(
    object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula     = default_formula(object, subgroupvar),
    drop        = varlevels_dont_clash(object, all.vars(formula)), 
    verbose     = FALSE, 
    ...
){
    create_design.data.frame(sdata(object), 
                            subgroupvar = subgroupvar, 
                            formula     = formula,
                            drop        = drop,
                            verbose     = verbose)
}

#' @rdname create_design
#' @export
create_design.data.frame <- function(
    object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL,
    formula     = default_formula(object, subgroupvar),
    drop        = varlevels_dont_clash(object, all.vars(formula)), 
    verbose     = FALSE, 
    ...
){
# Assert
    assert_is_subset(all.vars(formula), names(object))
    . <- NULL
# Ensure that subgroup vector is a factor to preserve order of levels
    for (var in all.vars(formula)){
        if (is.character(object[[var]])) object[[var]] %<>% factor() }
# Create design matrix
    if (verbose)   message('\t\tDesign: ', formula2str(formula))
    myDesign <- model.matrix(formula, data=object)
    colnames(myDesign) %<>% stri_replace_first_fixed('(Intercept)', 'Intercept')
    is_factor_var <- function(x, object) is.factor(object[[x]])
    if (drop){
        for (predictor in all.vars(formula)){
            if (is.factor(object[[predictor]]))  colnames(myDesign) %<>% 
                        stri_replace_first_fixed(predictor, '') }
            # Fails for e.g. T2D = YES/NO: a meaningless column "YES" is created
            # For other cases it works wonderfully, so I keep it for now.
            # If it gives too many issues, roll back to doing the dropping only
            # for "subgroup" levels:
            #colnames(myDesign) %<>% gsub('subgroup', '', ., fixed=TRUE)
    }
# Return
    return(myDesign)
}


#=============================================================================
#
#               contrast_coefs
#                   contrast_subgroup_cols
#                   contrast_subgroup_rows
#
#==============================================================================


#' Row/Col contrasts
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup svar
#' @return  matrix
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' subgroup_matrix(object, subgroupvar = 'Group')
#' contrast_subgroup_cols(object, subgroupvar = 'Group')
#' contrast_subgroup_rows(object, subgroupvar = 'Group')
#' @export
contrast_subgroup_cols <- function(object, subgroupvar){
    subgroupmat <- subgroup_matrix(object, subgroupvar)
    if (is_scalar(subgroupmat))  return(subgroupmat)
    if (ncol(subgroupmat)==1) return(matrix(, ncol=0, nrow=nrow(subgroupmat)))
    colcontrasts <- matrix(  sprintf('%s-%s',    # no space: as lm(contr.sdiff)
                                    subgroupmat[, -1],
                                    subgroupmat[, -ncol(subgroupmat)]),
                            nrow = nrow(subgroupmat),
                            ncol = ncol(subgroupmat)-1)
    rownames(colcontrasts) <- rownames(subgroupmat)
    colnames(colcontrasts) <- sprintf('%s-%s',   # no space: as lm(contr.sdiff)
                            colnames(subgroupmat)[-1],
                            colnames(subgroupmat)[-ncol(subgroupmat)])
    colcontrasts
}


#' @rdname contrast_subgroup_cols
#' @export
contrast_subgroup_rows <- function(object, subgroupvar){
    subgroupmat <- subgroup_matrix(object, subgroupvar)
    if (nrow(subgroupmat)==1) return(matrix(, nrow=0, ncol=ncol(subgroupmat)))
    rowcontrasts <- matrix(  sprintf('%s-%s',  # no space: as lm(contr.sdiff)
                                    subgroupmat[-nrow(subgroupmat), ],
                                    subgroupmat[-1, ]),
                            nrow = nrow(subgroupmat)-1,
                            ncol = ncol(subgroupmat))
    colnames(rowcontrasts) <- colnames(subgroupmat)
    rownames(rowcontrasts) <- sprintf('%s-%s', # no space: as lm(contr.sdiff)
                            rownames(subgroupmat)[-nrow(subgroupmat)],
                            rownames(subgroupmat)[-1])
    rowcontrasts
}


# contrast_coefs <- function(object, formula){
#     subgroupvar <- all.vars(formula)[1]
#     design <- create_design(object, formula = formula)
#     if (ncol(design)==1){
#         list(matrix(colnames(design), nrow=1, ncol=1), 
#             matrix(nrow=0, ncol=0))
#     } else if (all(design[, 1]==1)){
#         list(colnames(design)[-1][seq_len(nlevels(object, subgroupvar)-1)], 
#             matrix(nrow=0, ncol=0))
#     } else {
#         list(contrast_subgroup_cols(object, subgroupvar),
#             contrast_subgroup_rows( object, subgroupvar)) }
# }

contrast_coefs <- function(object, formula){
    subgroupvar <- all.vars(formula)[1]
    design <- create_design(object, formula = formula)
    if (ncol(design)==1){  
        colnames(design)
    } else if (all(design[, 1]==1)){ 
        colnames(design)[-1][seq_len(nlevels(object, subgroupvar)-1)]
    } else {
        c(contrast_subgroup_cols(object, subgroupvar), 
          contrast_subgroup_rows( object, subgroupvar))
    }
}


#==============================================================================
#
#                            fit_limma
#
#==============================================================================

contrvec2mat  <- function(contrasts)  matrix(
                    contrasts, nrow=1, dimnames=list("", contrasts))

contrmat2list <- function(contrasts)  list(colcontrasts = contrasts)

vectorize_contrasts <- function(contrasts){
    unname(unlist(lapply(contrasts, function(x) na.exclude(c(t(x))))))
}

add_fdr <- function(fitres){
    . <- NULL
    fdr <- fitres %>% extract(, stri_startswith_fixed(
                                    names(.), paste0('p', FITSEP)), with=FALSE)
    fdr[] %<>% lapply(p.adjust, method='fdr')
    names(fdr) %<>% stri_replace_first_fixed(
                paste0('p', FITSEP), paste0('fdr', FITSEP))
    fitres %<>% cbind(fdr)
    fitres
}

reset_fitres <- function(object, fit){
    . <- NULL
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_a_string(fit)
    fitcols <- fvars(object) %>% extract(stri_detect_fixed(., fit))
    fdata(object)[fitcols] <- NULL
    metadata(object)[[fit]] <- NULL
    object
}

#' Fit results separator
#' @examples
#' FITSEP
#' @export
FITSEP <- '~'

# object: SumExp
# fitres: data.table(p.contr1, p.contr2, effect.contr1, effect.contr2)
# stat:  'p', 'effect', 'fdr', 't'
# fit:   'limma', 'wilcoxon'
merge_fitres <- function(object, fitres, fit, statistic=NULL){
    . <- NULL
    fitresdt <- data.table::copy(fitres)   # dont change in original
    firstcols <- intersect(c('feature_id', 'Intercept'), names(fitresdt))
    fitresdt %<>% extract(,c(firstcols, 
                            sort(setdiff(names(.), firstcols))), with = FALSE)
    if (!is.null(statistic)) names(fitresdt)[-1] %<>% paste0(statistic,FITSEP,.)
    names(fitresdt)[-1] %<>% paste0(FITSEP, fit)
    object %<>% merge_fdt(fitresdt)
    object
}

mat2fdt <- function(mat)  mat2dt(mat, 'feature_id')


#' Fit model and test for differential expression
#'
#' @param object       SummarizedExperiment
#' @param subgroupvar  subgroup variable
#' @param formula      modeling formula
#' @param design       design matrix
#' @param coefs        NULL or character vector: model coefs to test
#' @param contrasts    NULL or character vector: coefficient contrasts to test
#' \itemize{
#' \item{c("t1-t0", "t2-t1", "t3-t2")}
#' \item{matrix(c("WT.t1-WT.t0", "WT.t2-WT.t1", "WT.t3-WT.t2"), \cr
#'      c("KD.t1-KD.t0", "KD.t2-KD.t1", "KD.t3-KD.t2"), nrow=2, byrow=TRUE)}
#' \item{list(matrix(c("WT.t1-WT.t0", "WT.t2-WT.t1", "WT.t3-WT.t2"), \cr
#'      c("KD.t1-KD.t0", "KD.t2-KD.t1", "KD.t3-KD.t2"), nrow=2, byrow=TRUE), \cr
#'      matrix(c("KD.t0-WT.t0", "KD.t1-WT.t1", "KD.t2-WT.t2", "KD.t3-WT.t3"),\cr
#'      nrow=1, byrow=TRUE))}}
#' @param block       block svar (or NULL)
#' @param weightvar   NULL or name of weight matrix in assays(object)
#' @param statvars  character vector: subset of c('effect', 'p', 'fdr', 't')
#' @param sep       string: pvar separator  ("~" in "p~t2~limma")
#' @param suffix    string: pvar suffix ("limma" in "p~t2~limma")
#' @param verbose   whether to msg
#' @param plot      whether to plot
#' @return Updated SummarizedExperiment
#' @examples
#' # classical: lm & limma
#'     require(magrittr)
#'     file <- download_data('atkin18.somascan.adat')
#'     object <- read_somascan(file, plot = FALSE, fit = NULL)
#'     object %<>% fit_limma()
#'     object %<>% fit_lm()
#'    #plot_contrast_venn(testmat(object, coef = 't2'))
#'     
#' # blocked: limma, lme, lmer
#'     object %<>% fit_limma(subgroupvar = 'subgroup', block = 'Subject_ID')
#'     object %<>% fit_lme(  subgroupvar = 'subgroup', block = 'Subject_ID')
#'    #object %<>% fit_lmer( subgroupvar = 'subgroup', block = 'Subject_ID') # slow
#'    #plot_contrast_venn(testmat(object, coef = 't3', fit = c('limma', 'lme')))
#'     
#' # flexible: limma contrasts
#'     object %<>% fit_limma(formula = ~subgroup,   block = 'Subject_ID', coefs = 't1') 
#'     object %<>% fit_limma(formula = ~0+subgroup, block = 'Subject_ID', contrasts = 't1-t0')
#'         # flexible, but only approximate
#'         # stat.ethz.ch/pipermail/bioconductor/2014-February/057682.html
#' 
#' # alternative coding: 
#'     stats::contrasts(object$subgroup) <- MASS::contr.sdif(levels(object$subgroup))
#'     object %<>% fit_limma(subgroupvar = 'subgroup', block = 'Subject_ID') # backward difs
#'     stats::contrasts(object$subgroup) <-stats::contr.treatment(levels(object$subgroup))
#'     object %<>% fit_limma(subgroupvar = 'subgroup', block = 'Subject_ID') # baseline difs
#'     
#' # non-parametric: wilcoxon
#'     object %<>% fit_limma(   subgroupvar = 'subgroup', block = 'Subject_ID')
#'     object %<>% fit_wilcoxon(subgroupvar = 'subgroup', block = 'Subject_ID')
#' @export
fit_limma <- function(
    object, 
    subgroupvar  = if ('subgroup' %in% svars(object))  'subgroup' else NULL,
    contrasts    = NULL,
    formula      = default_formula(object, subgroupvar, contrasts),
    drop         = varlevels_dont_clash(object, all.vars(formula)),
    design       = create_design(object, formula = formula, drop = drop),
    coefs        = if (is.null(contrasts))  colnames(design)     else NULL,
    block        = NULL,
    weightvar    = if ('weights' %in% assayNames(object)) 'weights' else NULL,
    statvars     = c('effect', 'p', 'fdr'),
    sep          = FITSEP,
    suffix       = paste0(sep, 'limma'),
    verbose      = TRUE, 
    plot         = FALSE
){
    limmadt <- .fit_limma(
        object       = object,        subgroupvar  = subgroupvar,
        contrasts    = contrasts,     formula      = formula, 
        drop         = drop,          design       = design, 
        coefs        = coefs,         block        = block,
        weightvar    = weightvar,     statvars     = statvars,
        sep          = sep,           suffix       = suffix,
        verbose      = verbose)
    object %<>% reset_fitres('limma')
    object %<>% merge_fdt(limmadt)
        #fdata(object)$F.limma   <- limmares$F
        #fdata(object)$F.p.limma <- limmares$F.p
    if (plot)  print(plot_volcano(object, fit = 'limma')) 
    object
}



#' Are varlevels unique
#' 
#' @param object SummarizedExperiment or data.table
#' @return TRUE or FALSE
#' @examples 
#' require(magrittr)
#' object1 <- expand.grid(genome = c('WT', 'MUT'), treat = c('control', 'drug'))
#' object2 <- expand.grid(mutant = c('YES', 'NO'), treated = c('YES', 'NO'))
#' varlevels_dont_clash(object1)
#' varlevels_dont_clash(object2)
#' @export
varlevels_dont_clash <- function(object, ...)  UseMethod('varlevels_dont_clash')

#' @rdname varlevels_dont_clash
#' @export
varlevels_dont_clash.data.table <- function(
    object, vars = names(object)
){
    object                         %>% 
    extract(, vars, with = FALSE)  %>%
    lapply(factor)                 %>% 
    lapply(levels)                 %>% 
    unlist()                       %>% 
    duplicated()                   %>% 
    any()                          %>%
    magrittr::not()
}

#' @rdname varlevels_dont_clash
#' @export
varlevels_dont_clash.SummarizedExperiment <- function(
    object, vars = svars(object)
){
    varlevels_dont_clash.data.table(sdt(object), vars)
}


#' @rdname fit_limma
#' @export
.fit_limma <- function(
    object, 
    subgroupvar  = if ('subgroup' %in% svars(object)) 'subgroup'  else  NULL,
    contrasts    = NULL,
    formula      = default_formula(object, subgroupvar, contrasts),
    drop         = varlevels_dont_clash(object, all.vars(formula)),
    design       = create_design(object, formula = formula, drop = drop),
    coefs        = if (is.null(contrasts))  colnames(design) else NULL,
    block        = NULL, 
    weightvar    = if ('weights' %in% assayNames(object)) 'weights' else NULL, 
    statvars     = c('effect', 'p', 'fdr'),
    sep          = FITSEP,
    suffix       = paste0(sep, 'limma'),
    verbose      = TRUE, 
    plot         = FALSE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_formula(formula)
    assert_is_subset(all.vars(formula), svars(object))
    assert_is_a_bool(drop)
    assert_is_matrix(design)
    if (!is.null(block))      assert_is_subset(block, svars(object))
    if (!is.null(weightvar))  assert_scalar_subset(weightvar, assayNames(object))
    assert_is_subset(statvars, c('effect', 'p', 'fdr', 't'))
# Design/contrasts/block/weights
    . <- NULL
    if (verbose)  message('\t\tlmFit(', formula2str(formula),
        if(is.null(block))     '' else paste0(' | ',block),
        if(is.null(weightvar)) '' else paste0(', weights = assays(object)$', 
                                            weightvar), ')')
    if (!is.null(block)){
        assert_is_subset(block, svars(object))
        blockvar <- block
        block <- sdata(object)[[block]]
        if (is.null(metadata(object)$dupcor)){
            if (verbose)  message('\t\t\t\tdupcor `', blockvar, '`')
            metadata(object)$dupcor <- duplicateCorrelation(
                values(object), design = design, block = block
            )$consensus.correlation }}
    design %<>% extract(intersect(snames(object), rownames(.)), , drop = FALSE) # required in mae
    exprmat <-  values(object)[, rownames(design)]
    weightmat <- if (is.null(weightvar)){ NULL 
            } else {assert_is_a_string(weightvar)
                    assert_is_subset(weightvar, assayNames(object))
                    assays(object)[[weightvar]][, rownames(design)] }
# Fit
    limmafit <- suppressWarnings(lmFit(
                object = exprmat, design = design, 
                block = block, correlation = metadata(object)$dupcor,
                weights = weightmat))
# Effect
    if (is.null(contrasts)){  limmafit %<>% contrasts.fit(coefficients = coefs) 
    } else {                  limmafit %<>% contrasts.fit(contrasts = makeContrasts(
                                                  contrasts = contrasts, levels = design)) }
    limmadt <- data.table(feature_id = rownames(limmafit))
    if ('effect' %in% statvars){
        dt0 <- data.table(limmafit$coefficients)
        names(dt0) %<>% paste0('effect', sep, ., suffix)
        limmadt %<>% cbind(data.table(dt0)) 
    }
# p/t/fdr
    if (!all(limmafit$df.residual==0)){
        limmafit %<>% eBayes()
        if ('p' %in% statvars){ 
            dt0 <- data.table(limmafit$p.value)
            names(dt0) %<>% paste0('p', sep, ., suffix)
            limmadt %<>% cbind(dt0)  
        } 
        if ('t' %in% statvars){ 
            dt0 <- data.table(limmafit$t)
            names(dt0) %<>% paste0('t', sep, ., suffix)
            limmadt %<>% cbind(dt0)  
        } 
        if ('fdr' %in% statvars){
            dt0 <- data.table(apply(limmafit$p.value, 2, p.adjust, 'fdr') )
            names(dt0) %<>% paste0('fdr', sep, .,suffix)
            limmadt %<>% cbind(dt0)  
        } 
    }
# Return
    if (verbose)  message_df('\t\t\t%s',  summarize_fit(limmadt, fit = 'limma'))
    limmadt

}

old_summarize_fit <- function(object, fit = fits(object)){
    . <- NULL
    downdt <- colSums(down(object, fit = fit)) %>% data.table(coef = names(.), ndown = .)
    downdt %<>% tidyr::separate(
                    col = .data$coef, into = c('contrast', 'fit'), sep = FITSEP)
    
    updt <- colSums(up(object)) %>% data.table(coef = names(.), nup = .)
    updt %<>% tidyr::separate(
                    col = .data$coef, into = c('contrast', 'fit'), sep = FITSEP)
    
    sumdt <- merge(downdt, updt, by = c('fit', 'contrast'))
    sumdt$contrast %<>% factor()
    sumdt$contrast %<>% pull_level('Intercept')
    setorderv(sumdt, c('fit', 'contrast'))
    sumdt %<>% extract(fit, on = 'fit')
    sumdt %<>% extract(coefs(object), on = 'contrast')
    sumdt
}

pull_level <- function(x, lev){
    assert_is_factor(x)
    if (lev %in% levels(x))  x %<>% 
        factor(levels = c(lev, setdiff(levels(x), lev)))
    x
}

#' Summarize fit
#' @param fdt  fdt(object)
#' @param fit  'limma', 'lme', 'lm', 'lme', 'wilcoxon' or NULL
#' @param coefs string vector
#' @return data.table(contrast, nup, ndown)
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, impute = TRUE, plot = FALSE)
#' object %<>% extract(!is_consistent_nondetect(.), )
#' object %<>% fit_limma()
#' object %<>% fit_lm()
#' summarize_fit(fdt(object), coefs = c('t1', 't2', 't3'))
#' @export
summarize_fit <- function(fdt, fit = NULL, coefs = NULL){
    assert_is_data.table(fdt)
    cols <- names(fdt) %>% extract(stri_detect_fixed(., FITSEP))
    fdt %<>% extract(, c('feature_id', cols), with = FALSE)
    
    longdt <- fdt %>% melt.data.table(id.vars = 'feature_id')
    longdt[, statistic := split_extract_fixed(variable, FITSEP, 1) %>% factor(unique(.))]
    longdt[,  contrast := split_extract_fixed(variable, FITSEP, 2) %>% factor(unique(.))]
    longdt[,       fit := split_extract_fixed(variable, FITSEP, 3) %>% factor(unique(.))]
    longdt[, variable := NULL]
    
    
    sumdt <- dcast.data.table(longdt, feature_id + contrast + fit ~ statistic, value.var = 'value')
    sumdt <- sumdt[, .(
        downfdr = sum(effect < 0  & fdr < 0.05, na.rm = TRUE), 
        upfdr   = sum(effect > 0  & fdr < 0.05, na.rm = TRUE),
        downp   = sum(effect < 0  &   p < 0.05, na.rm = TRUE), 
        upp     = sum(effect > 0  &   p < 0.05, na.rm = TRUE)), by = c('fit', 'contrast') ]
    if (!is.null(fit)){
        idx <- sumdt$fit %in% fit
        sumdt %<>% extract(idx)
    }
    if (!is.null(coefs)){
        sumdt <- sumdt[contrast %in% coefs]
    }
    sumdt
}

#' Plot fit summary
#' @param sumdt data.table
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, impute = TRUE, plot = FALSE)
#' object %<>% extract(!is_consistent_nondetect(.), )
#' object %<>% fit_lm()
#' object %<>% fit_limma(block = 'SUB')
#' sumdt <- summarize_fit(fdt(object), coefs = c('t1', 't2', 't3'))
#' plot_fit_summary(sumdt)
#' @export
plot_fit_summary <- function(sumdt, nrow = NULL, ncol = NULL, order = FALSE){
    if (order){
        sumdt <- sumdt[order(downfdr+upfdr, downp+upp)]
        sumdt[, contrast := factor(contrast, unique(contrast))]
    }
    ggplot(sumdt) + facet_wrap(vars(fit), nrow = nrow, ncol = ncol) + 
    geom_col(aes(y = contrast, x = -downp),   fill = 'firebrick',   alpha = 0.3) +
    geom_col(aes(y = contrast, x =    upp),   fill = 'forestgreen', alpha = 0.3) + 
    geom_col(aes(y = contrast, x = -downfdr), fill = 'firebrick',   alpha = 1) +
    geom_col(aes(y = contrast, x =    upfdr), fill = 'forestgreen', alpha = 1) + 
    geom_text(data = sumdt[  downp>0], aes(y = contrast, x = -max(downp), label = paste0(downp, ' | ', downfdr) ), hjust = +1) + 
    geom_text(data = sumdt[    upp>0], aes(y = contrast, x =    max(upp), label = paste0(upfdr, ' | ', upp) ), hjust = 0) + 
    xlab('count') + 
    ylab(NULL) + 
    xlim(c(-max(sumdt$downp)-100, max(sumdt$upp)+100)) + 
    #scale_x_continuous(n.breaks = 20) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


setna <- function(dt, value){
   for (j in seq_len(ncol(dt)))  set(dt, which(is.na(dt[[j]])), j, value)
    dt
}


#' formula to string
#' @param formula formula
#' @return string
#' @examples 
#' formula2str(~0+subgroup)
#' @export
formula2str <- function(formula)  Reduce(paste, deparse(formula))


