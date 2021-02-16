#=============================================================================
#
#                   guess_sep
#                   guess_sep.SummarizedExperiment
#                   guess_sep.factor
#                   guess_sep.character
#                       has_identical_values
#                       is_max
#                           cequals
#
#=============================================================================

#' Convenient equals operator
#'
#' Performs x == y, but returns FALSE rather than NA for NA elements of x.
#' @param x numeric vector or scalar
#' @param y numeric scalar
#' @return logical vector
#' @examples
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' y <- 3
#' cequals(x, y)
#' @noRd
cequals <- function(x,y){
    result <- rep(FALSE, length(x)) %>% set_names(names(x))
    if (is.na(y)){
        result[ is.na(x)] <- TRUE
        result[!is.na(x)] <- FALSE
    } else {
        result[ is.na(x)] <- FALSE
        result[!is.na(x)] <- x[!is.na(x)] == y
    }
    result
}


#' Is maximal
#' @param x numeric vector
#' @return logical vector
#' @examples
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' is_max(x)
#' @noRd
is_max <- function(x) cequals(x, max(x, na.rm = TRUE))


#' All elements of vector are identical
#' @param x vector
#' @return TRUE or FALSE
#' @examples
#' x <- c(2,2,1,2)
#' has_identical_values(x)
#' @noRd
has_identical_values <- function(x) length(unique(x))==1

#' Guess separator
#' @param x          character vector or SummarizedExperiment
#' @param var        svar or fvar
#' @param separators character vector: possible separators to look for
#' @param verbose    TRUE or FALSE
#' @param ...        used for proper S3 method dispatch
#' @return separator (string) or NULL (if no separator could be identified)
#' @examples
#' # charactervector
#'    x <- c('PERM_NON.R1[H/L]', 'PERM_NON.R2[H/L]', 'PERM_NON.R3[H/L]')
#'    guess_sep(x)
#'
#'    x <- c('WT untreated 1', 'WT untreated 2', 'WT treated 1')
#'    guess_sep(x)
#'
#'    x <- c('group1', 'group2', 'group3.R1')
#'    guess_sep(x)
#'
#' # SummarizedExperiment
#'    # file <- download_data('halama18.metabolon.xlsx')
#'    # object <- read_metabolon(file, plot=FALSE)
#'    # guess_sep(object)
#'
#'    # file <- download_data('billing16.proteingroups.txt')
#'    # object <- read_proteingroups(file, plot=FALSE)
#'    # guess_sep(object)
#' @export
guess_sep <- function (x, ...)  UseMethod("guess_sep", x)


#' @rdname guess_sep
#' @export
guess_sep.character <- function(
    x, separators = c('.', '_'), verbose = FALSE, ...
){
# Initialize
    . <- NULL
    sep_freqs <-Map(function(y) stri_split_fixed(x, y), separators) %>%
                lapply(function(y) vapply(y, length, integer(1)))            %>%
                extract( vapply(., has_identical_values, logical(1)))        %>%
                vapply(unique, integer(1))
# No separator detected - return NULL
    if (all(sep_freqs==1)){
        if (verbose) message(x[1],': no (consistent) separator. Returning NULL')
        return('NOSEP')   # no separator detected
    }
# Find best separator
    best_sep <- sep_freqs %>%
                extract(.!=1)  %>%
                extract(is_max(vapply(., extract, integer(1), 1)))   %>%
                names()
# Ambiguous separator - take first from tail
    if (length(best_sep)>1){
        pattern <- best_sep %>% paste0(collapse='') %>% paste0('[', ., ']')
        best_sep <- x[1] %>% stri_extract_last_regex(pattern)
    }
# Separator identified - return
    if (verbose) message("\t\tGuess sep: '", best_sep, "'")
    return(best_sep)
}


#' @rdname guess_sep
#' @export
guess_sep.factor <- function(x, ...)  guess_sep.character(levels(x))


#' @rdname guess_sep
#' @export
guess_sep.SummarizedExperiment <- function(
    x, var = 'sample_id', separators =  c('.', '_'), verbose = FALSE, ...
){
    assert_is_subset(var, c(svars(x), fvars(x)))
    (if (var %in% svars(x)) slevels(x, var) else flevels(x, var)) %>%
    guess_sep(separators = separators, verbose = verbose)
}


#=============================================================================
#
#                       nfactors
#                       split_extract
#
#=============================================================================

#' @export
#' @rdname split_extract
nfactors <- function(x, sep = guess_sep(x)){
    length(unlist(stri_split_fixed(x[1], sep)))
}

#' stri_split and extract
#' @param x string
#' @param i integer
#' @param sep string
#' @return character
#' @examples
#' require(magrittr)
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' x <- object$sample_id[1:5]
#' nfactors(x)
#' split_extract(x, 1:2)
#' split_extract(x, seq_len(nfactors(x)-1))
#' split_extract(x, nfactors(x))
#' @export
split_extract <- function(x, i, sep=guess_sep(x)){
    factors <- stri_split_fixed(x, sep)
    vapply(factors, function(y) paste0(y[i], collapse=sep), character(1))
}


#=============================================================================
#
#               merge_sfile
#                   file_exists
#                       default_sfile
#
#=============================================================================

# Deals properly with NULL values
# file.exists does not!
file_exists <- function(file){
    if (is.null(file))      return(FALSE)
    if (file.exists(file))  return(TRUE)
                            return(FALSE)
}

#' Default sfile
#' @param file data file
#' @return sample file
#' @examples
#' file <- download_data('billing19.proteingroups.txt')
#' default_sfile(file)
#' @export
default_sfile <- function(file){
    sfile <- tools::file_path_sans_ext(file)
    sfile %<>% paste0('.samples.txt')
    sfile
}


#==============================================================================
#
#                       create_design
#                           single_subgroup
#                           are_factor
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

singlelevel <- function(object, svar) nlevels(object, svar) ==1
multilevel  <- function(object, svar) nlevels(object, svar) > 1

#' Create design
#'
#'  Create design matrix  for statistical analysis
#' @param object  SummarizedExperiment
#' @param formula formula with svars
#' @param verbose whether to message
#' @return design matrix
#' @examples
#' file <- download_data('billing16.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' create_design(object)
#'
#' object$subgroup <- 'billing16'
#' create_design(object)
#'
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' create_design(object)
#' create_design(object, ~ 0 + subgroup + Sex + T2D + age + bmi)
#' object$subgroup <- 'atkin18'
#' create_design(object)
#' @export
create_design <- function(object, formula = NULL, verbose = TRUE){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    if (is.null(formula))  formula <- 
        if (multilevel(object, 'subgroup')) ~ 0 + subgroup else ~ 1
    assert_is_subset(all.vars(formula), svars(object))
    . <- NULL
# Ensure that subgroup vector is a factor to preserve order of levels
    object$subgroup %<>% factor() # required for '(Intercept)' -> factor1.level1
    for (var in setdiff(all.vars(formula), 'subgroup')){
        if (is.character(sdata(object)[[var]])){
            sdata(object)[[var]] %<>% factor() }}
# Create design matrix
    if (verbose)  cmessage('\t\tDesign: %s', deparse(formula))
    myDesign <- model.matrix(formula,  data = sdata(object))
# Rename "(Intercept)" column
    if (ncol(myDesign)==1){ 
        if (singlelevel(object, 'subgroup'))  colnames(myDesign) <- 
                                                slevels(object, 'subgroup')
        return(myDesign) }
    factors <- svars(object)[are_factor(sdata(object))] # Rename intercept
    factor1 <- factors[1]                               # to factor1level1
    level1  <- levels(sdata(object)[[factor1]])[1]
    colnames(myDesign) %<>% gsub('(Intercept)', level1, ., fixed = TRUE)
# Rename regressors
    for (var in factors) colnames(myDesign) %<>% gsub(var, '', ., fixed = TRUE)
        # Fails for e.g. T2D = YES/NO: a meaningless column "YES" is created
        # For other cases it works wonderfully, so I keep it for now.
        # If it gives too many issues, roll back to doing the dropping only
        # for "subgroup" levels:
        #colnames(myDesign) %<>% gsub('subgroup', '', ., fixed=TRUE)
# Validify names
    colnames(myDesign) %<>% gsub(':', '..', ., fixed = TRUE)
    colnames(myDesign) %<>% make.names()
# Return
    return(myDesign)
}


#=============================================================================
#
#               subgroup_matrix
#                   split_subgroup_levels
#                       split_subgroup_values
#                           split_values
#
#=============================================================================

split_values <- function(x){
    sep <- guess_sep(x)
    dt <- data.table::data.table(x = x)
    dt[, data.table::tstrsplit(x, sep) ]
}

split_subgroup_values <- function(object){
    subgroupvalues <- subgroup_values(object)
    cbind(subgroup = subgroupvalues, split_values(subgroupvalues))
}

split_subgroup_levels <- function(object){
    subgrouplevels <- subgroup_levels(object)
    cbind(subgroup = subgrouplevels, split_values(subgrouplevels))
}


#' @rdname subgroup_matrix
#' @export
subgroup_array <- function(object){
    x <- subgroup_levels(object)
    sep <- guess_sep(object)
    #x %<>% sort()
    dt <- data.table(subgroup = factor(x, x))
    components <- dt[, tstrsplit(subgroup, sep, fixed=TRUE)]
    for (i in seq_len(ncol(components)))   components[[i]] %<>%
                                    factor(., levels=unique(.))
    dt %<>% cbind(components)
    data.table::setorderv(dt, rev(names(components)))
    levels  <- dt[, -1] %>% lapply(unique)
    #levels[1:2] %<>% rev()
    nlevels <- levels %>% vapply(length, integer(1))
    array(dt$subgroup, dim = nlevels, dimnames = levels)
}


#' Get subgroup matrix
#'
#' Arrange (subgroup)levels in matrix
#'
#' @param object SummarizedExperiment
#' @return matrix
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' subgroup_matrix(object)
#' @export
subgroup_matrix <- function(object){
    subgroup_array <- subgroup_array(object)
    if (length(dim(subgroup_array))==1)  return(matrix(subgroup_array,
        byrow=TRUE, nrow=1, dimnames=list(NULL, subgroup_array)))
    otherdims <- names(dim(subgroup_array)) %>% setdiff('V1')
    ncol1   <- Reduce('*', dim(subgroup_array)[otherdims])
    colnames1 <- dimnames(subgroup_array)[otherdims] %>%
                expand.grid()                        %>%
                apply(1, paste0, collapse='.')
    subgroupmat <- matrix(subgroup_array,
                        nrow = nrow(subgroup_array), ncol = ncol1,
                        dimnames=list(rownames(subgroup_array), colnames1))
    subgroupmat %>% extract(nrow(.):1, )
    #dt <- split_subgroup_levels(object)
    #subgroupmat <- as.matrix(data.table::dcast(
    #    dt, V1 ~ V2, value.var = 'subgroup'), rownames = 'V1')
    #subgroupmat %>% extract(rev(order(rownames(.))), order(colnames(.)))
}





#=============================================================================
#
#               contrast_subgroup_cols
#               contrast_subgroup_rows
#               col_contrast_defs
#               row_contrast_defs
#
#==============================================================================


#' Row/Col contrasts
#' @param object SummarizedExperiment
#' @return  matrix
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' subgroup_matrix(object)
#' contrast_subgroup_cols(object)
#' contrast_subgroup_rows(object)
#' @export
contrast_subgroup_cols <- function(object){
    subgroupmat <- subgroup_matrix(object)
    if (is_scalar(subgroupmat))  return(subgroupmat)
    if (ncol(subgroupmat)==1) return(matrix(, ncol=0, nrow=nrow(subgroupmat)))
    colcontrasts <- matrix(  sprintf('%s - %s',
                                    subgroupmat[, -1],
                                    subgroupmat[, -ncol(subgroupmat)]),
                            nrow = nrow(subgroupmat),
                            ncol = ncol(subgroupmat)-1)
    rownames(colcontrasts) <- rownames(subgroupmat)
    colnames(colcontrasts) <- sprintf('%s - %s',
                            colnames(subgroupmat)[-1],
                            colnames(subgroupmat)[-ncol(subgroupmat)])
    colcontrasts
}


#' @rdname contrast_subgroup_cols
#' @export
contrast_subgroup_rows <- function(object){
    subgroupmat <- subgroup_matrix(object)
    if (nrow(subgroupmat)==1) return(matrix(, nrow=0, ncol=ncol(subgroupmat)))
    rowcontrasts <- matrix(  sprintf('%s - %s',
                                    subgroupmat[-nrow(subgroupmat), ],
                                    subgroupmat[-1, ]),
                            nrow = nrow(subgroupmat)-1,
                            ncol = ncol(subgroupmat))
    colnames(rowcontrasts) <- colnames(subgroupmat)
    rownames(rowcontrasts) <- sprintf('%s - %s',
                            rownames(subgroupmat)[-nrow(subgroupmat)],
                            rownames(subgroupmat)[-1])
    rowcontrasts
}

contrast_subgroups <- function(object, design){
    if (ncol(design)==1){
        list(matrix(colnames(design), nrow=1, ncol=1), 
            matrix(nrow=0, ncol=0))
    } else {
        list(contrast_subgroup_cols(object),
            contrast_subgroup_rows(object)) }
}

#==============================================================================
#
#                            lm / lme / lmer
#
#==============================================================================

.cbindstats <- function(
    fitres, fname = 'F value', f.pname = 'Pr(>F)', pname = 'Pr(>|t|)', 
    tname = 't value', effectname = 'Estimate', sename = 'Std. Error'
){
    fval <- stats::anova(fitres)['subgroup', fname]
    f.p <- stats::anova(fitres)['subgroup', f.pname]
    fitres %<>% summary()
    fitres %<>% stats::coefficients()
    pvalues <- data.table::as.data.table(as.list(fitres[, pname]))
    tvalues <- data.table::as.data.table(as.list(fitres[, tname]))
    effects <- data.table::as.data.table(as.list(fitres[, effectname]))
    stderrs <- data.table::as.data.table(as.list(fitres[, sename]))
    names(pvalues) %<>% paste0('p.', .)
    names(tvalues) %<>% paste0('t.', .)
    names(effects) %<>% paste0('effect.', .)
    names(stderrs) %<>% paste0('se.', .)
    cbind(pvalues, tvalues, effects, stderrs, F = fval, F.p = f.p )
}

.lmtest <- function(sd, formula, block=NULL){
    fitres <- lm(  formula    = formula, 
                    data      = sd,
                    na.action = stats::na.omit)
    .cbindstats(fitres)
}

.lmetest <- function(sd, formula, block){
    fitres <- lme(  fixed     = formula, 
                    random    = block, 
                    data      = sd,
                    na.action = stats::na.omit)
    .cbindstats( fitres, fname='F-value', f.pname='p-value', pname='p-value', 
                tname = 't-value', effectname = 'Value', sename = 'Std.Error')
}

.lmertest <- function(sd, formula, block = NULL){
    fitres <- lme4::lmer(
        formula   = formula,
        data      = sd,
        na.action = stats::na.omit, 
        control   = lme4::lmerControl(check.conv.singular = lme4::.makeCC(
                                            action = "ignore",  tol = 1e-4)))
                    # https://stackoverflow.com/a/55367171
    fitres %<>% lmerTest::as_lmerModLmerTest()
    .cbindstats(fitres)
}

.extractstat <- function(lmeres, quantity){
    idx <- stri_startswith_fixed(names(lmeres), paste0(quantity, '.'))
    lmemat <- as.matrix(lmeres[, idx, with=FALSE])
    rownames(lmemat) <- lmeres$feature_id
    colnames(lmemat) %<>% stri_replace_first_fixed(paste0(quantity, '.'), '')
    lmemat
}

#' Statistical tests supported in autonomics
#' @export
TESTS <- c('limma','lm','lme','lmer', 'wilcoxon')


lmxtest <- function(
    object, 
    test, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = as.formula(sprintf(
                'value ~ %s', if (is.null(subgroupvar)) 1 else subgroupvar)),
    block = NULL, contrastdefs = NULL, verbose = TRUE, plot =  TRUE
){
# Assert
    assert_is_a_string(test);  assert_is_subset(test, TESTS)
# Prepare dt
    allx <- c(setdiff(all.vars(formula), 'value'), all.vars(block))
    dt <- sumexp_to_long_dt(object, svars = allx)
    fixedx <- setdiff(allx, all.vars(block))
    for (x in fixedx){  
        dt[[x]] %<>% factor()
        stats::contrasts(dt[[x]]) <- MASS::contr.sdif(levels(dt[[x]])) }
# lmefit
    fitmethod <- get(paste0('.', test, 'test'))
    lmeres <- dt[, fitmethod(.SD, formula=formula, block=block),by='feature_id']
    fdt <- lmeres[, c('feature_id', 'F', 'F.p'), with=FALSE]
    setnames(fdt, c('F', 'F.p'), c('F.lme', 'F.p.lme'))
    object %<>% merge_fdata(fdt)
# lmeextract
    for (x in setdiff(all.vars(formula), 'value'))  names(lmeres) %<>% 
        stri_replace_first_fixed(x, '')
    quantities <- c('p', 't', 'effect', 'se')
    lmeres <- mapply(.extractstat, quantity=quantities, 
                            MoreArgs = list(lmeres=lmeres), SIMPLIFY = FALSE)
    lmeres$fdr  <- apply(lmeres$p, 2, p.adjust, 'fdr')
    lmeres$bonf <- apply(lmeres$p, 2, p.adjust, 'bonf')
    metadata(object)$lmx <- do.call(abind::abind, c(lmeres, along=3))
    metadata(object)$lmx %<>% extract(rownames(object),,)
    names(metadata(object)) %<>% stri_replace_first_fixed('lmx', test)
    
    if (verbose) cmessage_df('\t\t\t%s', extract_test_summary(object, test))
    if (is.null(contrastdefs)) contrastdefs <-colnames(metadata(object)[[test]])
    if (length(contrastdefs) > 1) contrastdefs %<>% setdiff('(Intercept)')
    if (plot)  print(plot_volcano(
        object, test=test, contrastdefs = contrastdefs)) 
    object
}


#' @rdname limmatest
#' @export
lmtest <- function(
    object,
    test, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = as.formula(sprintf(
                'value ~ %s', if (is.null(subgroupvar)) 1 else subgroupvar)),
    block = NULL, contrastdefs = NULL, verbose = TRUE, plot =  TRUE
){
    if (verbose)  cmessage('\t\tlm(%s)', Reduce(paste, deparse(formula)))
    lmxtest(object, test = 'lm', subgroupvar = subgroupvar, 
            formula = formula, block = block, contrastdefs = contrastdefs, 
            verbose = verbose, plot = plot)
}


#' @rdname limmatest
#' @export
lmetest <- function(
    object, 
    test, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = as.formula(sprintf(
                'value ~ %s', if (is.null(subgroupvar)) 1 else subgroupvar)),
    block = NULL, contrastdefs = NULL, verbose = TRUE, plot =  TRUE
){
    assert_is_not_null(block)
    if (is_a_string(block)){ 
        if (is_subset(block, svars(object)))  block %<>% sprintf('~1|%s', .)
        block %<>% as.formula() }
    if (verbose)  cmessage('\t\tlme(%s, random = %s)', 
                        Reduce(paste, deparse(formula)), 
                        Reduce(paste, deparse(block)))
    lmxtest(object, test = 'lme', subgroupvar = subgroupvar, 
            formula = formula, block = block, contrastdefs = contrastdefs, 
            verbose = verbose, plot = plot)
}


#' @rdname limmatest
#' @export
lmertest <- function(
    object, 
    test, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = as.formula(sprintf(
                'value ~ %s', if (is.null(subgroupvar)) 1 else subgroupvar)),
    block = NULL, contrastdefs = NULL, verbose = TRUE, plot =  TRUE
){
    if (is_formula(block))  block <- Reduce(paste, deparse(block))
    if (is_a_string(block)){ 
        if (is_subset(block, svars(object)))  block %<>% sprintf('1|%s', .)
        formula <- Reduce(paste, deparse(formula))
        formula %<>% sprintf('%s + (%s)', ., block)
        formula %<>% as.formula()
    }
    if (verbose)  cmessage('\t\tlmer(%s)', Reduce(paste, deparse(formula)))
    lmxtest(object, test = 'lmer', subgroupvar = subgroupvar, 
            formula = formula, block = NULL, contrastdefs = contrastdefs, 
            verbose = verbose, plot = plot)
}


#==============================================================================
#
#                                limma
#
#==============================================================================

contrvec2mat  <- function(contrastdefs)  matrix(
                    contrastdefs, nrow=1, dimnames=list("", contrastdefs))

contrmat2list <- function(contrastdefs)  list(colcontrasts = contrastdefs)

vectorize_contrastdefs <- function(contrastdefs){
    unname(unlist(lapply(contrastdefs, function(x) na.exclude(c(t(x))))))
}

.limmacontrast <- function(object, fit){
    # compute contrasts
    contrastmat <- makeContrasts(
        contrasts = vectorize_contrastdefs(contrastdefs(object)),
        levels    = design(object))
    fit %<>% contrasts.fit(contrasts = contrastmat)
    limma_quantities <- if (all(fit$df.residual==0)){ c('effect', 'rank')
    } else { c('effect','rank','t','se','p','fdr','bonf') }
    limma(object) <- array( dim=c(nrow(fit),ncol(fit),length(limma_quantities)),
                            dimnames = list(feature  = rownames(fit),
                                            contrast = colnames(fit),
                                            quantity = limma_quantities))
    limma(object)[,,'effect'] <- fit$coefficients
    limma(object)[,,'rank'  ] <- apply(-abs(fit$coefficients), 2, rank)
    # perform moderated t test
    if (!all(fit$df.residual==0)){
        fit %<>% eBayes()
        pp <- fit$p.value
        limma(object)[,,'t' ] <- fit$t
        limma(object)[,,'se'] <- sqrt(fit$s2.post) * fit$stdev.unscaled
        limma(object)[,,'p' ] <- pp
        limma(object)[,,'rank'] <- apply(pp, 2, rank)
        limma(object)[,,'fdr' ] <- apply(pp, 2, p.adjust, 'fdr')
        limma(object)[,,'bonf'] <- apply(pp, 2, p.adjust, 'bonferroni')
        fdata(object)$F.limma   <- fit$F
        fdata(object)$F.p.limma <- fit$F.p
    }
    return(object)
}

#' Test for differential expression
#'
#' Limma results can be easily accessed with limma(object).
#' @param object       SummarizedExperiment
#' @param contrastdefs contrastdef vector / matrix / list
#' \itemize{
#' \item{c("t1-t0", "t2-t1", "t3-t2")}
#' \item{matrix(c("WT.t1-WT.t0", "WT.t2-WT.t1", "WT.t3-WT.t2"), \cr
#'      c("KD.t1-KD.t0", "KD.t2-KD.t1", "KD.t3-KD.t2"), nrow=2, byrow=TRUE)}
#' \item{list(matrix(c("WT.t1-WT.t0", "WT.t2-WT.t1", "WT.t3-WT.t2"), \cr
#'      c("KD.t1-KD.t0", "KD.t2-KD.t1", "KD.t3-KD.t2"), nrow=2, byrow=TRUE), \cr
#'      matrix(c("KD.t0-WT.t0", "KD.t1-WT.t1", "KD.t2-WT.t2", "KD.t3-WT.t3"),\cr
#'      nrow=1, byrow=TRUE))}}
#' @param formula   designmat formula
#' @param block     block svar (or NULL)
#' @param verbose   whether to msg
#' @param plot      whether to plot
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' lmtest(object)
#' wilcoxtest(object)
#' limmatest(object, block = 'Subject_ID')
#' lmetest(  object, block = 'Subject_ID')
#' lmertest( object, block = 'Subject_ID')
#'
#' file <- download_data('billing19.proteingroups.txt')
#' select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
#' select %<>% paste0('_STD')
#' object <- read_proteingroups(file, select_subgroups = select, plot = FALSE)
#' object %<>% impute_systematic_nondetects(plot=FALSE)
#' object %<>% limmatest()
#' object %<>% lmtest()
#'
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE, voom=TRUE)
#' object %<>% limmatest(plot=FALSE)
#' weights(object) <- NULL; object %<>% limmatest(plot=FALSE)
#' object %<>% lmtest(plot=FALSE)
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% impute_systematic_nondetects(plot=FALSE)
#' object %<>% limmatest(plot=FALSE)
#' object %<>% lmtest(plot=FALSE)
#' @export
limmatest <- function(object, 
    subgroupvar = if ('subgroup' %in% svars(object)) 'subgroup' else NULL, 
    formula = as.formula(sprintf(
                '~ 0 + %s', if (is.null(subgroupvar)) 1 else subgroupvar)),
    contrastdefs = NULL, block = NULL, verbose = TRUE, plot =  TRUE
){
# Initialize
    assert_is_all_of(object, 'SummarizedExperiment')
    design <- create_design(object, formula=formula, verbose = verbose)
    design(object)    <- design
    if (is.null(contrastdefs)) contrastdefs <- contrast_subgroups(object,design)
    if (is.character(contrastdefs)) contrastdefs %<>% contrvec2mat()
    if (is.matrix(contrastdefs))    contrastdefs %<>% contrmat2list()
    contrastdefs(object) <- contrastdefs
# Block
    if (verbose)  cmessage('\t\tLmfit: %s%s',
        if(is.null(block))            '' else paste0('block on `', block, '`'),
        if(is.null(weights(object)))  '' else paste0(' (weights)'))
    if (!is.null(block)){
        assert_is_subset(block, svars(object))
        blockvar <- block
        block <- sdata(object)[[block]]
        if (is.null(metadata(object)$blockcor)){
            if (verbose)  cmessage('\t\t\t\tdupcor `%s`', blockvar)
            metadata(object)$blockcor <- duplicateCorrelation(
                exprs(object), design=design, block=block
            )$consensus.correlation }}
# Lmfit
    fit <- suppressWarnings(lmFit(object = exprs(object[, rownames(design)]),
        design = design, block = block, correlation = metadata(object)$blockcor,
        weights = weights(object)))
# Contrast
    if (verbose)  cmessage('\t\tContrast:')
    object %<>% .limmacontrast(fit)
    if (plot)  print(plot_volcano(object, test='limma')) 
                    # plot_contrastogram(object)
    if (verbose) cmessage_df('\t\t\t%s', extract_test_summary(object, 'limma'))
    return(object)
}


#==============================================================================
#
#                                wilcoxon
#
#==============================================================================


wilcoxtest <- function(object, formula = ~ subgroup, 
    subgroupvar = all.vars(formula), 
    contrastdefs = contrast_subgroups(object, create_design(object, formula)), 
    block = NULL, verbose = TRUE, plot = TRUE
){
# fit
    dt <- sumexp_to_long_dt(object, svars = c(subgroupvar, block))
    results <- lapply(vectorize_contrastdefs(contrastdefs), 
                    .wilcoxtest, dt, subgroupvar, block)
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
    if (plot)  print(plot_volcano(object, test='wilcoxon')) 
                    # plot_contrastogram(object)
    if (verbose) cmessage_df('\t\t\t%s', extract_test_summary(object, 'wilcoxon'))
    object
}


.wilcoxtest <- function(contrastdef, dt, subgroupvar, block){
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
    

#==============================================================================
#
#                    limma & limma<-
#                    extract_test_dt
#
#==============================================================================

#' @title Get/set design matrix
#' @param object SummarizedExperiment
#' @param value list
#' @return design (get) or SummarizedExperiment (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'               file, invert_subgroups=inv, test='limma', plot=FALSE)
#' design(object)
#' @export
setGeneric("design", function(object)   standardGeneric("design") )


#' @rdname design
#' @export
setGeneric("design<-",
function(object, value)  standardGeneric("design<-") )


#' @rdname design
setMethod("design", signature("SummarizedExperiment"),
function(object) metadata(object)$design)


#' @rdname design
setReplaceMethod("design", signature("SummarizedExperiment", "matrix"),
function(object, value){
    metadata(object)$design <- value
    object  })


#' @title Get/set contrastdefs
#' @param object SummarizedExperiment
#' @param value list
#' @return contrastdefs (get) or SummarizedExperiment (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'             file, invert_subgroups=inv, test='limma', plot=FALSE)
#' contrastdefs(object)
#' @export
setGeneric("contrastdefs", function(object)   standardGeneric("contrastdefs") )


#' @rdname contrastdefs
setMethod("contrastdefs", signature("SummarizedExperiment"),
function(object) metadata(object)$contrastdefs)


#' @rdname contrastdefs
#' @export
setGeneric("contrastdefs<-",
function(object, value)  standardGeneric("contrastdefs<-") )


#' @rdname contrastdefs
setReplaceMethod("contrastdefs", signature("SummarizedExperiment", "list"),
function(object, value){
    metadata(object)$contrastdefs <- value
    object  })


#' @title Get/set limma results
#' @description Get/Set limma results
#' @param object SummarizedExperiment
#' @param value list
#' @return limma results (get) or updated object (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'            file, invert_subgroups=inv, test='limma', plot=FALSE)
#' dim(limma(object))
#' dim(limma(object[1:5, ]))
#' @export
setGeneric("limma", function(object)   standardGeneric("limma") )


#' @rdname limma
setMethod("limma", signature("SummarizedExperiment"),
function(object){
    limma_array <- metadata(object)$limma
    if (is.null(limma_array)) NULL else limma_array[
                                            fnames(object), , , drop=FALSE] })

#' @rdname limma
#' @export
setGeneric("limma<-", function(object, value)  standardGeneric("limma<-") )


#' @rdname limma
setReplaceMethod("limma", signature("SummarizedExperiment", "array"),
function(object, value){
    metadata(object)$limma <- value
    object  })


#' @rdname limma
setReplaceMethod("limma", signature("SummarizedExperiment", "NULL"),
function(object, value) object)



#' Extract limma quantity
#' @param object SummarizedExperiment
#' @param quantity 'effect', 'p', 'fdr', 'bonf'
#' @return melted data.table
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'             file, invert_subgroups=inv, test='limma', plot=FALSE)
#' extract_test_quantity(object, test = 'limma')
#' @noRd
extract_test_quantity <- function(object, test, quantity='p'){
    assert_is_subset(test, TESTS)
    fvars0 <- c('feature_id','feature_name','imputed')
    fvars0 %<>% intersect(fvars(object))
    fdt <- data.table(fdata(object)[, fvars0, drop=FALSE])
    testdt <- metadata(object)[[test]][, , quantity, drop=FALSE]
    testdt %<>% adrop(drop=3)
    fdt %<>% cbind(testdt)
    data.table::melt.data.table(
        fdt, id.vars = fvars0, variable.name = 'contrast', value.name=quantity)
}


merge_test_quantities <- function(x, y){
    names0 <- c('feature_id','feature_name','imputed', 'contrast')
    names0 %<>% intersect(names(x)) %>% intersect(names(y))
    merge(x, y, by = names0)
}


#' Extract limma datatable
#' @param object SummarizedExperiment
#' @param test 'limma', 'lme', 'lm', 'wilcoxon'
#' @return melted data.table
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, test='limma', plot=FALSE)
#' extract_test_dt(object, test = 'limma')
#'
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'            file, invert_subgroups=inv, test='limma', plot=FALSE)
#' extract_test_dt(object, test = 'limma')
#' @noRd
extract_test_dt <- function(object, test){
    Reduce( merge_test_quantities,
            mapply( extract_test_quantity,
                    quantity = c('effect', 'p', 'fdr', 'bonf'),
                    MoreArgs = list(object=object, test=test),
                    SIMPLIFY = FALSE ))
}

#' Extract contrast analysis summary
#' @param object SummarizedExperiment
#' @param test 'limma', 'lme', 'lm', 'wilcoxon'
#' @return data.table(contrast, nup, ndown)
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, test='limma', plot=FALSE)
#' extract_test_summary(object, 'limma')
#' @export
extract_test_summary <- function(object, test){
    effect <- fdr <- NULL
    extract_test_dt(object, test = test)[,
        .(  ndown = sum(effect<0 & fdr<0.05, na.rm=TRUE),
            nup   = sum(effect>0 & fdr<0.05, na.rm=TRUE)),
        by='contrast']
}

#=============================================================================
#
#             plot_contrastogram
#                 compute_connections
#
#=============================================================================


true_names <- function(x) names(x[x])

compute_connections <- function(
    object, colors = make_colors(subgroup_levels(object), guess_sep(object))
){
# subgroup matrix, difference contrasts, limma
    pvalues <- limma(object)[, , 'p',      drop=FALSE]
    effects <- limma(object)[, , 'effect', drop=FALSE]
    nsignif <- apply(pvalues < 0.05, 2, sum, na.rm=TRUE)
                #colSums( pvalues < 0.05, na.rm=TRUE)  # BREAKS ON SINGLE CONTR!
    nup     <- apply(pvalues < 0.05 & effects>0, 2, sum, na.rm=TRUE)
                #colSums((pvalues < 0.05) & (effects > 0), na.rm=TRUE)
    ndown   <- apply(pvalues < 0.05 & effects<0, 2, sum, na.rm=TRUE)
                #colSums((pvalues < 0.05) & (effects < 0), na.rm=TRUE)
# Create diagram
    sep <- guess_sep(object)
    subgroupmatrix <- subgroup_matrix(object)
    subgrouplevels <- c(t(subgroupmatrix))
    arrowsizes <- arrowcolors <- matrix(0,
        nrow = length(subgrouplevels), ncol = length(subgrouplevels),
        dimnames = list(subgrouplevels, subgrouplevels))
    arrowlabels <- matrix("0", nrow = nrow(arrowsizes), ncol = ncol(arrowsizes),
                        dimnames = dimnames(arrowsizes))
# Add contrast numbers
    design    <- design(object)
    colcontrasts <- contrastdefs(object)[[1]]
    rowcontrasts <- contrastdefs(object)[[2]]
    contrastdefs <-  c( c(t(colcontrasts)), c(t(rowcontrasts)))
    contrastmat  <- makeContrasts(contrasts = contrastdefs, levels = design)
    for (contrastname in colnames(contrastmat)){
        contrastvector <- contrastmat[, contrastname]
        to   <- true_names(contrastvector>0)
        from <- if (any(contrastvector<0)) true_names(contrastvector<0) else to
        ns <- nsignif[[contrastname]]
        nu <- nup[[contrastname]]
        nd <- ndown[[contrastname]]
        arrowsizes[ to, from] <- nu#ns
        arrowsizes[ from, to] <- nd#ns
        arrowcolors[to, from] <- colors[[to]]
        arrowcolors[from, to] <- colors[[from]]
        arrowlabels[to, from] <- if (nu>0) nu else 0
                            #paste0(nu,  " %up% phantom(.)") else "phantom(.)"
        arrowlabels[from, to] <- if (nd>0) nd else 0
                            #paste0(nd," %down% phantom(.)") else "phantom(.)"
    }
# Return
    #arrowlabels[arrowcolors==0] <- "0"
    list(arrowsizes = arrowsizes,
        arrowcolors = arrowcolors,
        arrowlabels = arrowlabels)
}


#' Plot contrastogram
#' @param object SummarizedExperiment
#' @param colors named color vector (names = subgroups)
#' @param curve  arrow curvature
#' @return list returned by \code{\link[diagram]{plotmat}}
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, test='limma', plot=FALSE)
#' plot_contrastogram(object)
#' @export
plot_contrastogram <- function(
    object,
    colors = make_colors(subgroup_levels(object), guess_sep(object)),
    curve  = 0.1
){
# Initialize
    V2 <- N <- NULL
# Prepare
    contrastogram_matrices <- compute_connections(object, colors = colors)
    arrowsizes  <- contrastogram_matrices$arrowsizes
    arrowcolors <- contrastogram_matrices$arrowcolors
    arrowlabels <- contrastogram_matrices$arrowlabels
    widths <- scales::rescale(arrowsizes, c(0.01,30))
# Plot
    dt <- split_subgroup_levels(object)
    nrow <- dt[, data.table::uniqueN(V2)]
    nperrow <- dt[, .N, by = 'V1'][, N]
    if (all(nperrow==1)) nperrow %<>% length()
    #dir.create('~/autonomicscache/contrastogram')
    #pdf('~/autonomicscache/contrastogram/directed_contrastogram.pdf',
    #width = 9, height = 9)
    arrowlabels %<>% as.data.frame()
    diagram::plotmat(A          = arrowlabels,
                    pos         = nperrow,
                    curve       = curve,
                    name        = rownames(arrowsizes),
                    relsize     = 1,
                    box.size    = 0.05,
                    box.col     = colors[rownames(arrowsizes)],
                    box.type    = 'square',
                    box.prop    = 0.8,
                    arr.lwd     = widths,
                    shadow.size = 0, # sqrt(arrowsizes)
                    arr.lcol    = arrowcolors,
                    arr.col     = arrowcolors,
                    arr.type    = 'triangle')
    #, arr.lcol = log2(1+diagram_matrix))
    #dev.off()
}



#==============================================================================
#
#          plot_volcano
#              make_volcano_dt
#                  top_down/top_up
#                      nmax/nmin
#
#==============================================================================

invwhich <- function(indices, totlength) is.element(seq_len(totlength), indices)

#' Return nth max (min) value in vector
#'
#' Orders a vector and returns n'th ordered value.
#' When vector length is smaller than n, returns last value.
#'
#' @param x numeric vector
#' @param n integer
#' @return value
#' @examples
#' nmax(c(1,2,3,4,5), 2)
#' nmin(c(1,2,3,4,5), 2)
#' @noRd
nmax <- function(x, n) sort(x, decreasing = TRUE) %>% extract(min(length(.), n))
nmin <- function(x, n) sort(x) %>% extract(min(length(.), n))

top_down <- function(effect, fdr, mlp, ntop){
    fdr_ok   <- fdr  < 0.05
    coef_ok  <- effect < -1
    coef_top <- if (any(fdr_ok)){  effect < nmin(effect[fdr_ok], ntop+1)
                } else {           rep(FALSE, length(effect))            }
    mlp_top  <- if (any(coef_ok)){ mlp  > nmax(mlp[coef_ok], ntop+1)
                } else {           rep(FALSE, length(effect))            }
    fdr_ok & coef_ok & (coef_top | mlp_top)
}

#' @examples
#' file <- download_data("billing16.proteingroups.txt")
#' invert_subgroups <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'           file, invert_subgroups=invert_subgroups, test='limma', plot=FALSE)
#' effect <-      limma(object)[,1,'effect']
#' fdr    <-      limma(object)[,1,'fdr']
#' mlp    <- -log(limma(object)[,1,'p'])
#' ntop   <- 3
#' table(top_up(effect, fdr, mlp, ntop))
#' table(top_down(effect, fdr, mlp, ntop))
#' @noRd
top_up <- function(effect, fdr, mlp, ntop){
    fdr_ok   <- fdr  < 0.05
    coef_ok  <- effect >  0 # currently no filter 
    coef_top <- if(any(fdr_ok)){  effect > nmax(effect[fdr_ok], ntop+1)
                } else {          rep(FALSE, length(effect)) }
    mlp_top  <- if (any(coef_ok)){ mlp > nmax(mlp[coef_ok], ntop+1)
                } else {           rep(FALSE, length(effect)) }
    fdr_ok & coef_ok & (coef_top | mlp_top)
}


melt_contrastdefs <- function(contrastdefmat){
    facetrow <- NULL
    contrastdefdt <- data.table(contrastdefmat, facetrow = "")
    if (!is.null(rownames(contrastdefmat))) contrastdefdt[,
                                        facetrow := rownames(contrastdefmat)]
    data.table::melt.data.table(
        contrastdefdt,
        id.vars       = 'facetrow',
        variable.name = 'facetcol',
        value.name    = 'contrast', by = 'contrast')
}


#' Create volcano datatable
#' @param object          SummarizedExperiment
#' @param test            'limma', 'lme', 'lm', 'wilcoxon'
#' @param contrastdefmat  contrastdef matrix
#' @param ntop            no of top features to be annotated
#' @return data.table
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, test='limma', plot=FALSE)
#' make_volcano_dt(object, test = 'limma')
#' @export
make_volcano_dt <- function(
    object, test, contrastdefmat = contrastdefs(object)[[1]], ntop = 3
){
    effect <- p <- mlp <- topdown <- topup <- significance <- fdr <- NULL
    dt <- extract_test_dt(object, test)
    dt %<>% merge(melt_contrastdefs(contrastdefmat), by = 'contrast')
    dt %<>% extract(!is.na(effect) & !is.na(p))
    dt[, mlp  := -log10(p)]

    # Prepare volcano datatable
    # Note: Using effect <= 0 (rather than effect <0) is required.
    # Otherwise the (very few) features with effect=0 will have no effect for
    # 'significance'
    by <- intersect(c('contrast', 'imputed'), names(dt))
    dt[,topdown := top_down(effect, fdr, mlp, ntop), by=by]
    dt[,topup   := top_up(  effect, fdr, mlp, ntop), by=by]
    dt[effect<=0,            significance := 'down']
    dt[effect> 0,            significance :=   'up']
    dt[effect<=0 & fdr<0.05, significance := 'down: fdr<0.05']
    dt[effect> 0 & fdr<0.05, significance :=   'up: fdr<0.05']
    dt[topdown==TRUE,        significance := 'down: top']
    dt[topup  ==TRUE,        significance :=   'up: top']
    dt[,significance := factor(significance, c(
        'down: top','down: fdr<0.05','down','up','up: fdr<0.05','up: top'))]
    dt[]
}

#' Plot volcano
#' @param object         SummarizedExperiment
#' @param test           'limma', 'lme', 'lm', 'wilcoxon'
#' @param contrastdefs   contrastdef vector / matrix / list
#' @param label          fvar for labeling top features
#' @param ntop           number: n top features to be annotated
#' @return ggplot object
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, test='limma', plot=FALSE)
#' plot_volcano(object)
#' @export
plot_volcano <- function(object, 
    test = intersect(names(metadata(object)), TESTS)[1], 
    contrastdefs = autonomics::contrastdefs(object)[[1]], 
    label = feature_name, ntop = 1
){
# Assert/Process
    assert_is_all_of(object, "SummarizedExperiment")
    assert_is_a_string(test)
    assert_is_subset(test, TESTS)
    assert_is_subset(test, names(metadata(object)))
    if (is.null(contrastdefs)) contrastdefs<-colnames(metadata(object)[[test]])
    if (is.character(contrastdefs)) contrastdefs %<>% contrvec2mat()
    if (is.list(contrastdefs))      contrastdefs %<>% extract2(1)
    assert_is_matrix(contrastdefs)
    topup <- topdown <- effect <- mlp <- facetrow <- facetcol <- NULL
    label <- enquo(label)
# Prepare
    plotdt <- make_volcano_dt(object, test, contrastdefs, ntop = ntop)
    txtdt  <- copy(plotdt)[topup==TRUE | topdown==TRUE]
    colorvalues <-c(hcl(h=  0, l=c(20, 70, 100), c=100), # 20 70 100
                    hcl(h=120, l=c(100, 70, 20), c=100)) # 100 70 20
    names(colorvalues) <- levels(plotdt$significance)
# Plot
    imputed <- NULL # fallback when plotdt misses "imputed"
    significance <- NULL
    p <- ggplot(plotdt) + facet_grid(
        rows = if (all(stri_isempty(plotdt$facetrow))) NULL else vars(facetrow),
        cols = if (all(stri_isempty(plotdt$facetcol))) NULL else vars(facetcol),
        scales = 'fixed') +
    geom_point(aes(x=effect,y=mlp,color=significance,shape=imputed),na.rm=TRUE)
    if (!quo_is_null(label)) p <- p + geom_text_repel(
                        data = txtdt,
                        aes(x=effect, y=mlp, label=!!label, color=significance),
                        #hjust = 'outward',
                        na.rm = TRUE,
                        show.legend = FALSE)#,
                        #direction = 'x'
    p + theme_bw() +
        scale_color_manual(values = colorvalues, name = NULL) +
        xlab('log2(FC)') +
        ylab('-log10(p)') +
        ggtitle(test)#+
        #guides(color=FALSE)
}
