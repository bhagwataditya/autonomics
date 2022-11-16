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
guess_sep.numeric <- function(x) NULL


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
#' @rdname split_extract_fixed
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
#' split_extract_fixed(x, '_', 1:2)
#' split_extract_fixed(x, '_', seq_len(nfactors(x)-1))
#' split_extract_fixed(x, '_', nfactors(x))
#' 
#' # With NA values
#' split_extract_fixed(fdt(object)$PUBCHEM, ';', 1)
#' @export
split_extract_fixed <- function(x, sep, i){
    if (is.factor(x))  x %<>% as.character()
    assert_is_a_string(sep)
    assert_is_numeric(i)
    factors <- stri_split_fixed(x, sep)
    extracted <- rep(NA_character_, length(factors))
    idx <- !is.na(factors)
    extracted[idx] <- vapply(factors[idx], 
                     function(y) paste0(y[i], collapse=sep), character(1))
    extracted
}

#' @export
#' @rdname split_extract_fixed
split_extract_regex <- function(x, sep, i){
    if (is.factor(x))  x %<>% as.character()
    assert_is_a_string(sep)
    assert_is_numeric(i)
    factors <- stri_split_regex(x, sep)
    extracted <- rep(NA_character_, length(factors))
    idx <- !is.na(factors)
    extracted[idx] <- vapply(factors[idx], 
                             function(y) paste0(y[i], collapse=sep), character(1))
    extracted
}

#' @export
#' @rdname split_extract_fixed
split_extract <- function(x, i, sep=guess_sep(x)){
    .Deprecated('split_extract_fixed')
    factors <- stri_split_fixed(x, sep)
    extracted <- rep(NA_character_, length(factors))
    idx <- !is.na(factors)
    extracted[idx] <- vapply(
                        factors[idx], 
                        function(y) paste0(y[i], collapse=sep), character(1))
    extracted
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

split_subgroup_values <- function(object, subgroupvar){
    subgroupvalues <- svalues(object, subgroupvar)
    cbind(subgroup = subgroupvalues, split_values(subgroupvalues))
}

split_subgroup_levels <- function(object, subgroupvar){
    subgrouplevels <- slevels(object, subgroupvar)
    cbind(subgroup = subgrouplevels, split_values(subgrouplevels))
}


#' @rdname subgroup_matrix
#' @export
subgroup_array <- function(object, subgroupvar){
    . <- NULL
    x <- slevels(object, subgroupvar)
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
#' @param subgroupvar subgroup svar
#' @return matrix
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' subgroup_matrix(object, 'Group')
#' @export
subgroup_matrix <- function(object, subgroupvar){
    . <- NULL
    subgroup_array <- subgroup_array(object, subgroupvar)
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

#==============================================================================
#
#                   pvars /  fdrvars  / effectvars
#                 pvalues / fdrvalues / effectvalues
#                             fits
#
#==============================================================================

#' Get p/fdr/effect/down/up/sign matrix
#' @param object    SummarizedExperiment
#' @param fit       string
#' @param coef      string
#' @return matrix (feature x coef)
#' @examples 
#' # Read
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     object <- read_metabolon(file, impute = TRUE, plot = FALSE)
#'     object %<>% fit_limma(subgroupvar = 'SET')
#'     object %<>% fit_lm(   subgroupvar = 'SET')
#' # Values
#'     effect(object)[1:3, ]
#'          p(object)[1:3, ]
#'        fdr(object)[1:3, ]
#'       down(object)[1:3, ]
#'         up(object)[1:3, ]
#'       sign(object)[1:3, ]
#' # Variables
#'     effectvar(object)
#'          pvar(object)
#'        fdrvar(object)
#' @export
p <- function(object, coef = coefs(object), fit = fits(object)){
    var <- pvar(object, coef = coef, fit = fit)
    dt <- fdt(object)[, var, with = FALSE]
    names(dt) %<>% stri_replace_first_fixed(paste0('p', FITSEP), '')
    mat <- as.matrix(dt)
    rownames(mat) <- rownames(object)
    mat
}


#' @rdname p
#' @export
pvar <- function(object, coef = coefs(object), fit = fits(object)){
    x <- expand.grid(var = 'p', fit = fit, coef = coef)
    x <-  paste(x$var, x$coef, x$fit, sep = FITSEP)
    x %<>% intersect(fvars(object))
    x
}


#' @rdname p
#' @export
effect <- function(object, coef = coefs(object), fit = fits(object)){
    var <- effectvar(object, coef = coef, fit = fit)
    dt <- fdt(object)[, var, with = FALSE]
    names(dt) %<>% stri_replace_first_fixed(paste0('effect', FITSEP), '')
    mat <- as.matrix(dt)
    rownames(mat) <- rownames(object)
    mat
}


#' @rdname p
#' @export
effectvar <- function(object, coef = coefs(object), fit = fits(object)){
    x <- expand.grid(var = 'effect', fit = fit, coef = coef)
    x <- paste(x$var, x$coef, x$fit, sep = FITSEP)
    x %<>% intersect(fvars(object))
    x
}


#' @rdname p
#' @export
fdr <- function(object, coef = coefs(object), fit = fits(object)){
    var <- fdrvar(object, coef = coef, fit = fit)
    dt <- fdt(object)[, var, with = FALSE]
    names(dt) %<>% stri_replace_first_fixed(paste0('fdr', FITSEP), '')
    mat <- as.matrix(dt)
    rownames(mat) <- rownames(object)
    mat
}

#' @rdname p
#' @export
fdrvar <- function(object, coef = coefs(object), fit = fits(object)){
    x <- expand.grid(var = 'fdr', fit = fit, coef = coef)
    x <- paste(x$var, x$coef, x$fit, sep = FITSEP)
    x %<>% intersect(fvars(object))  # all fits not necessarily all coefs, e.g:
    x               # limma(contrasts = .) generally run without intercept
}                   # lm(coefs = .)        typically run with intercept


bonvar <- function(object, coef = coefs(object), fit = fits(object)){
    x <- expand.grid(var = 'bonferroni', fit = fit, coef = coef)
    x <- paste(x$var, x$coef, x$fit, sep = FITSEP)
    x %<>% intersect(fvars(object))  # all fits not necessarily all coefs, e.g:
    x               # limma(contrasts = .) generally run without intercept
}                   # lm(coefs = .)        typically run with intercept

#' @rdname p
#' @export
up <- function(
    object, 
    coef   = coefs(object),
    fit    = fits(object),
    var    = 'fdr', 
    cutoff = 0.05
){
    fdrmat  <- get(var)(object, coef = coef, fit = fit)
    effectmat <- effect(object, coef = coef, fit = fit)
    y <- fdrmat < cutoff  &  effectmat > 0
    y[is.na(y)] <- FALSE
    mode(y) <- 'numeric'
    y
}

#' @rdname p
#' @export
down <- function(
    object, 
    coef    = coefs(object), 
    fit     = fits(object), 
    var     = 'fdr', 
    cutoff  = 0.05
){
    fdrmat  <- get(var)(object, coef = coef, fit = fit)
    effectmat <- effect(object, coef = coef, fit = fit)
    y <- fdrmat < cutoff  &  effectmat < 0
    y[is.na(y)] <- FALSE
    mode(y) <- 'numeric'
    y
}


#' Get sign
#' @export
sign.SummarizedExperiment <- function(
    object, 
    coef   = coefs(object),
    fit    = fits(object),
    var    = 'fdr',
    cutoff = 0.05
){ 
    if (length(coef)==0) return(matrix(0, nrow = nrow(object), ncol = ncol(object), dimnames = dimnames(object)))
    up(   object, coef = coef, fit = fit, var = var, cutoff = cutoff)  -
    down( object, coef = coef, fit = fit, var = var, cutoff = cutoff)
}

# dont rm - its the lower-level function used by fits() and coefs() !
.pvars <- function(object){
    . <- NULL
    fvars(object) %>% extract(stri_startswith_fixed(., paste0('p', FITSEP)))
}

#' Get fit models
#' 
#' @param object SummarizedExperiment
#' @return  character vector
#' @examples 
#' require(magrittr)
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, fit='limma', plot=FALSE)
#' fits(object)
#' @export
fits <- function(object){
    x <- .pvars(object)
    x %<>% split_extract_fixed(FITSEP, 3)
    x %<>% unique()
    x
}

#' Get fit coefs
#' 
#' @param object  SummarizedExperiment
#' @param fit    string: 'limma', 'lm', 'lme', 'lmer'
#' @param svars   NULL/charactervec : retain only coefficients relevant for this svar
#' @return  character vector
#' @examples 
#' require(magrittr)
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, fit='limma', plot=FALSE)
#' coefs(object)
#' @export
coefs <- function(object, fit = fits(object), svars = NULL){
    . <- NULL
    coefs0 <- split_extract_fixed(.pvars(object), FITSEP, 2)
    fits0  <- split_extract_fixed(.pvars(object), FITSEP, 3)
    coefs0 %<>% extract(fits0 %in% fit)
    coefs0 %<>% unique()
    if (!is.null(svars))  coefs0 %<>% extract(Reduce('|', lapply(svars, grepl, .)))
    coefs0 
}

#==============================================================================
#
#                    summarize_fit
#                    old_summarize_fit
#
#==============================================================================

#' Summarize fit
#' @param object SummarizedExperiment
#' @param fit 'limma', 'lme', 'lm', 'lme', 'wilcoxon'
#' @return data.table(contrast, nup, ndown)
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, subgroupvar = 'SET', impute=TRUE, plot=FALSE)
#' object %<>% fit_limma(subgroupvar = 'SET')
#' object %<>% fit_lm(subgroupvar = 'SET')
#' summarize_fit(object)
#' @export
summarize_fit <- function(object, fit = fits(object)){
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

#' Extract fit quantity
#' @param object SummarizedExperiment
#' @param quantity 'effect', 'p', 'fdr', 'bonf'
#' @return melted data.table
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'             file, invert_subgroups=inv, fit='limma', plot=FALSE)
#' extract_fit_quantity(object, fit = 'limma')
#' @noRd
extract_fit_quantity <- function(object, fit, quantity='p'){
    fvars0 <- c('feature_id','feature_name','imputed', 'control')
    fvars0 %<>% intersect(fvars(object))
    fdt <- data.table(fdata(object)[, fvars0, drop=FALSE])
    fitdt <- metadata(object)[[fit]][, , quantity, drop=FALSE]
    fitdt %<>% adrop(drop=3)
    fdt %<>% cbind(fitdt)
    melt.data.table(
        fdt, id.vars = fvars0, variable.name = 'contrast', value.name=quantity)
}


merge_fit_quantities <- function(x, y){
    names0 <- c('feature_id','feature_name','imputed', 'control', 'contrast')
    names0 %<>% intersect(names(x)) %>% intersect(names(y))
    merge(x, y, by = names0)
}


#' Extract fit datatable
#' @param object SummarizedExperiment
#' @param fit 'limma', 'lme', 'lm', 'lmer', 'wilcoxon'
#' @return melted data.table
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, fit='limma', plot=FALSE)
#' extract_fit_dt(object, fit = 'limma')
#'
#' file <- download_data('atkin18.metabolo.xlsx')
#' object <- read_proteingroups(
#'            file, invert_subgroups=inv, fit='limma', plot=FALSE)
#' extract_fit_dt(object, fit = 'limma')
#' @noRd
extract_fit_dt <- function(object, fit){
    Reduce( merge_fit_quantities,
            mapply( extract_fit_quantity,
                    quantity = c('effect', 'p', 'fdr'),
                    MoreArgs = list(object=object, fit=fit),
                    SIMPLIFY = FALSE ))
}

.summarize_fit <- function(object, fit){
    
    effect <- fdr <- . <- NULL
    extract_fit_dt(object, fit = fit)[,
        .(  ndown = sum(effect<0 & fdr<0.05, na.rm=TRUE),
            nup   = sum(effect>0 & fdr<0.05, na.rm=TRUE)),
        by='contrast']
}


old_summarize_fit <- function(
    object, 
    fit = fits(object)
){
    . <- NULL
    if (is_scalar(fit)) return(.summarize_fit(object, fit))
    res <- mapply(.summarize_fit, 
                fit = fit, 
                MoreArgs = list(object = object), 
                SIMPLIFY = FALSE)
    res %<>% mapply(function(dt, model) cbind(model = model, dt), 
                    dt = ., model = names(res), SIMPLIFY = FALSE)
    res %<>% data.table::rbindlist()
    res
}

#============================================================================
#
#                   is_sig  .is_sig
#                   plot_contrast_venn
#
#============================================================================

.is_sig <- function(object, fit, contrast, quantity = 'fdr'){
    issig  <- get(quantity)(object) < 0.05
    isdown <- issig & effect(object) < 0
    isup   <- issig & effect(object) > 0
    isdown[is.na(isdown)] <- FALSE
    isup[  is.na(isup)  ] <- FALSE
    testmat <- matrix(0, nrow(issig), ncol(issig), dimnames=dimnames(issig))
    testmat[isdown] <- -1
    testmat[isup]   <-  1
    testmat[, paste0(contrast, FITSEP, fit), drop=FALSE]
}


#' Is significant?
#' @param object    SummarizedExperiment
#' @param fit       subset of autonomics::TESTS
#' @param contrast  subset of colnames(metadata(object)[[fit]])
#' @param quantity  value in dimnames(metadata(object)[[fit]])[3]
#' @return matrix: -1 (downregulated), +1 (upregulatd), 0 (not fdr significant)
#' @examples
#' require(magrittr)
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' object %<>% fit_lm()
#' object %<>% fit_limma()
#' issig <- is_sig(object, fit = c('lm','limma'), contrast = 'Adult')
#' plot_contrast_venn(issig)
#' @export
is_sig <- function(
    object,
    fit = fits(object)[1],
    contrast = coefs(object),
    quantity = 'fdr'
){
# Assert
    . <- NULL
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_character(fit)
    assert_is_subset(fit, fits(object))
    if (is.character(contrast))  for (fi in fit){
        assert_is_subset(contrast, coefs(object, fit=fi)) }
# Run across models
    res <-  mapply(.is_sig, fit, 
                    MoreArgs = list(object = object, contrast = contrast, quantity = quantity),
                    SIMPLIFY = FALSE)
    add_model_names <- function(isfdrmat, model){
                        colnames(isfdrmat) %<>% paste(model, sep='.')
                        isfdrmat }
    if (length(fit)>1){
        res %<>% mapply(add_model_names , ., names(.), SIMPLIFY=FALSE)
    }
    res %<>% do.call(cbind, .)
    res
}


