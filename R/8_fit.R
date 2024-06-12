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
#' x <- c(A = 1, B = 3, C = 2, D = 3, E = NA)
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
#' x <- c(A = 1,B = 3,C = 2,D = 3, E = NA)
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
#'    guess_sep(c('PERM_NON.R1[H/L]', 'PERM_NON.R2[H/L]'))
#'    guess_sep(c('WT_untreated_1', 'WT_untreated_2'))
#'    guess_sep(c('group1', 'group2.R1'))
#' # SummarizedExperiment
#'    file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'    object <- read_metabolon(file)
#'    guess_sep(object)
#' @export
guess_sep <- function (x, ...)  UseMethod("guess_sep", x)


#' @rdname guess_sep
#' @export
guess_sep.numeric <- function(x, ...) NULL


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
#' @param x    character vector
#' @param sep  string
#' @param i    integer
#' @return character vector
#' @examples
#' # Read
#'     file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'     object <- read_metabolon(file)
#'     x <- object$sample_id[1:5]
#'     nfactors(x)
#' # Split
#'     split_extract_fixed(x, '.', 1:2)
#'     split_extract_fixed(x, '.', seq_len(nfactors(x)-1))
#'     split_extract_fixed(x, '.', nfactors(x))
#'     split_extract_fixed(fdt(object)$PUBCHEM, ';', 1) # with NA values
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
#                   default_sfile
#
#=============================================================================

#' Default sfile
#' @param file data file
#' @return sample file
#' @examples
#' file <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
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
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' object$subgroup <- paste0(object$Diabetes, '.', object$subgroup)
#' subgroup_matrix(object, 'subgroup')
#' @export
subgroup_matrix <- function(object, subgroupvar){
    . <- NULL
    subgroup_array <- subgroup_array(object, subgroupvar)
    if (length(dim(subgroup_array)) == 1){
        return(matrix(subgroup_array, 
                      byrow = TRUE, 
                      nrow = 1, 
                      dimnames = list(NULL, subgroup_array)))  
    }
    otherdims <- names(dim(subgroup_array)) %>% setdiff('V1')
    ncol1   <- Reduce('*', dim(subgroup_array)[otherdims])
    colnames1 <- dimnames(subgroup_array)[otherdims] %>%
                 expand.grid()                       %>%
                 apply(1, paste0, collapse = '.')
    subgroupmat <- matrix(subgroup_array,
                        nrow = nrow(subgroup_array), ncol = ncol1,
                        dimnames=list(rownames(subgroup_array), colnames1))
    subgroupmat %>% extract(nrow(.):1, )
    #dt <- split_subgroup_levels(object)
    #subgroupmat <- as.matrix(data.table::dcast(
    #    dt, V1 ~ V2, value.var = 'subgroup'), rownames = 'V1')
    #subgroupmat %>% extract(rev(order(rownames(.))), order(colnames(.)))
}

#------------------------------------------------------------------------------
#
#                   modelvar
#
#------------------------------------------------------------------------------


#' Get model variable
#' @param object          data.table or SummarizedExperiment
#' @param quantity        'p', 'effect', 'fdr', 't', or 'se'
#' @param fit             string (vector)
#' @param coef            string (vector)
#' @param fvar            'feature_id' or other fvar for values (pvec) or names (upfeatures)
#' @param significancevar 'p' or 'fdr'
#' @param significance     p or fdr cutoff (fractional number)
#' @param effectdirection  '<>', '<' or '>'
#' @param effectsize      effectsize cutoff (positive number)
#' @param ...             S3 dispatch
#' @return string (tvar), matrix (tmat), numeric vector (tvec), character vector (tfeatures)
#' @examples 
#'     file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'     object <- read_metabolon(file)
#'     object %<>% fit_limma(statvars = c('effect', 't', 'p'))
#'     object %<>% fit_lm(   statvars = c('effect', 't', 'p'))
#' 
#'     effectvar(object)
#'     effectvec(object)[1:3]
#'      effectdt(object)[1:3, ]
#'     effectmat(object)[1:3, ]
#' 
#'          tvar(object)
#'          tvec(object)[1:3]
#'           tdt(object)[1:3, ]
#'          tmat(object)[1:3, ]
#' 
#'          pvar(object)
#'          pvec(object)[1:3]
#'           pdt(object)[1:3, ]
#'          pmat(object)[1:3, ]
#' 
#' modelfeatures(object)
#'  downfeatures(object)
#'    upfeatures(object)
#' @export    
modelvar <- function(object, ...) UseMethod('modelvar')

#' @rdname modelvar
#' @export
modelvar.data.table <- function(
    object, quantity, fit = fits(object), coef = default_coefs(object, fit = fit), ...
){
# Assert
    assert_is_subset(quantity, c('fdr', 'p', 't', 'effect', 'se', 'abstract'))
    assert_is_subset(fit,   fits(object))
    assert_is_subset(coef, coefs(object, fit = fit))
    assert_has_no_duplicates(fit)   # Avoid duplicate columns
    assert_has_no_duplicates(coef)  # Downstream functionality melts and dcasts
# Return                            # Which requires the cast to be unique
    sep <- guess_fitsep(object)  # Else length rather than data values are entered
    x <- expand.grid(quantity = quantity, fit = fit, coef = coef)
    x <-  paste(x$quantity, x$coef, x$fit, sep = sep)
    x %<>% intersect(names(object))     # fits dont always contain same coefs: 
    if (length(x)==0)  x <- NULL           # `limma(contrasts)` mostly without intercept
    x   # NULL[1] and c('a', NULL) work!   # `lm(coefs)` mostly with intercept               
}


#' @rdname modelvar
#' @export
modelvar.SummarizedExperiment <- function(
    object, quantity, fit = fits(object), coef = default_coefs(object, fit = fit), ...
){  
    modelvar.data.table(fdt(object), quantity = quantity, fit = fit, coef = coef )
}


#' @rdname modelvar
#' @export
effectvar <- function( 
    object, fit = fits(object), coef = default_coefs(object, fit = fit)
){
    modelvar(object, quantity = 'effect', fit = fit, coef = coef)
}


#' @rdname modelvar
#' @export
tvar <- function(
    object, fit = fits(object), coef = default_coefs(object, fit = fit)
){
    modelvar(object, quantity = 't', fit = fit, coef = coef)
}


#' @rdname modelvar
#' @export
pvar <- function( object, fit = fits(object), coef = default_coefs(object, fit = fit) ){
    modelvar(object, quantity = 'p', fit = fit, coef = coef )
}


#' @rdname modelvar
#' @export
fdrvar <- function( object, fit = fits(object), coef = default_coefs(object, fit = fit) ){
    modelvar(object, quantity = 'fdr', fit = fit, coef = coef) 
}


#' @rdname modelvar
#' @export
abstractvar <- function(object, ...)  UseMethod('abstractvar')


#' @rdname modelvar
#' @export
abstractvar.data.table <- function(
    object, fit = fits(object), coef = default_coefs(object, fit = fit), ...
){
    sep <- guess_fitsep(object)             # cant use modelvar because its
    y <- paste(  rep(coef, length(fit)),         # t1~limma
                 rep(fit, each = length(coef)),
                 sep = sep )
    if (y %in% names(object))  y  else NULL
}

#' @rdname modelvar
#' @export 
abstractvar.SummarizedExperiment <- function(
    object, fit = fits(object), coef = default_coefs(object, fit = fit), ...
){
    abstractvar.data.table(fdt(object), fit = fit, coef = coef)
}


#------------------------------------------------------------------------------
#
#                   modelvec
#
#------------------------------------------------------------------------------

#' @rdname modelvar
#' @export
modelvec <- function(object, ...)  UseMethod('modelvec')


#' @rdname modelvar
#' @export
modelvec.data.table <- function(
    object, quantity, fit = fits(object)[1], coef = default_coefs(object, fit = fit)[1], 
    fvar = 'feature_id', ...
){
    valuevar <- modelvar(object, quantity = quantity, fit = fit, coef = coef)
    if (is.null(valuevar))  return(NULL)
    y        <- object  %>%  extract2(valuevar)
    names(y) <- object  %>%  extract2(fvar)
    y
}


#' @rdname modelvar
#' @export
modelvec.SummarizedExperiment <- function(
    object, quantity, fit = fits(object)[1], coef = default_coefs(object, fit = fit)[1], 
    fvar = 'feature_id', ...
){
    modelvec.data.table(fdt(object), quantity = quantity, fit = fit, coef = coef, fvar = fvar)
}


#' @rdname modelvar
#' @export
effectvec <- function(
    object, fit = fits(object)[1], coef = default_coefs(object)[1], fvar = 'feature_id'
){
    modelvec(object, quantity = 'effect', fit = fit, coef = coef, fvar = fvar)
}

#' @rdname modelvar
#' @export
tvec <- function(
    object, fit = fits(object)[1], coef = default_coefs(object, fit = fit)[1], fvar = 'feature_id'
){
    modelvec(object, quantity = 't', fit = fit, coef = coef, fvar = fvar)
}


#' @rdname modelvar
#' @export
pvec <- function(
    object, fit = fits(object)[1], coef = default_coefs(object, fit = fit)[1], fvar = 'feature_id'
){
    modelvec(object, quantity = 'p', fit = fit, coef = coef, fvar = fvar)
}


#' @rdname modelvar
#' @export
fdrvec <- function(
     object, fit = fits(object)[1], coef = default_coefs(object, fit = fit)[1],  fvar = 'feature_id'
){
    modelvec(object, quantity = 'fdr', fit = fit, coef = coef, fvar = fvar)
}


#' @rdname modelvar
#' @export
abstractvec <- function(object, ...)  UseMethod('abstractvec')


#' @rdname modelvar
#' @export
abstractvec.data.table <- function(
    object, fit = fits(object)[1], coef = default_coefs(object, fit = fit)[1], fvar = 'feature_id', ...
){
    var <- abstractvar(object, fit = fit, coef = coef)
    if (is.null(var))  return(NULL)
    object[[var]]
}


#' @rdname modelvar
#' @export
abstractvec.SummarizedExperiment <- function(
    object, fit = fits(object)[1], coef = default_coefs(object, fit = fit)[1], fvar = 'feature_id', ...
){
    y <- abstractvec.data.table(fdt(object), fit = fit, coef = coef, fvar = fvar)
    names(y) <- fdt(object)[[fvar]]
    y
}



#------------------------------------------------------------------------------
#
#                   modeldt
#
#------------------------------------------------------------------------------


#' @rdname modelvar
#' @export
modeldt <- function(object, ...)  UseMethod('modeldt')


#' @rdname modelvar
#' @export
modeldt.data.table <- function(
    object, quantity, fit = fits(object), coef = default_coefs(object, fit = fit), ...
){
    var <- modelvar(object, quantity, coef = coef, fit = fit)
    if (is.null(var))  return(NULL)
    dt <- object[, c('feature_id', var), with = FALSE]
    #names(dt) %<>% stri_replace_first_fixed(paste0(var, sep), '')
    dt
}

#' @rdname modelvar
#' @export
modeldt.SummarizedExperiment <- function(
    object, quantity, fit = fits(object), coef = default_coefs(object, fit = fit), ...
){
    modeldt.data.table(fdt(object), quantity = quantity, fit = fit, coef = coef)
}


#' @rdname modelvar
#' @export
effectdt <- function( object, quantity, fit = fits(object), coef = default_coefs(object, fit = fit) ){
    modeldt(object, quantity = 'effect', fit = fit, coef = coef)
}


#' @rdname modelvar
#' @export
tdt <- function( object, quantity, fit = fits(object), coef = default_coefs(object, fit = fit) ){
    modeldt(object, quantity = 't', fit = fit, coef = coef)
}


#' @rdname modelvar
#' @export
pdt <- function( object, quantity, fit = fits(object), coef = default_coefs(object, fit = fit)){
    modeldt(object, quantity = 'p', fit = fit, coef = coef)
}



#------------------------------------------------------------------------------
#
#                   modelmat
#
#------------------------------------------------------------------------------

#' @rdname modelvar
#' @export
modelmat <- function(object, ...)  UseMethod('modelmat')

#' @rdname modelvar
#' @export
modelmat <- function(
    object, quantity, fit = fits(object), coef = default_coefs(object, fit = fit)
){
    dt2mat(modeldt(object, quantity = quantity, fit = fit, coef = coef))
}


#' @rdname modelvar
#' @export
effectmat <- function(object, fit = fits(object), coef = default_coefs(object, fit = fit)){
    modelmat(object, quantity = 'effect', fit = fit, coef = coef)
}


#' @rdname modelvar
#' @export
effectsizemat <- function(object, fit = fits(object), coef = default_coefs(object, fit = fit)){
  # dont rm: used in ..extract_statistic_features : 
  # getFromNamespace(sprintf('%smat', statistic), 'autonomics')
    abs(modelmat(object, quantity = 'effect', fit = fit, coef = coef))
}


#' @rdname modelvar
#' @export
tmat <- function(object, fit = fits(object), coef = default_coefs(object, fit = fit)){
    modelmat(object, quantity = 't', fit = fit, coef = coef)
}


#' @rdname modelvar
#' @export
pmat <- function(object, fit = fits(object), coef = default_coefs(object, fit = fit)){
    modelmat(object, quantity = 'p', fit = fit, coef = coef)
}


#' @rdname modelvar
#' @export
fdrmat <- function(object, fit = fits(object), coef = default_coefs(object, fit = fit)){
    modelmat(object, quantity = 'fdr', fit = fit, coef = coef)
}


#==============================================================================
#
#                   modelfeatures
#
#==============================================================================

#' @rdname modelvar
#' @export
modelfeatures <- function(object, ...) UseMethod('modelfeatures')


#' @rdname modelvar
#' @export
modelfeatures.data.table <- function(
             object,
                fit = fits(object)[1],
               coef = default_coefs(object, fit = fit)[1], 
               fvar = 'feature_id', 
    significancevar = 'p',
       significance = 0.05,
    effectdirection = '<>',
         effectsize = 0, 
    ...
){
    significancevalues <- modelvec(object, quantity =  significancevar, fit = fit, coef = coef, fvar = fvar)
          effectvalues <- modelvec(object, quantity = 'effect',         fit = fit, coef = coef, fvar = fvar)
    idx <- switch(effectdirection, 
                  `<>` = (effectvalues < -effectsize) | (effectvalues > +effectsize), 
                  `<`  = (effectvalues < -effectsize), 
                  `>`  = (effectvalues > +effectsize) )
    idx <- idx & (significancevalues < significance)
    y <- object[[fvar]][idx]
    unique(y)
}


#' @rdname modelvar
#' @export
modelfeatures.SummarizedExperiment <- function(object, ...)   modelfeatures.data.table(fdt(object), ...)


#' @rdname modelvar
#' @export
upfeatures <- function(
             object, 
                fit = fits(object)[1], 
               coef = default_coefs(object, fit = fit)[1], 
               fvar = 'feature_id',
    significancevar = 'p',
       significance = 0.05,
         effectsize = 0
){
    modelfeatures(    object = object, 
                            fit = fit,
                           coef = coef,
                           fvar = fvar, 
                significancevar = significancevar,
                   significance = significance, 
                effectdirection = '>',
                     effectsize = effectsize )
}

#' @rdname modelvar
#' @export
downfeatures <- function(
             object,
                fit = fits(object)[1], 
               coef = default_coefs(object, fit = fit)[1], 
               fvar = 'feature_id',
    significancevar = 'p',
       significance = 0.05,
         effectsize = 0
){
    modelfeatures(  object = object,
                          fit = fit,
                         coef = coef,
                         fvar = fvar, 
              significancevar = significancevar,
                 significance = significance, 
              effectdirection = '<', 
                   effectsize = effectsize )
}


#----------------------------------------------------------------------------------------


# dont rm - its the lower-level function used by fits() and coefs() !
.effectvars <- function(featuredt){
    . <- NULL
    sep <- guess_fitsep(featuredt)
    names(featuredt) %>% extract(stri_startswith_fixed(., paste0('effect', sep)))
}


#' Get fit models
#' 
#' @param object SummarizedExperiment or data.table
#' @param ...    S3 dispatch
#' @return  character vector
#' @examples 
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file, fit = 'limma')
#' fits(object)
#' @export
fits <- function(object, ...)  UseMethod('fits')


#' @rdname fits
#' @export
fits.data.table <- function(object, ...){
    sep <- guess_fitsep(object)
    if (is.null(sep))  return(NULL)
    x <- .effectvars(object)
    x %<>% split_extract_fixed(sep, 3)
    x %<>% unique()
    x
}

#' @rdname fits
#' @export
fits.SummarizedExperiment <- function(object, ...){
    fits.data.table(fdt(object))
}

#' Get coefs
#' 
#' @param object  factor, data.table, SummarizedExperiment
#' @param fit     'limma', 'lm', 'lme', 'lmer', 'wilcoxon'
#' @param svars   NULL or charactervector (svar for which to return coefs)
#' @param ...     required for s3 dispatch
#' @return  character vector
#' @examples
#' # Factor
#'     x <- factor(c('A', 'B', 'C'))
#'     coefs(x)
#'     coefs(code(x, contr.treatment.explicit))
#'     coefs(code(x, code_control))
#'     
#' # SummarizedExperiment
#'     file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'     object <- read_metabolon(file, fit = 'limma')
#'     coefs(object)
#' @export
coefs <- function(object, ...)  UseMethod('coefs')

#' @rdname coefs
#' @export
coefs.factor <- function(object, ...)   colnames(contrasts(object))

#' @rdname coefs
#' @export
coefs.data.table <- function(
    object, fit = fits(object), svars = NULL, ...
){
    sep <- guess_fitsep(object)
    if (is.null(sep))  return(NULL)
    if (is.null(fit))  return(NULL)
    . <- NULL
    coefs0 <- split_extract_fixed(.effectvars(object), sep, 2)
    fits0  <- split_extract_fixed(.effectvars(object), sep, 3)
    coefs0 %<>% extract(fits0 %in% fit)
    coefs0 %<>% unique()
    #if (!is.null(svars))  coefs0 %<>% extract(Reduce('|', lapply(svars, grepl, .)))
    coefs0 
}

#' @rdname coefs
#' @export
coefs.SummarizedExperiment <- function(object, fit = fits(object), ...){
    coefs.data.table(fdt(object), fit = fit)
}

#============================================================================
#
#                   is_sig  .is_sig
#                   plot_contrast_venn
#
#============================================================================

.is_sig <- function(object, fit, contrast, quantity = 'fdr'){
    sigfun <- paste0(quantity, 'mat')
    issig  <- get(sigfun)(object) < 0.05
    isdown <- issig & effectmat(object) < 0
    isup   <- issig & effectmat(object) > 0
    isdown[is.na(isdown)] <- FALSE
    isup[  is.na(isup)  ] <- FALSE
    testmat <- matrix(0, nrow(issig), ncol(issig), dimnames=dimnames(issig))
    testmat[isdown] <- -1
    testmat[isup]   <-  1
    col <- modelvar(object, fit = fit, coef = contrast, quantity = quantity)
    testmat[, col, drop = FALSE]
}


#' Is significant?
#' @param object    SummarizedExperiment
#' @param fit       subset of autonomics::TESTS
#' @param contrast  subset of colnames(metadata(object)[[fit]])
#' @param quantity  value in dimnames(metadata(object)[[fit]])[3]
#' @return matrix: -1 (downregulated), +1 (upregulatd), 0 (not fdr significant)
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' object %<>% fit_lm()
#' object %<>% fit_limma()
#' issig <- is_sig(object, fit = c('lm','limma'), contrast = 'Adult')
#' plot_contrast_venn(issig)
#' @export
is_sig <- function(
      object,
         fit = fits( object)[1],
    contrast = coefs(object),
    quantity = 'fdr'
){
# Assert
    . <- NULL
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_character(fit)
    assert_is_subset(fit, fits(object))
    if (is.character(contrast))  for (fi in fit)  assert_is_subset(contrast, coefs(object, fit = fi))
# Run across models
    if (quantity == 'fdr')  fdt(object) %<>% add_adjusted_pvalues('fdr')
    res <-  mapply(.is_sig, fit, 
                    MoreArgs = list(object = object, contrast = contrast, quantity = quantity),
                    SIMPLIFY = FALSE)
    #add_model_names <- function(isfdrmat, model){
    #                    colnames(isfdrmat) %<>% paste(model, sep='.')
    #                    isfdrmat }
    #if (length(fit)>1){
    #    res %<>% mapply(add_model_names , ., names(.), SIMPLIFY=FALSE)
    #}
    res %<>% do.call(cbind, .)
    res
}

#' Linear Modeling Engines
#' @examples
#' LINMOD_ENGINES
#' @export
LINMOD_ENGINES <- c('limma', 'lm', 'lme', 'lmer', 'wilcoxon')

