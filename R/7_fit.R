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
subgroup_array <- function(object, subgroupvar){
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
#' subgroup_matrix(object)
#' @export
subgroup_matrix <- function(object, subgroupvar){
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
#                    limma & limma<-
#                    extract_fit_dt
#
#==============================================================================

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
    assert_is_subset(fit, TESTS)
    fvars0 <- c('feature_id','feature_name','imputed')
    fvars0 %<>% intersect(fvars(object))
    fdt <- data.table(fdata(object)[, fvars0, drop=FALSE])
    fitdt <- metadata(object)[[fit]][, , quantity, drop=FALSE]
    fitdt %<>% adrop(drop=3)
    fdt %<>% cbind(fitdt)
    data.table::melt.data.table(
        fdt, id.vars = fvars0, variable.name = 'contrast', value.name=quantity)
}


merge_fit_quantities <- function(x, y){
    names0 <- c('feature_id','feature_name','imputed', 'contrast')
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
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(
#'            file, invert_subgroups=inv, fit='limma', plot=FALSE)
#' extract_fit_dt(object, fit = 'limma')
#' @noRd
extract_fit_dt <- function(object, fit){
    Reduce( merge_fit_quantities,
            mapply( extract_fit_quantity,
                    quantity = c('effect', 'p', 'fdr', 'bonf'),
                    MoreArgs = list(object=object, fit=fit),
                    SIMPLIFY = FALSE ))
}

#' Extract fit summary
#' @param object SummarizedExperiment
#' @param fit 'limma', 'lme', 'lm', 'lme', 'wilcoxon'
#' @return data.table(contrast, nup, ndown)
#' @examples
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, fit='limma', plot=FALSE)
#' summarize_fit(object, 'limma')
#' @export
summarize_fit <- function(object, fit){
    effect <- fdr <- NULL
    extract_fit_dt(object, fit = fit)[,
        .(  ndown = sum(effect<0 & fdr<0.05, na.rm=TRUE),
            nup   = sum(effect>0 & fdr<0.05, na.rm=TRUE)),
        by='contrast']
}


#============================================================================
#
#                   is_fdr  .is_fdr
#                   plot_venn
#
#============================================================================

.is_fdr <- function(object, fit, contrast){
    fitres <- metadata(object)[[fit]]
    isfdr  <- adrop(fitres[, , 'fdr', drop = FALSE] < 0.05, 3)
    isdown <- isfdr & adrop(fitres[, , 'effect', drop=FALSE]<0, 3)
    isup   <- isfdr & adrop(fitres[, , 'effect', drop=FALSE]>0, 3)
    testmat <- matrix(0, nrow(isfdr), ncol(isfdr), dimnames=dimnames(isfdr))
    testmat[isdown] <- -1
    testmat[isup]   <-  1
    testmat[, contrast, drop=FALSE]
}

#' @rdname plot_venn
#' @export
is_fdr <- function(
     object,
     fit = intersect(names(metadata(object)), TESTS),
     contrast = if (is_scalar(fit)) colnames(metadata(object)[[fit]]) else 1
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    assert_is_character(fit)
    assert_is_subset(fit, names(metadata(object)))
    if (is.character(contrast))  for (fi in fit)  assert_is_subset(
                        contrast, colnames(metadata(object)[[fi]]))
# Run across models
    res <-  mapply(.is_fdr, fit, 
                    MoreArgs = list(object=object, contrast=contrast),
                    SIMPLIFY = FALSE)
    add_model_names <- function(isfdrmat, model){
                        colnames(isfdrmat) %<>% paste(model, sep='.')
                        isfdrmat }
    res %<>% mapply(add_model_names , ., names(.), SIMPLIFY=FALSE)
    res %<>% do.call(cbind, .)
    res
}


#' Plot venn
#' @param isfdr matrix(nrow, ncontrast): -1 (down), +1 (up)
#' @return NULL
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' object %<>% fit_wilcoxon() # fails on unpaired imputed data
#' object %<>% fit_lm()
#' object %<>% fit_limma()
#' isfdr <- is_fdr(object, fit = c('lm','limma'), contrast = 'X30dpt-Adult')
#' plot_venn(isfdr)
#' @export
plot_venn <- function(isfdr){
    layout(matrix(1:2, nrow=2))
    vennDiagram(isfdr, include='up',   mar = rep(0,4), show.include=TRUE)
    vennDiagram(isfdr, include='down', mar = rep(0,4), show.include=TRUE)
}
    
