#==============================================================================
#
#                   Switch between nondetect representations
#
#==============================================================================


#' Switch between nondetect representations
#' @param object    SummarizedExperiment
#' @param verbose   logical(1)
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#'
#' # 0 -> NA (proteinGroups LFQ intensities)
#'
#' # NaN -> NA (proteinGroups ratios)
#'     file <- load_data('stemcells.proteinGroups.txt')
#'     object <- read_proteingroups(file)
#'     nan_to_na(object, verbose=TRUE)
#'
#' # -Inf -> NA (log2 transformed proteinGroups LFQ intensity)
#'
#' # NA -> 0
#'     file <- load_data('glutaminase.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'     na_to_zero(object, verbose = TRUE)
#' @noRd
zero_to_na <- function(object, verbose = FALSE){
    selector <- exprs(object) == 0
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0('\t\tReplace 0 -> NA for %d/%d values ',
                            '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector)), nrow(object),
                    sum(colAnys(selector)), ncol(object))
        exprs(object)[selector] <- NA_real_
    }
    object
}


#' @noRd
nan_to_na <- function(object, verbose = FALSE){
    selector <- is.nan(exprs(object))
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0( '\t\tReplace NaN -> NA for %d/%d values ',
                            '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector)), nrow(object),
                    sum(colAnys(selector)), ncol(object))
        exprs(object)[selector] <- NA_real_
    }
    object
}


na_to_zero <- function(object, verbose = FALSE){
    selector <- is.na(exprs(object))
    if (any(selector)){
        if (verbose) cmessage(
                        paste0( '\t\tReplace NA -> 0 for %d/%d values ',
                                '(in %d/%d features and %d/%d samples)'),
                        sum(selector), nrow(selector)*ncol(selector),
                        sum(rowAnys(selector)), nrow(object),
                        sum(colAnys(selector)), ncol(object))
        exprs(object)[selector] <- 0
    }
    object
}


inf_to_na <- function(object, verbose){
    selector <- is.infinite(exprs(object))
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0(
                        '\t\tReplace -Inf -> NA for %d/%d values ',
                        '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector)), nrow(object),
                    sum(colAnys(selector)), ncol(object))
        exprs(object)[selector] <- NA_real_
    }
    object
}


minusinf_to_na <- function(object, verbose = FALSE){
    selector <- exprs(object)==-Inf
    if (any(c(selector), na.rm = TRUE)){
        if (verbose) cmessage(
                    paste0( '\t\tReplace -Inf -> NA for %d/%d values ',
                            '(in %d/%d features and %d/%d samples)'),
                    sum(selector, na.rm=TRUE), nrow(selector)*ncol(selector),
                    sum(rowAnys(selector)), nrow(object),
                    sum(colAnys(selector)), ncol(object))
        exprs(object)[selector] <- NA_real_
    }
    object
}


na_to_string <- function(x){
    x[is.na(x)] <- ''
    x
}


#============================================================================
#
#                            impute_common_nas
#
#=============================================================================


impute_common_nas <- function(
    object, imputefun = function(x) impute_around_zero(x), verbose = FALSE
){
    # Assert no NaN or -Inf (nondetects should be represented as NA)
    assert_any_are_not_nan(exprs(object))
    assert_all_are_finite(na.exclude(c(exprs(object))))

    # Typify nondetects
    is_na <- is.na(exprs(object))
    is_consistent_na <- rowAlls(is_na)
    if (verbose)  cmessage(
                "\t\tImpute %d/%d features with NA value in all %d samples",
                sum(is_consistent_na), length(is_consistent_na), ncol(object))
    is_inconsistent_na <- is_na &
                        !matrix(is_consistent_na, nrow(object), ncol(object))

    # Replace inconsistent nondetects by median
    exprs1 <- exprs(object)
    for (sample in seq_len(ncol(object))){
        exprs1[, sample] %<>% (function(x){
            x[is_inconsistent_na[, sample]] <- stats::median(x, na.rm=TRUE)
            x})
    }

    # Impute consistent nondetects
    is_imputed(object) <- is.na(exprs1)
    exprs1 %<>% imputefun()

    # Re-NA inconsistent nondetects
    exprs1[is_inconsistent_na] <- NA_real_

    # Update & Return
    exprs(object) <- exprs1
    object
}

#' Impute around zero
#' @param x exprs matrix
#' @return exprs matrix
#' @examples
#' # Example data
#' #-------------
#' x <- matrix(c(20, 21, 22, 23,
#'               30, 31, 32, 33,
#'               NA, NA, NA, NA,
#'               40, 41, NA, NA,
#'               NA, NA, NA, NA,
#'               50, 51, NA, NA), ncol=4, byrow=TRUE)
#' object <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs=x))
#' object$subgroup <- c('EV', 'EV', 'MS', 'MS')
#'
#'# Impute
#'#-------
#' x
#' impute_around_zero(x)
#'
#' exprs(object)
#' exprs(impute_consistent_nas(object, imputefun = impute_around_zero))
#' exprs(impute_common_nas(object, imputefun = impute_around_zero))
#' @noRd
impute_around_zero <- function(x){
    if (ncol(x)==1) return(x) # single replicate -> no sd!
    meansd <- mean(rowSds(x), na.rm = TRUE)
    is_common_na_feature <- rowAlls(is.na(x))
    imputed_values  <-  x[is_common_na_feature, , drop=FALSE]       %>%
                        apply(1, function(y) impute_row(y, meansd)) %>%
                        t()
    x[is_common_na_feature, ] <- imputed_values
    x
}

impute_row <- function(x, sd){
    abs(rnorm(length(x), sd = sd * 2))
}


#=============================================================================
#
#                       impute_consistent_nas
#
#=============================================================================

#' Impute consistent NA values
#'
#' Impute values missing in all (subgroup) samples.
#'
#' @param object     SummarizedExperiment
#' @param imputefun  imputation function
#' @param svar       string
#' @param verbose    TRUE or FALSE
#' @examples
#' # Read object
#'     file <- download_data('glutaminase.metabolon.xlsx')
#'     object <- read_metabolon(file)
#'
#' # Common NA values - missing in all samples
#'    impute_common_nas(split_by_svar(object)[[2]], verbose=TRUE)
#'
#' # Consistent NA values - missing in all subgroup samples
#'    impute_consistent_nas(object, verbose = TRUE)
#'
#' # Use a different imputation method
#'    impute_consistent_nas(object, verbose = TRUE,
#'       imputefun = function(x) imputeLCMD::impute.QRILC(x)[[1]])
#' @return SummarizedExperiment with updated exprs
#' @noRd
impute_consistent_nas <- function(
    object, imputefun = function(x) impute_around_zero(x), svar = 'subgroup',
    verbose      = FALSE
){
    # Check
    . <- NULL
    if (!svar %in% svars(object)){
        cmessage("\t\tNo svar '%s' => no imputation performed", svar)
        return(object)
    }
    if (any(is_missing_or_empty_character(svalues(object, svar)))){
        cmessage(
            "\t\tSome '%s' values are missing => no imputation performed", svar)
        return(object)
    }
    selector <- is_available_for_some_feature(object)
    if (!all(selector)) stop('First rm ',
                            collapse_words(names(selector)[!selector]),
                            ', which contain only NA/0 values')

    # Impute
    imputed_object <- split_by_svar(object, svar) %>%
                        lapply(impute_common_nas, imputefun=imputefun) %>%
                        do.call(SummarizedExperiment::cbind, .)
    imputed_object %<>% extract(, snames(object))

    # Message
    nimputed <- sum(rowAnys(is_imputed(imputed_object)))
    if (verbose & nimputed>0) cmessage(
                    '\t\tImpute consistent NA values among %d/%d features',
                    nimputed, nrow(object))
    # Return
    return(imputed_object)
}

is_missing_or_empty_character <- function(x){
    x[is.na(x)] <- ''
    x == ''
}


#' Split SummarizedExperiment by svar
#' @param object SummarizedExperiment
#' @param svar 'subgroup'
#' @return list of SummarizedExperiments
#' @examples
#' file <- download_data('stemcells.proteinGroups.txt')
#' object <- read_proteingroups(file)
#' split_by_svar(object)
#' @noRd
split_by_svar <- function(object, svar = 'subgroup'){
# Return object if null svar
    if (is.null(svar)) return(list(object))
# Split
    extract_samples <- function(sg){
                        idx <- sdata(object)[[svar]] == sg
                        object[, idx]
                    }
    Map(extract_samples, slevels(object, svar))
}



# Linguistic collapse
# x <- c('a', 'b', 'c')
# collapse_words(x)
collapse_words <- function(x, collapsor = 'and'){
    if (length(x) == 1) return(x)
    if (length(x) == 2) return(sprintf('%s %s %s', x[[1]], collapsor, x[[2]]))
    paste0(x[-length(x)], collapse = ', ') %>% paste0(' and ', x[[3]])
}


