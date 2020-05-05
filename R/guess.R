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

#' All elements of vector are identical
#' @param x vector
#' @return TRUE or FALSE
#' @examples
#' x <- c(2,2,1,2)
#' has_identical_values(x)
#' @noRd
has_identical_values <- function(x) length(unique(x))==1

#' Is maximal
#' @param x numeric vector
#' @return logical vector
#' @examples
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' is_max(x)
#' @noRd
is_max <- function(x) cequals(x, max(x, na.rm = TRUE))


#=============================================================================
#' Guess separator
#' @param x                   character vector or SummarizedExperiment
#' @param var                 svar or fvar
#' @param possible_separators character vector: possible separators to look for
#' @param verbose             TRUE or FALSE
#' @param ...                 used for proper S3 method dispatch
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
#'    # file <- download_autonomics_data('glutaminase_metab.xlsx')
#'    # object <- read_metabolon(file)
#'    # guess_sep(object)
#'
#'    # file <- download_autonomics_data('stemcells_proteinGroups.txt')
#'    # object <- read_proteingroups(object)
#'    # guess_sep(object)
#' @export
guess_sep <- function (x, ...) {
   UseMethod("guess_sep", x)
}


#' @rdname guess_sep
#' @export
guess_sep.character <- function(
   x,
   possible_separators = c('.', ' ', '_'),
   verbose = FALSE,
   ...
){
   . <- NULL
   sep_freqs <- Map(function(y) stri_split_fixed(x, y), possible_separators) %>%
                lapply(function(y) vapply(y, length, integer(1)))            %>%
                extract( vapply(., has_identical_values, logical(1)))        %>%
                vapply(unique, integer(1))

   # No separator detected - return NULL
   if (all(sep_freqs==1)){
      if (verbose)  message(x[1], ': no (consistent) separator. Returning NULL')
      return(NULL)   # no separator detected
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
   if (verbose) cmessage("\t\tGuess sep: '", best_sep, "'")
   return(best_sep)
}


#' @rdname guess_sep
#' @export
guess_sep.factor <- function(x, ...)  guess_sep.character(levels(x))


#' @rdname guess_sep
#' @importFrom magrittr %>%
#' @export
guess_sep.SummarizedExperiment <- function(
   x,
   var = 'sample_id',
   possible_separators = c('.', '_', ' '),
   # if (contains_ratios(x)) c('.', ' ') else c('.', '_', ' '),
   verbose = FALSE,
   ...
){
   assert_is_subset(var, c(svars(x), fvars(x)))
  (if (var %in% svars(x)) slevels(x, var) else flevels(x, var)) %>%
   guess_sep(possible_separators = possible_separators,
             verbose             = verbose)
}



#=======================================================================

extract_first_components <- function(x, sep){
   x                          %>%
   stri_split_fixed(sep)      %>%
   vapply(function(y) y       %>%
   extract(1:(length(y)-1))   %>%
   paste0(collapse = sep), character(1))
}

extract_last_component <- function(x, sep){
   x                          %>%
   stri_split_fixed(sep)      %>%
   vapply(function(y) y       %>%
   extract(length(y))         %>%
   paste0(collapse = sep), character(1))
}


#' Guess subgroup values
#' @param x       character vector or SummarizedExperiment
#' @param sep     string
#' @param invert  FALSE (default) or TRUE: guess "non-subgroup" component?
#' @param verbose logical(1)
#' @param ...     used for proper S3 method dispatch
#' @return character(n)
#' @examples
#' require(magrittr)
#'
#' # charactervector
#'    # No sep: subgroup = x
#'       x <- c("EM00", "EM01", "EM02")
#'       guess_subgroup_values(x)
#'
#'    # Sep: subgroup = head components of x
#'       x <- c("UT_10h_R1", "UT_10h_R2", "UT_10h_R3")
#'       guess_subgroup_values(x)
#'       guess_subgroup_values(x, invert = TRUE)
#'
#'       x <- c("EM00_STD.R1", "EM01_STD.R1", "EM01_EM00.R1")
#'       guess_subgroup_values(x)
#'       guess_subgroup_values(x, invert = TRUE)
#'
#' @export
guess_subgroup_values <- function (x, ...) {
   UseMethod("guess_subgroup_values", x)
}

#' @rdname guess_subgroup_values
#' @export
guess_subgroup_values.character <- function(
   x,
   sep     = guess_sep(x),
   invert  = FALSE,
   verbose = FALSE,
   ...
){
   # Guess
   subgroup_values <- if (is.null(sep)){  x
                      } else if (invert){ extract_last_component(x, sep)
                      } else {            extract_first_components(x, sep)
                      }
   # Inform
   if (verbose)   message(
      '\t\tGuess subgroup values: ', x[1], ' => ', subgroup_values[1])

   # Return
   return(subgroup_values)
}


#' @rdname guess_subgroup_values
#' @export
guess_subgroup_values.SummarizedExperiment <- function(
   x,
   sep      = guess_sep(x),
   invert   = FALSE,
   verbose  = FALSE,
   ...
){

   # already in x
   if ('subgroup' %in% svars(x)){
      if (verbose) message("\t\tUse 'subgroup' values in x ")
      return(sdata(x)$subgroup)
   }

   # guess from sampleid values
    guess_subgroup_values(
       sampleid_values(x), sep = sep, invert = invert, verbose = verbose)
}


#==========================================================

# guess_subject_values <- function (x, ...) {
#    UseMethod("guess_subject_values", x)
# }
#
# guess_subject_values.character(
#    x,
#    sep     = guess_sep(x),
#    verbose = FALSE
# ){
#    NULL
# }


