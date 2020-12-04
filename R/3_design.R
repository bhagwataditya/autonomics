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
#'    # object <- read_metabolon(file)
#'    # guess_sep(object)
#'
#'    # file <- download_data('billing16.proteingroups.txt')
#'    # object <- read_proteingroups(object)
#'    # guess_sep(object)
#' @export
guess_sep <- function (x, ...) {
    UseMethod("guess_sep", x)
}


#' @rdname guess_sep
#' @export
guess_sep.character <- function(
    x, separators = c('.', ' ', '_'), verbose = FALSE, ...
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
    x, var = 'sample_id', separators =  c('.', '_', ' '),
    verbose = FALSE, ...
){
    assert_is_subset(var, c(svars(x), fvars(x)))
    (if (var %in% svars(x)) slevels(x, var) else flevels(x, var)) %>%
    guess_sep(separators = separators, verbose = verbose)
}


#=============================================================================
#
#                   guess_subgroup_values
#                   guess_replicate_values
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

#' Guess subgroup/replicate values
#' @param x       sampleid values (string vector)
#' @param sep     subfactor separator (string)
#' @param verbose boolean
#' @return character(n)
#' @examples
#' require(magrittr)
#' x <- c("EM00", "EM01", "EM02")
#' guess_subgroup_values(x)
#'
#' x <- c("EM_R1", "EM_R2", "EM_R3")
#' guess_subgroup_values(x)
#'
#' x <- c("UT_10h_R1", "UT_10h_R2", "UT_10h_R3")
#' guess_subgroup_values(x)
#' guess_replicate_values(x)
#'
#' x <- c("EM00_STD.R1", "EM01_STD.R1", "EM01_EM00.R1")
#' guess_subgroup_values(x)
#' guess_replicate_values(x)
#' @export
guess_subgroup_values <- function(x, sep = guess_sep(x), verbose=TRUE){
    y <-  if (is.null(sep)){  x
            } else {         split_extract(x, seq_len(nfactors(x)-1), sep) }
    if (verbose)   message('\t\tGuess subgroup values: ', x[1], ' => ', y[1])
    y
}


#' @rdname guess_subgroup_values
#' @export
guess_replicate_values <- function(x, sep = guess_sep(x), verbose=TRUE){
    y  <-  if (is.null(sep)){   x
            } else {            split_extract(x, nfactors(x,sep), sep) }
    if (verbose)   message('\t\tGuess replicate values: ', x[1], ' => ', y[1])
    y
}


#=============================================================================
#
#                     create_subgroup_values
#                     create_replicate_values
#
#=============================================================================

#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' create_subgroup_values(object, subgroup_var=NULL,  verbose=TRUE)
#' create_replicate_values(object, subgroup_var=NULL, verbose=TRUE)
#' @noRd
create_subgroup_values <- function(object, subgroup_var, verbose){
    values <- svalues(object, subgroup_var)
    if (is.null(values) | all(is.na(values)) | all(values=="")) values <-
         guess_subgroup_values(object$sample_id, verbose=verbose)
    if (all(is.na(values)) | all(values=="")) values[] <- 'subgroup1'
    values
}


create_replicate_values <- function(object, subgroup_var, verbose){
    sampleid_values <- sdata(object)$sample_id
    if (is.null(subgroup_var)){
        replicate_values <- guess_replicate_values(
                                sampleid_values, verbose=FALSE)
    } else if (all(is.na(slevels(object, subgroup_var)))){
        replicate_values <- guess_replicate_values(
                                sampleid_values, verbose=FALSE)
    } else {
        subgroup_values <- sdata(object)[[subgroup_var]]
        replicate_values <- rep('', length(subgroup_values))
        for (i in seq_along(subgroup_values)){
            replicate_values[i] <-
                sampleid_values[i] %>%
                stri_replace_first_fixed(subgroup_values[i], '') %>%
                stri_replace_all_regex('^[._ ]', '') %>%
                stri_replace_all_regex('[._ ]$', '')
        }
    }
    if (verbose)   message('\t\tGuess replicate values from sampleids: ',
                            sampleid_values[1], ' => ', replicate_values[1])
    replicate_values
}


#=============================================================================
#
#               add_designvars
#                   file_exists
#                   get_default_designfile
#                       default_designfile
#
#=============================================================================

# Deals properly with NULL values
# file.exists does not!
file_exists <- function(file){
    if (is.null(file))      return(FALSE)
    if (file.exists(file))  return(TRUE)
                            return(FALSE)
}

default_designfile <- function(file, platform = NULL, quantity = NULL){

    # Initialize
    if (is.null(platform))  platform <- ''
    if (is.null(quantity))  quantity <- ''

    # No designfile for SOMASCAN and METABOLON
    if (platform %in% c('metabolon', 'somascan')) return(NULL)

    # Take basename file
    designfile <- tools::file_path_sans_ext(file)

    # Append quantity for MaxQuant files
    if (platform == 'maxquant'){
        designfile %<>% paste(make.names(quantity), sep = '.')
    }

    # Add .design.tx
    designfile %<>% paste0('.design.txt')
    designfile
}


#' @param object        SummarizedExperiment
#' @param subgroup_var  subgroup svar or NULL
#' @param designfile    design file path (to read/write) or NULL (don't write)
#' @param verbose       TRUE (default) or FALSE
#'@examples
#'# PROTEINGROUPS
#'    file <- download_data('billing19.proteingroups.txt')
#'    object <- read_proteingroups(file)
#'    get_default_designfile(object)
#'
#' # SOMASCAN
#'     inputfile <- download_data('atkin18.somascan.adat')
#'     default_designfile(inputfile, platform = 'somascan')
#'
#' # METABOLON
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     default_designfile(inputfile, platform = 'metabolon')
#'
#' # RNACOUNTS
#'@noRd
get_default_designfile <- function(object){
    file     <- metadata(object)$file
    platform <- metadata(object)$platform
    quantity <- metadata(object)$quantity

    default_designfile(file, platform, quantity)
}


#' @param object        SummarizedExperiment
#' @param subgroup_var  subgroup svar or NULL
#' @param designfile    design file path (to read/write) or NULL (don't write)
#' @param verbose       TRUE (default) or FALSE
#'@examples
#'# PROTEINGROUPS
#'    file <- download_data('billing19.proteingroups.txt')
#'    rm_subgroups <- c('BLANK_BM00', 'BLANK_STD', 'BM00_BM00', 'EM01_EM00',
#'                        'EM05_EM02', 'EM30_EM15')
#'    object <- read_proteingroups(file, rm_subgroups = rm_subgroups)
#'    add_designvars(object)
#'
#'    file <- download_data('billing16.proteingroups.txt')
#'    invert_subgroups <- c('EM_E', 'E_BM', 'EM_BM')
#'    object <- read_proteingroups(file, invert_subgroups = invert_subgroups)
#'    add_designvars(object)
#'
#' # SOMASCAN
#'     file <- download_data('atkin18.somascan.adat')
#'     read_somascan(file)
#'
#' # METABOLON
#'     file <- download_data('atkin18.metabolon.xlsx')
#'     read_metabolon(file)
#'
#' # RNACOUNTS
#'@noRd
add_designvars <- function(object, subgroup_var = NULL,
    designfile = get_default_designfile(object),
    verbose = TRUE
){
# Read
    dt <- if (file_exists(designfile)){
            if (verbose) message(
                '\t\tRead design (update if required!): ', designfile)
            fread(designfile)
        } else {
            if (verbose) message(
                '\t\tInfer design from sample_ids')
            data.table(
            sample_id = object$sample_id,
            subgroup  = create_subgroup_values( object, subgroup_var, verbose),
            replicate = create_replicate_values(object, subgroup_var, verbose))
        }
# Add to object
    object %<>% merge_sdata(dt)
    sdata(object) %<>% pull_columns(c('sample_id', 'subgroup', 'replicate'))
# Write
    if (!is.null(designfile)){
        if (!file.exists(designfile)){
            if (verbose) message('\t\tWrite design (update if required!): ',
                            designfile)
            fwrite(dt, designfile, sep = '\t', row.names = FALSE)
        }
    }
# Return
    object
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


#' Create design
#'
#'  Create design matrix  for statistical analysis
#' @param object      SummarizedExperiment
#' @param formula     formula with svars
#' @param ...         backward compatibility
#' @return design matrix
#' @examples
#' file <- download_data('billing16.rnacounts.txt')
#' object <- read_counts(file, plot=FALSE)
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
create_design <- function(
    object,
    formula = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup
){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    formula <- enquo(formula)
    formula <- eval_tidy(formula, sdata(object))
    assert_is_subset(all.vars(formula), svars(object))
    . <- NULL
# Ensure that subgroup vector is a factor to preserve order of levels
    object$subgroup %<>% factor() # required for '(Intercept)' -> factor1.level1
    for (var in setdiff(all.vars(formula), 'subgroup')){
        if (is.character(sdata(object)[[var]])){
            sdata(object)[[var]] %<>% factor()
        }
    }
# Create design matrix
    myDesign <- model.matrix(formula,  data = sdata(object))
# Rename: intercept -> factor1level1
    factors <- svars(object)[are_factor(sdata(object))]
    factor1 <- factors[1]
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


