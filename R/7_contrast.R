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


#==============================================================================
#
#                is_valid_contrast
#                select_valid_contrast
#                validify_contrasts
#
#==============================================================================


#' Is a valid contrast?
#' @param contrast contrast
#' @param design   design
#' @return logical
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' design <- create_design(object)
#' contrast <- default_contrasts(object)[1]
#' is_valid_contrast(contrast, design)
#' @noRd
is_valid_contrast <- function(contrast, design){
  contrast %<>% as.character()
  contrast %<>% stri_replace_all_fixed(' ', '')
  contrast %<>% stri_replace_all_regex('(/[0-9]+)', '') %>%
                stri_replace_all_fixed('(', '')         %>%
                stri_replace_all_fixed(')', '')
  terms <- strsplit(contrast, '[-+ ]+') %>% unlist() %>% unname()
  all(terms %in% colnames(design))
}


#' Select valid contrasts
#' @param  contrasts vector with contrast definitions
#' @param  design    design matrix
#' @param  verbose   whether or not to report
#' @return subset of contrasts
#' @noRd
select_valid_contrasts <- function(contrasts, design, verbose = TRUE){
  selector <- contrasts %>% vapply(is_valid_contrast, logical(1), design)
  if (verbose){
    cmessage('\t\tKeep %d valid contrasts (out of %d)',
            sum(selector), length(selector))
  }
  contrasts[selector]
}


#' Validy subgroups and contrasts
#' @param design    design matrix
#' @param contrasts contrasts vector
#' @return validified contrasts
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object  <- read_proteingroups(file)
#' design  <- create_design(object)
#' contrasts <- default_contrasts(object)
#' validify_contrasts(contrasts, design)
#' @noRd
validify_contrasts <- function(contrasts, design){

  # Assert valid inputs
   assert_is_matrix(design)
   assert_is_numeric(design)

  # Validify contrast names
  if (!has_names(contrasts)){
    names(contrasts) <- make.names(contrasts)
  }
  assert_has_no_duplicates(names(contrasts))

  # Select valid contrasts
  contrasts %<>% select_valid_contrasts(design)
  contrasts
}


#==============================================================================
#
#                       create_contrastmat
#
#==============================================================================


#' Create contrast matrix
#' @param object        SummarizedExperiment
#' @param contrasts  vector of contrast definitions
#' @param design     design matrix
#' @return contrast matrix
#' @examples
#' # PROTEINGROUPS
#' file <- download_data('billing16.proteingroups.txt')
#' invert_subgroups <- c('BM_EM', 'EM_E', 'BM_E')
#' object <- read_proteingroups(file, invert_subgroups = invert_subgroups)
#' create_contrastmat(object)
#'
#' # RNACOUNTS
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_counts(file)
#' create_contrastmat(object, c(EM0.8_0 = 'EM.8 - EM00'))
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' create_contrastmat(object)
#' @export
create_contrastmat <- function(
    object, contrasts = default_contrasts(object),
    design = create_design(object) # be explicit to disambiguate!
){
    if (!has_names(contrasts)) names(contrasts) <-make.names(contrasts)
    if (length(contrasts) == 0)  return(NULL)
    assert_has_no_duplicates(names(contrasts))
    makeContrasts(contrasts = contrasts, levels = design) %>%
    set_colnames(names(contrasts))
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
    for (i in 1:ncol(components)) components[[i]] %<>% factor(., levels=unique(.))
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
#' # Matrix
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot=FALSE)
#'     subgroup_matrix(object)
#' # Vector
#'     require(magrittr)
#'     file <-  download_data('billing19.proteingroups.txt')
#'     rm_subgroups <-  c('BLANK_BM00', 'BLANK_STD', 'BM00_BM00', 'EM01_EM00',
#'                        'EM05_EM02', 'EM30_EM15')
#'     object <- read_proteingroups(file, rm_subgroups=rm_subgroups, plot=FALSE)
#'     object$subgroup %<>% gsub('_STD', '', .)
#'     subgroups <- c('EM00','EM01','EM02','EM05','EM15','EM30','BM00')
#'     object$subgroup %<>% factor(subgroups)
#'     contrasts <- c(col_contrasts(object), row_contrasts(object))
#'     plot_contrastogram(object, contrasts=c(contrasts))
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
#               col_contrast_matrix
#               row_contrast_matrix
#
#==============================================================================


#' @export
#' @rdname col_contrasts
col_contrast_matrix <- function(object, symbol = ' - '){
    subgroupmat <- subgroup_matrix(object)
    contrastmat <- matrix(  sprintf('%s%s%s',
                                    subgroupmat[, -1],
                                    symbol,
                                    subgroupmat[, -ncol(subgroupmat)]),
                            nrow = nrow(subgroupmat),
                            ncol = ncol(subgroupmat)-1)
    rownames(contrastmat) <- rownames(subgroupmat)
    colnames(contrastmat) <- sprintf('%s - %s',
                            colnames(subgroupmat)[-1],
                            colnames(subgroupmat)[-ncol(subgroupmat)])
    contrastmat
}


#' Row/Col contrasts
#' @param object SummarizedExperiment
#' @param symbol contrazst symbol (default " - ")
#' @return  matrix (col_contrast_matrix) or vector (col_contrasts)
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' subgroup_matrix(object)
#' col_contrast_matrix(object)
#' col_contrasts(object)
#' row_contrast_matrix(object)
#' row_contrasts(object)
#' @export
col_contrasts <- function(object){
    col_contrmat   <- col_contrast_matrix(object)
    col_contrnames <- col_contrast_matrix(object, '__')
    c(structure(c(t(col_contrmat)), names = c(t(col_contrnames))))
}


#' @rdname col_contrasts
#' @export
row_contrast_matrix <- function(object, symbol = ' - '){
    subgroupmat <- subgroup_matrix(object)
    contrastmat <- matrix(  sprintf('%s%s%s',
                                  subgroupmat[-nrow(subgroupmat), ],
                                  symbol,
                                  subgroupmat[-1, ]),
                            nrow = nrow(subgroupmat)-1,
                            ncol = ncol(subgroupmat))
    colnames(contrastmat) <- colnames(subgroupmat)
    rownames(contrastmat) <- sprintf('%s - %s',
                            rownames(subgroupmat)[-nrow(subgroupmat)],
                            rownames(subgroupmat)[-1])
    contrastmat
}


#' @rdname col_contrasts
#' @export
row_contrasts <- function(object){
    row_contrmat   <- row_contrast_matrix(object)
    row_contrnames <- row_contrast_matrix(object, '__')
    c(structure(c(t(row_contrmat)), names = c(t(row_contrnames))))
}


#=============================================================================
#
#               aggregate_col_contrasts
#               aggregate_row_contrasts
#                   aggregate_contrasts
#
#==============================================================================


#' Aggregate contrasts
#'
#' Aggregate row contrasts across columns (or column contrasts across rows)
#'
#' @param contrastmat contrast matrix
#' @param dim 1 (aggregate across rows) or 2 (aggregate across columns)
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' aggregate_contrasts(row_contrast_matrix(object), 1)    # some conc across t
#' aggregate_contrasts(row_contrast_matrix(object), 2)    # concentrations at t
#' aggregate_contrasts(col_contrast_matrix(object), 2) # some time across conc
#' aggregate_contrasts(col_contrast_matrix(object), 1) # times at conc
#' @noRd
aggregate_contrasts <- function(contrastmat, dim){
    apply(contrastmat,
          dim,
          function(x){
              paste0(sprintf('(%s)/%d', x, length(x)), collapse = ' + ')})
}

aggregate_col_contrasts <- function(contrastmat){
    aggregate_contrasts(contrastmat, 2)
}

aggregate_row_contrasts <- function(conc_contrastmat){
    aggregate_contrasts(conc_contrastmat, 1)
}


#=============================================================================
#
#               default_contrasts
#
#==============================================================================


validify_contrast_names <- function(x){
   x %>% gsub(' ', '',  ., fixed = TRUE) %>%
         gsub('-', '_', ., fixed = TRUE) %>%
         make.names()
}


#' Default contrasts
#' @param object SummarizedExperiment
#' @return named character vector: contrast definitions
#' @examples
#' # Ratios
#'     require(magrittr)
#'     file <- download_data('billing16.proteingroups.txt')
#'     inv <- c('EM_BM', 'EM_E', 'E_BM')
#'     object <- read_proteingroups(file, invert_subgroups=inv, plot=FALSE)
#'     default_contrasts(object)
#'
#' # GLUTAMINASE
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot=FALSE)
#'     default_contrasts(object)
#' @export
default_contrasts <- function(object){

    # Ratios themselves for ratio data
    if (contains_ratios(object)){
        message('\tGenerate contrasts from ratios')
        return(
            subgroup_levels(object) %>% set_names(validify_contrast_names(.)))

    # Difference contrasts for abundance data
    } else {
        message('\tGenerate difference contrasts')
        sep <- guess_sep(object)
        contrasts <- c(col_contrasts(object), row_contrasts(object))
    }

    # Return
    contrasts
}


#==============================================================================
#
#                            add_limma
#
#==============================================================================


#' Add limma
#'
#' Run limma and add results
#'
#' Limma results can be easily accessed with limma(object).
#' This returns a list with components:
#' \itemize{
#'    \item {effect} matrix (ngene x ncontrast): effect sizes
#'    \item {rank}   matrix (ngene x ncontrast): effect ranks
#'    \item {t}      matrix (ngene x ncontrast): t    values (moderated t test)
#'    \item {se}     matrix (ngene x ncontrast): se   values (moderated t test)
#'    \item {p}      matrix (ngene x ncontrast): p    values (moderated t test)
#'    \item {fdr}    matrix (ngene x ncontrast): fdr  values (moderated t test)
#'    \item {bonf}   matrix (ngene x ncontrast): bonf values (moderated t test)
#'    \item {F}      vector (ngene)            : F    values (moderated F test)
#'    \item {F.p}    vector (ngene)            : p    values (moderated F test)
#' }
#' @param object        SummarizedExperiment
#' @param contrasts  contrast vector, preferably named
#'                   (automatically generated names are not always intuitive)
#' @param formula   formula to create design matrix (using svars)
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file)
#' add_limma(object)
#'
#' file <- download_data('billing19.proteingroups.txt')
#' rm_subgroups<-c('BLANK_BM00','BM00_BM00','EM01_EM00','EM05_EM02','EM30_EM15')
#' object <- read_proteingroups(file, rm_subgroups=rm_subgroups, plot=FALSE)
#' contrasts <- c(EM01_EM00 = 'EM01_STD - EM00_STD')
#' object %<>% add_limma(contrasts)
#' sum(limma(object)[,'EM01_EM00','bonf'] < 0.05, na.rm = TRUE)
#'
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_counts(file, plot=FALSE)
#' contrasts <- c(EM01_EM00 = 'EM01 - EM00')
#' object %<>% add_limma(contrasts)
#' sum(limma(object)[,'EM01_EM00','bonf'] < 0.05, na.rm = TRUE)
#' @export
add_limma <- function(object, contrasts = default_contrasts(object),
    formula = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup){
# Assert
    design <- create_design(object, formula=!!enquo(formula))
    assert_is_matrix(design)
    assert_is_numeric(design)
    assert_is_identical_to_true(unname(ncol(object)) == nrow(design))
    if (is.null(contrasts))    return(object)
    contrasts %<>% validify_contrasts(design)
    if (length(contrasts)==0) return(object)
# Set block and correlation if required
    cmessage('\t\tRun limma')
    block <- NULL; correlation <- NULL
    my_sdata <- sdata(object)
    if (has_complete_block_values(object)){
        cmessage("\t\tBlock on svar 'block'")
        block <- my_sdata$block
        correlation  <- duplicateCorrelation(exprs(object), design,
                                     block = block)[['consensus.correlation']]}
# Fit lm and compute contrasts
    fit <- suppressWarnings(lmFit(object = exprs(object), design = design,
                                  block = block, correlation = correlation,
                                  weights = weights(object)))
    object %<>% add_contrast_results(fit, contrasts)
    return(object)
}


add_contrast_results <- function(object, fit, contrasts){
    contrastmat <- create_contrastmat(object, contrasts, design)
    metadata(object)$contrastmat <- contrastmat
    fit %<>% contrasts.fit(contrasts = contrastmat)
    limma_quantities <- if (all(fit$df.residual==0)){ c('effect', 'rank')
                        } else { c('effect','rank','t','se','p','fdr','bonf')}
    limma(object) <- array( dim=c(nrow(fit),ncol(fit),length(limma_quantities)),
                            dimnames = list(feature  = rownames(fit),
                                            contrast = colnames(fit),
                                            quantity = limma_quantities))
    limma(object)[,,'effect'] <- fit$coefficients
    limma(object)[,,'rank'  ] <- apply(-abs(fit$coefficients), 2, rank)
# Perform moderated t test
    if (!all(fit$df.residual==0)){
        fit %<>% eBayes()
        pp <- fit$p.value
        limma(object)[,,'t' ] <- fit$t
        limma(object)[,,'se'] <- sqrt(fit$s2.post) * fit$stdev.unscaled
        limma(object)[,,'p' ] <- pp
        limma(object)[,,'rank'] <- apply(pp, 2, rank)
        limma(object)[,,'fdr' ] <- apply(pp, 2, p.adjust, 'fdr')
        limma(object)[,,'bonf'] <- apply(pp, 2, p.adjust, 'bonferroni')
        fdata(object)$F   <- fit$F
        fdata(object)$F.p <- fit$F.p
    }
    return(object)
}


#==============================================================================
#
#                    limma & limma<-
#                    extract_limma_dt
#
#==============================================================================


#' @title Get/set limma results
#' @description Get/Set limma results
#' @param object SummarizedExperiment
#' @param value list
#' @return limma results (get) or updated object (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(file, invert_subgroups=inv, plot=FALSE)
#' dim(limma(object))
#' dim(limma(object[1:5, ]))
#' @export
setGeneric("limma", function(object)   standardGeneric("limma") )


#' @rdname limma
setMethod("limma", signature("SummarizedExperiment"),
function(object){
    limma_array <- S4Vectors::metadata(object)$limma
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
#' object <- read_proteingroups(file, invert_subgroups=inv, plot=FALSE)
#' extract_limma_quantity(object)
#' @noRd
extract_limma_quantity <- function(object, quantity='p'){
    fvars0 <- c('feature_id','feature_name','imputed')
    fvars0 %<>% intersect(fvars(object))
    dt <- data.table(fdata(object)[, fvars0, drop=FALSE])
    dt %<>% cbind(adrop(limma(object)[, , quantity, drop=FALSE], drop=3))
    data.table::melt.data.table(
      dt, id.vars = fvars0, variable.name = 'contrast', value.name = quantity)
}

merge_limma_quantities <- function(x, y){
   names0 <- c('feature_id','feature_name','imputed', 'contrast')
   names0 %<>% intersect(names(x)) %>% intersect(names(y))
   merge(x, y, by = names0)
}


#' Extract limma datatable
#' @param object SummarizedExperiment
#' @return melted data.table
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file)
#' extract_limma_dt(object)
#'
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(file, invert_subgroups=inv, plot=FALSE)
#' extract_limma_dt(object)
#' @noRd
extract_limma_dt <- function(object){
    Reduce( merge_limma_quantities,
            mapply( extract_limma_quantity,
                    quantity = c('effect', 'p', 'fdr', 'bonf'),
                    MoreArgs = list(object=object),
                    SIMPLIFY = FALSE ))
}


#=============================================================================
#
#             plot_contrastogram
#                 compute_connections
#                     default_color_values2
#                         default_color_values
#                             default_color_var
#
#=============================================================================


#' default color_var
#' @param object SummarizedExperiment
#' @return default value of color_var
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' rm_subgroups <- c('EM_E', 'BM_E', 'EM_BM')
#' object <- read_proteingroups(file, rm_subgroups = rm_subgroups)
#' default_color_var(object)
#'
#' file <- download_data('billing16.somascan.adat')
#' object <- read_somascan(file)
#' default_color_var(object)
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' default_color_var(object)
#' @export
default_color_var <- function(object){
   if (     'block'    %in% svars(object))  'block'
   else if ('subgroup' %in% svars(object))  'subgroup'
   else                                      NULL
}


#' Default color values
#' @param object     SummarizedExperiment
#' @param color_var  string: svar mapped to color
#' @param show       logical
#' @return default color values vector
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' rm_subgroups <- c('EM_E', 'BM_E', 'EM_BM')
#' object <- read_proteingroups(file, rm_subgroups = rm_subgroups, plot=FALSE)
#' default_color_values(object, show = TRUE)
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' default_color_values(object, show = TRUE)
#' @export
default_color_values <- function(
   object,
   color_var = default_color_var(object),
   show = FALSE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(color_var, svars(object))
# sep
    sep              <- guess_sep(object, color_var)
    color_var_levels <- slevels(object, color_var)
# Make colors
    color_values <- make_colors(color_var_levels, sep, show)
    color_values[color_var_levels]
}


default_color_values2 <- function(object){
    subgrouplevels <- subgroup_levels(object)
    subgroupmatrix <- subgroup_matrix(object)
    default_color_values(object, 'subgroup')[c(t(subgroupmatrix))]
}

true_names <- function(x) names(x[x])


compute_connections <- function(
    object, subgroup_colors = default_color_values2(object)
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
    sizes <- colors <- matrix(0,
        nrow = length(subgrouplevels), ncol = length(subgrouplevels),
        dimnames = list(subgrouplevels, subgrouplevels))
    labels <- matrix("0", nrow = nrow(sizes), ncol = ncol(sizes),
                     dimnames = dimnames(sizes))
# Add contrast numbers
    contrastmat <- metadata(object)$contrastmat
    for (contrastname in colnames(contrastmat)){
        contrastvector <- contrastmat[, contrastname]
        to   <- true_names(contrastvector>0)
        from <- if (any(contrastvector<0)) true_names(contrastvector<0) else to
        ns <- nsignif[[contrastname]]
        nu <- nup[[contrastname]]
        nd <- ndown[[contrastname]]
        sizes[ to, from] <- nu#ns
        sizes[ from, to] <- nd#ns
        colors[to, from] <- subgroup_colors[[to]]
        colors[from, to] <- subgroup_colors[[from]]
        labels[to, from] <- if (nu>0) nu else 0#paste0(nu,  " %up% phantom(.)") else "phantom(.)"
        labels[from, to] <- if (nd>0) nd else 0#paste0(nd," %down% phantom(.)") else "phantom(.)"
    }
# Return
    #labels[colors==0] <- "0"
    list(sizes = sizes, colors = colors, labels = labels)
}


#' Plot contrastogram
#' @param object SummarizedExperiment
#' @param subgroup_colors named color vector (names = subgroups)
#' @examples
#' # subgroup matrix
#'    file <- download_data('halama18.metabolon.xlsx')
#'    object <- read_metabolon(file, plot=FALSE)
#'    plot_contrastogram(object)
#'
#' # subgroup vector
#'     require(magrittr)
#'     file <-  download_data('billing19.proteingroups.txt')
#'     rm_subgroups <-  c('BLANK_BM00', 'BLANK_STD', 'BM00_BM00', 'EM01_EM00',
#'                        'EM05_EM02', 'EM30_EM15')
#'     object <- read_proteingroups(file, rm_subgroups=rm_subgroups, plot=FALSE)
#'     object$subgroup %<>% gsub('_STD', '', .)
#'     object$subgroup %<>% factor(c('EM00','EM01','EM02','EM05','EM15','EM30','BM00'))
#'     contrasts <- c(col_contrasts(object), row_contrasts(object))
#'     plot_contrastogram(object, contrasts=c(contrasts))
#'
#' # subgroup scalar
#'    contrasts <- c(col_contrasts(object), row_contrasts(object))
#'    plot_contrastogram(object, contrasts=contrasts[1])
#'
#' # Ratios: self-contrasts
#'    file <- download_data('billing16.proteingroups.txt')
#'    invert <- c('EM_E', 'BM_E', 'BM_EM')
#'    object <- read_proteingroups(file, invert_subgroups=invert, plot=FALSE)
#'    plot_contrastogram(object)
#' @export
plot_contrastogram <- function(
    object, subgroup_colors = default_color_values2(object)
){
# Initialize
    V2 <- N <- NULL

# Perform limma
    #object %<>% add_limma(contrasts = contrasts)
    contrastogram_matrices <- compute_connections(
            object, subgroup_colors = subgroup_colors)
    sizes  <- contrastogram_matrices$sizes
    colors <- contrastogram_matrices$colors
    labels <- contrastogram_matrices$labels
    widths <- scales::rescale(sizes, c(0.01,30))

# Plot diagram
    dt <- split_subgroup_levels(object)
    nrow <- dt[, data.table::uniqueN(V2)]
    nperrow <- dt[, .N, by = 'V1'][, N]
    if (all(nperrow==1)) nperrow %<>% length()
    #dir.create('~/importomicscache/contrastogram')
    #pdf('~/importomicscache/contrastogram/directed_contrastogram.pdf',
    #width = 9, height = 9)
    labels %<>% as.data.frame()
    diagram::plotmat(labels, nperrow, relsize = 1, box.size = 0.05,
    #diagram::plotmat(sizes, nperrow, relsize = 1, box.size = 0.05,
        name = rownames(sizes), box.col = subgroup_colors,
        box.type = 'square', arr.lwd = widths, shadow.size =0, # sqrt(sizes)
        arr.lcol = colors, arr.col = colors)
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
#' object <- read_proteingroups(file, invert_subgroups=invert_subgroups)
#' effect <-      limma(object)[,1,'effect']
#' fdr    <-      limma(object)[,1,'fdr']
#' mlp    <- -log(limma(object)[,1,'p'])
#' ntop   <- 3
#' table(top_up(effect, fdr, mlp, ntop))
#' table(top_down(effect, fdr, mlp, ntop))
#' @noRd
top_up <- function(effect, fdr, mlp, ntop){
    fdr_ok   <- fdr  < 0.05
    coef_ok  <- effect >  1
    coef_top <- if(any(fdr_ok)){  effect > nmax(effect[fdr_ok], ntop+1)
                } else {          rep(FALSE, length(effect)) }
    mlp_top  <- if (any(coef_ok)){ mlp > nmax(mlp[coef_ok], ntop+1)
                } else {           rep(FALSE, length(effect)) }
    fdr_ok & coef_ok & (coef_top | mlp_top)
}


#' Create volcano datatable
#' @param object SummarizedExperiment
#' @param ntop   number: how many top features (either FC wise or p wise) to be annotated
#' @return data.table
#' @examples
#' file <- download_data("billing16.proteingroups.txt")
#' invert_subgroups <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(file, invert_subgroups=invert_subgroups)
#' print(make_volcano_dt(object))
#' @export
make_volcano_dt <- function(object, ntop = 3){

    effect <- p <- mlp <- topdown <- topup <- significance <- fdr <- NULL

    dt <- extract_limma_dt(object)
    dt %<>% extract(!is.na(effect) & !is.na(p))
    dt[, mlp  := -log10(p)]

    # Prepare volcano datatable
    # Note: Using effect <= 0 (rather than effect <0) is required.
    #       Otherwise (the very few) features with effect=0 will have no effect for 'significance'
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
#' @param object          SummarizedExperiment
#' @param label           fvar for labeling top features
#' @param contrastnames   character vector: contrasts for which to plot volcano
#' @param nrow            number: n rows in faceted plot
#' @param ntop            number: n top features to be annotated
#' @return ggplot object
#' @examples
#' # proteingroup group ratios
#'     file <- download_data("billing16.proteingroups.txt")
#'     inv <- c('EM_E', 'BM_E', 'BM_EM')
#'     object <- read_proteingroups(file, invert_subgroups=inv, plot=FALSE)
#'     plot_volcano(object)
#'
#' # proteingroup LFQ intensities
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     object <- read_proteingroups(file, plot=FALSE)
#'     plot_volcano(object)
#'
#' # metabolon intensities: complex design
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot=FALSE)
#'     plot_volcano(object)
#'     colcontrasts <- names(col_contrasts(object))
#'     rowcontrasts <- names(row_contrasts(object))
#'     plot_volcano(object, contrastnames=colcontrasts, nrow=4, ntop=1)
#'     plot_volcano(object, contrastnames=rowcontrasts, nrow=3, ntop=1)
#'
#' # proteingroup internalstandard ratios
#'     file <-  download_data('billing19.proteingroups.txt')
#'     rm_subgroups <-  c('BLANK_BM00', 'BLANK_STD', 'BM00_BM00', 'EM01_EM00',
#'                        'EM05_EM02', 'EM30_EM15')
#'     object <- read_proteingroups(file, rm_subgroups=rm_subgroups, plot=FALSE)
#'     require(magrittr)
#'     object$subgroup %<>% gsub('_STD', '', .)
#'     subgroups <- c('EM00','EM01','EM02','EM05','EM15','EM30','BM00')
#'     object$subgroup %<>% factor(subgroups)
#'     contrasts <- c(col_contrasts(object), row_contrasts(object))
#'     object %<>% add_limma(contrasts = contrasts)
#'     plot_volcano(object)
#' @export
plot_volcano <- function(
    object, contrastnames = colnames(limma(object)),
    label = feature_name,
    nrow  = 1, ntop = 3
){
# Assert
    assert_is_all_of(object, "SummarizedExperiment")
    assert_is_character(contrastnames)
    assert_is_subset(contrastnames, colnames(limma(object)))
    topup <- topdown <- effect <- mlp <- NULL
    label <- enquo(label)
# Prepare
    limma(object) %<>% extract(, contrastnames, , drop=FALSE)
    plotdt <- make_volcano_dt(object, ntop = ntop)
    txtdt  <- copy(plotdt)[topup==TRUE | topdown==TRUE]
    colorvalues <-c(hcl(h=  0, l=c(20, 70, 100), c=100),
                    hcl(h=120, l=c(100, 70, 20), c=100))
    names(colorvalues) <- levels(plotdt$significance)
# Plot
    imputed <- NULL # fallback when plotdt misses "imputed"
    significance <- NULL
    p <- ggplot(plotdt) +
        facet_wrap(~ contrast, nrow = nrow, scales = 'fixed') +
        geom_point( aes(x=effect, y=mlp, color = significance, shape=imputed),
                    na.rm = TRUE)
    if (!quo_is_null(label)){
        p <- p + geom_text_repel(
                        data = txtdt,
                        aes(x=effect, y=mlp, label=!!label, color=significance),
                        #hjust = 'outward',
                        na.rm = TRUE,
                        show.legend = FALSE)}#,
                        #direction = 'x'
    p + theme_bw() +
        scale_color_manual(values = colorvalues, name = NULL) +
        xlab(expression(log[2](FC))) +
        ylab(expression(-log[10](p))) #+
        #guides(color=FALSE)
}
