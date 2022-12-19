#=============================================================================
#
#                    filter_sample_type
#                    filter_sample_quality
#                    filter_feature_type
#                    filter_feature_quality
#                    rm_na_columns
#                    rm_single_value_cols
#
#=============================================================================

filter_sample_type <- function(object, sample_type, verbose){
    if ('SampleType' %in% svars(object)){ # missing in older versions
        SampleType <- NULL
        object %<>% filter_samples(
            SampleType %in% !!enquo(sample_type), verbose = verbose)
    }
    object
}


filter_sample_quality <- function(object, sample_quality, verbose){
    if ('RowCheck'   %in% svars(object)){ # sample quality
        RowCheck <- NULL
        object %<>% filter_samples(
            RowCheck %in% !!enquo(sample_quality), verbose = verbose)
    }
    object
}

filter_feature_type <- function(object, feature_type, verbose){
    if ('Type'       %in% fvars(object)){ # feature type
        Type <- NULL
        object %<>% filter_features(
            Type %in% !!enquo(feature_type), verbose = verbose)
    }
    object

}

filter_feature_quality <- function(object, feature_quality, verbose){
    if ('ColCheck'   %in% fvars(object)){ # feature quality
        ColCheck <- NULL
        object %<>% filter_features(
            ColCheck %in% !!enquo(feature_quality), verbose = verbose)
    }
    object
}

#' Rm columns with only nas
#' @param df dataframe
#' @return dataframe with re-ordered columns
#' @examples
#' df <- data.frame(
#'    symbol    = c('A1BG', 'A2M'),
#'    id        = c('1',    '2'),
#'    name      = c('alpha-1-B glycoprotein', 'alpha-2-macroglobulin'),
#'    relevance = c(NA_character_, NA_character_),
#'    version   = c('NA', 'NA'),
#'    type      = c('proteincoding', 'proteincoding'))
#' rm_na_columns(df)
#' @noRd
rm_na_columns <- function(df){
    Filter(function(x) !all(is.na(x)|x=='NA'), df) #
}


#'Rm single value columns
#'@param df dataframe
#'@return dataframe with informative columns
#'@examples
#' df <- data.frame(
#'    symbol    = c('A1BG', 'A2M'),
#'    id        = c('1',    '2'),
#'    name      = c('alpha-1-B glycoprotein', 'alpha-2-macroglobulin'),
#'    relevance = c(NA_character_, NA_character_),
#'    type      = c('proteincoding', 'proteincoding'))
#' rm_single_value_columns(df)
#' @noRd
rm_single_value_columns <- function(df){
    Filter(function(x) length(unique(x))>1, df)
}


#==============================================================================
#
#                         .read_somascan
#                          read_somascan
#
#=============================================================================


#' @rdname read_somascan
#' @export
.read_somascan <- function(
    file, fidvar = 'SeqId', sidvar = 'SampleId', 
    sfile = NULL, by.x = NULL, by.y = NULL, subgroupvar = 'SampleGroup', 
    verbose = TRUE
){
# Assert
    assert_all_are_existing_files(file)
    assert_is_a_string(fidvar)
    assert_is_a_string(sidvar)
    . <- NULL
# Understand file structure
    content <- readLines(file)
    n_row <- length(content)
    n_col <- max(count.fields(file, quote='', sep='\t'))
    f_row <- 1 + which(stri_detect_fixed(content, '^TABLE_BEGIN'))
    s_row <- which(stri_detect_regex(content[f_row:length(content)],
                                    '^\t+', negate=TRUE))[1] + f_row-1
    f_col <- content %>% extract(f_row) %>%
        stri_extract_first_regex('^\\t+') %>% stri_count_fixed('\t') %>% add(1)
    fid_cols <- (1+f_col):n_col
    fid_rows <- content[f_row:(s_row-1)] %>% stri_extract_first_words() %>%
        equals(fidvar) %>% which() %>% add(f_row-1)
    sid_rows <- (1+s_row):n_row
    sid_cols <- content %>% extract(s_row) %>% stri_extract_all_words() %>%
                unlist() %>% equals(sidvar) %>% which()
# Read
    object <- read_rectangles(file,
        fid_rows   =  fid_rows,         fid_cols   =  fid_cols,
        sid_rows   =  sid_rows,         sid_cols   =  sid_cols,
        expr_rows  = (s_row+1):n_row,   expr_cols  = (f_col+1):n_col,
        fvar_rows  =  f_row:(s_row-1),  fvar_cols  =  f_col,
        fdata_rows =  f_row:(s_row-1),  fdata_cols = (f_col+1):n_col,
        svar_rows  =  s_row,            svar_cols  =  seq_len(f_col-1),
        sdata_rows = (s_row+1):n_row,   sdata_cols =  seq_len(f_col-1),
        transpose  =  TRUE,             verbose    =  verbose)
# Add metadata/subgroup
    assayNames(object) <- 'somascan'
    object %<>% merge_sample_file(sfile = sfile, by.x = by.x, by.y = by.y)
    if (is.null(subgroupvar))  subgroupvar <- 'SampleGroup'
    object %<>% add_subgroup(subgroupvar)
    object
}


#' Read somascan adatfile
#' @param file                  somascan (adat) file
#' @param fidvar                featureid var
#' @param sidvar                sampleid  var
#' @param sfile                 sample file
#' @param by.x                 `file`  mergeby column
#' @param by.y                 `sfile` mergeby column
#' @param subgroupvar           subgroup var: string
#' @param fname_var             featurename var: string
#' @param sample_type           subset of c('Sample','QC','Buffer','Calibrator')
#' @param feature_type          subset of c('Protein', 'Hybridization Control Elution','Rat Protein')
#' @param sample_quality        subset of c('PASS', 'FLAG', 'FAIL')
#' @param feature_quality       subset of c('PASS', 'FLAG', 'FAIL')
#' @param rm_na_svars           TRUE/FALSE: rm NA svars?
#' @param rm_single_value_svars TRUE/FALSE: rm single value svars?
#' @param plot                  TRUE/FALSE
#' @param pca                   TRUE/FALSE: run pca?
#' @param fit                   model engine: 'limma', 'lm', 'lme(r)','wilcoxon' or NULL
#' @param formula               model formula
#' @param block                 model blockvar
#' @param coefs                 model coefficients          of interest: string vector or NULL
#' @param contrasts             model coefficient contrasts of interest: string vector or NULL
#' @param verbose               TRUE/FALSE: message?
#' @return Summarizedexperiment
#' @examples
#' file <- download_data('atkin18.somascan.adat')
#' read_somascan(file, pca = TRUE, fit = 'limma', block = 'Subject_ID')
#' @export
read_somascan <- function(file, fidvar = 'SeqId', sidvar = 'SampleId',
    sfile = NULL, by.x = NULL, by.y = NULL, subgroupvar = 'SampleGroup', 
    fname_var    = 'EntrezGeneSymbol',
    sample_type = 'Sample', feature_type = 'Protein',
    sample_quality  = c('FLAG', 'PASS'), feature_quality = c('FLAG', 'PASS'),
    rm_na_svars = FALSE, rm_single_value_svars = FALSE, plot = FALSE, 
    pca = plot, fit = if (plot) 'limma' else NULL, formula = NULL, block = NULL, 
    coefs = NULL, contrasts = NULL, verbose = TRUE
){
# Read
    object <- .read_somascan(
        file, fidvar = fidvar, sidvar = sidvar, sfile = sfile, 
        by.x = by.x, by.y = by.y, subgroupvar = subgroupvar, verbose = verbose)
    object$sample_id %<>% make.unique()
# Prepare
    assert_is_subset(fname_var, fvars(object))
    fdata(object)$feature_name <- fdata(object)[[fname_var]]
    fdata(object) %<>% pull_columns(c('feature_id', 'feature_name'))
    SampleType <- RowCheck <- Type <- ColCheck <- NULL
    assert_is_character(sample_type);
    assert_is_character(feature_type)
    assert_is_character(sample_quality)
    assert_is_character(feature_quality)
    assert_is_a_bool(rm_na_svars)
    assert_is_a_bool(rm_single_value_svars)
    object %<>% filter_sample_type(    {{sample_type}},     verbose)
    object %<>% filter_sample_quality( {{sample_quality}},  verbose)
    object %<>% filter_feature_type(   {{feature_type}},    verbose)
    object %<>% filter_feature_quality({{feature_quality}}, verbose)
    if (rm_na_svars)            sdata(object) %<>% rm_na_columns()
    if (rm_single_value_svars)  sdata(object) %<>% rm_single_value_columns()
    object %<>% log2transform(verbose = verbose)
# Analyze
    object %<>% analyze(
        pca         = pca,          fit          = fit, 
        subgroupvar = subgroupvar,  formula      = formula, 
        block       = block,        coefs        = coefs, 
        contrasts   = contrasts,    verbose      = verbose, 
        plot=plot)
# Return
    object
}


