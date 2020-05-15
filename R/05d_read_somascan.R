#=======================================================
# SOMASCAN
#=======================================================

#' Read somascan
#'
#' Read data from somascan adat file
#'
#' @param file         string: path to *.adat file
#' @param fid_var      string: feature_id   variable
#' @param sid_var      string: sample_id    variable
#' @param subgroup_var string: subgroup     variable
#' @param fname_var    string: feature_name variable
#' @return Summarizedexperiment
#' @seealso prepare_somascan
#' @examples
#' file <- download_data('hypoglycemia.somascan.adat')
#' read_somascan(file)
#' @export
read_somascan <- function(file, fid_var = 'SeqId', sid_var = 'SampleId',
    subgroup_var = 'SampleGroup', fname_var    = 'EntrezGeneSymbol'
){
    # Assert
    assert_all_are_existing_files(file)
    assert_is_a_string(fid_var)
    assert_is_a_string(sid_var)

    # Understand file structure
    content <- readLines(file)
    n_row <- length(content)
    n_col <- max(count.fields(file, quote='', sep='\t'))
    f_row <- 1 + which(stri_detect_fixed(content, '^TABLE_BEGIN'))
    s_row <- content %>% extract(f_row:length(.)) %>%
        detect_index(function(y) stri_detect_regex(y, '^\t+', negate=TRUE)) %>%
        add(f_row-1)
    f_col <- content %>% extract(f_row) %>%
        stri_extract_first_regex('^\\t+') %>% stri_count_fixed('\t') %>% add(1)
    fid_cols <- (1+f_col):n_col
    fid_rows <- content %>% extract(f_row:(s_row-1)) %>%
        stri_extract_first_words() %>% equals(fid_var) %>% which() %>%
        add(f_row-1)
    sid_rows <- (1+s_row):n_row
    sid_cols <- content %>% extract(s_row) %>% stri_extract_all_words() %>%
                unlist() %>% equals(sid_var) %>% which()
    # Read
    object <- read_omics(
        file,
        fid_rows   =  fid_rows,         fid_cols   =  fid_cols,
        sid_rows   =  sid_rows,         sid_cols   =  sid_cols,
        expr_rows  = (s_row+1):n_row,   expr_cols  = (f_col+1):n_col,
        fvar_rows  =  f_row:(s_row-1),  fvar_cols  =  f_col,
        fdata_rows =  f_row:(s_row-1),  fdata_cols  = (f_col+1):n_col,
        svar_rows  =  s_row,            svar_cols  = seq_len(f_col-1),
        sdata_rows = (s_row+1):n_row,   sdata_cols = seq_len(f_col-1),
        transpose  = TRUE,
        verbose    = TRUE)

    # Add sdata/fdata
    sdata(object) %<>% (function(y){ y$subgroup <- y[[subgroup_var]]
    y %>% pull_columns(c('sample_id', 'subgroup'))})

    assert_is_subset(fname_var, fvars(object))
    fdata(object) %<>% (function(y){ y$feature_name <- y[[fname_var]]
    y %>% pull_columns(c('feature_id', 'feature_name'))})

    # Return
    object
}


