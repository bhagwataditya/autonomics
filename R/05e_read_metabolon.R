#======================================================
# METABOLON
#======================================================

find_origscale_sheet <- function(file){
    readxl::excel_sheets(file) %>%
    extract(stringi::stri_detect_fixed(., 'OrigScale')) %>%
    extract2(1)
}


#' Read metabolon
#' @param file         string: path to metabolon xlsx file
#' @param sheet        number/string: xls sheet number or name
#' @param fid_var      string: feature_id variable (ideally transcends dataset)
#' @param sid_var      string: sample_id variable
#' @param subgroup_var string: subgroup variable (human comprehensible)
#' @param fname_var    string: feature_name variable
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('hypoglycemia.metabolon.xlsx')
#' read_metabolon(file)
#' @export
read_metabolon <- function(file, sheet = find_origscale_sheet(file),
    fid_var      = '(COMP|COMP_ID)', sid_var = '(CLIENT_IDENTIFIER|Client ID)',
    subgroup_var = 'Group', fname_var    = 'BIOCHEMICAL'){
# Assert
    assert_all_are_existing_files(file)
# Initial read
    d_f <- read_excel(file, sheet, col_names = FALSE, .name_repair = 'minimal')

    fvar_rows <- which(!is.na(d_f %>% extract_dt_col(1))) %>% extract(1)
    svar_cols <- which(!is.na(d_f %>% extract_dt_row(1))) %>% extract(1)
    fvar_cols <- fdata_cols <- seq_len(svar_cols)
    svar_rows <- sdata_rows <- seq_len(fvar_rows)
    fvar_names <- d_f  %>%  extract_dt_row(fvar_rows)   %>%
                            extract(seq_len(svar_cols))
    svar_names <- d_f  %>%  extract_dt_col(svar_cols)   %>%
                            extract(seq_len(fvar_rows))
    fid_var <- fvar_names %>% extract(stri_detect_regex(., fid_var))
    sid_var <- svar_names %>% extract(stri_detect_regex(., sid_var))
    fid_rows  <- fdata_rows <- expr_rows <- (fvar_rows+1):nrow(d_f)
    sid_cols  <- sdata_cols <- expr_cols <- (svar_cols+1):ncol(d_f)

    fid_cols  <-  fvar_names %>% equals(fid_var) %>% which()
    sid_rows  <-  svar_names %>% is_in(sid_var) %>% which() %>% extract(1)
# Systematic read
    object <- read_omics(file,
                        sheet      = sheet,
                        fid_rows   = fid_rows,      fid_cols   = fid_cols,
                        sid_rows   = sid_rows,      sid_cols   = sid_cols,
                        expr_rows  = expr_rows,     expr_cols  = expr_cols,
                        fvar_rows  = fvar_rows,     fvar_cols  = fvar_cols,
                        svar_rows  = svar_rows,     svar_cols  = svar_cols,
                        fdata_rows = fdata_rows,    fdata_cols = fdata_cols,
                        sdata_rows = svar_rows,     sdata_cols = sdata_cols,
                        transpose  = FALSE,
                        verbose    = TRUE)
# Add sdata
    is_subgroup_col <- stri_detect_regex(svars(object), subgroup_var)
    subgroup_var <- if (any(is_subgroup_col)){  svars(object)[is_subgroup_col]
                    } else {                    sid_var }
    sdata(object) %<>% (function(y){
                            y$subgroup <- y[[subgroup_var]]
                            y %>% pull_columns(c('sample_id', 'subgroup'))})
# Add fdata
    assert_is_subset(fname_var, fvars(object))
    fdata(object) %<>% (function(y){ y$feature_name <- y[[fname_var]]
    y %>% pull_columns(c('feature_id', 'feature_name'))})
# return
    object
}

