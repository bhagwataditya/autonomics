#==========================================================
# GENEX
#==========================================================

#' Read genex file
#' @param file string: path to exiqon genex file
#' @return SummarizedExperiment
#' @noRd
read_genex <- function(file){
    assert_all_are_existing_files(file)
    dt <- extract_rectangle(file, sheet=1)
    read_omics(
        file,
        sheet = 1,
        fid_rows   = 1,                     fid_cols   = 2:(ncol(dt)-2),
        sid_rows   = 2:(nrow(dt)-3),        sid_cols   = 1,
        expr_rows  = 2:(nrow(dt)-3),        expr_cols  = 2:(ncol(dt)-2),
        fvar_rows  = (nrow(dt)-2):nrow(dt), fvar_cols  = 1,
        svar_rows  = 1,                     svar_cols  = (ncol(dt)-1):ncol(dt),
        fdata_rows = (nrow(dt)-2):nrow(dt), fdata_cols = 2:(ncol(dt)-2),
        sdata_rows = 2:(nrow(dt)-3),        sdata_cols = (ncol(dt)-1):ncol(dt),
        transpose  = TRUE,
        verbose    = TRUE)
}


