
#=============================================================================
#
#                   add_kegg_pathways
#                       kegg_entry_to_pathways
#
#=============================================================================

#' Add KEGG Pathways
#'
#' @param object      SummarizedExperiment
#' @param entry_var   kegg entry fvar
#' @param pathway_var kegg pathway fvar
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' object %<>% add_kegg_pathways()
#' @references http://www.kegg.jp/kegg/rest/keggapi.html
#' @noRd
add_kegg_pathways <- function(
    object, entry_var = 'KEGG', pathway_var = 'KEGGPATHWAYS'
){
# Add KEGG Pathways
    fdata(object)[[pathway_var]] <- fdata(object)[[entry_var]]     %>%
                                    extract_first_from_collapsed() %>%
                                    kegg_entry_to_pathways()
# Report
    nkeggid  <- sum(!is.na(fdata(object)[[entry_var]]))
    npathway <- sum(!is.na(fdata(object)[[pathway_var]]))
    message('\t\tAdd KEGG Pathways: ', nrow(object), ' features ',
            '-> ', nkeggid,  ' map to KEGG IDS ',
            '-> ', npathway, ' map to KEGG Pathways')
# Return
    object
}


#' Map kegg entry values to kegg pathways
#'
#' @param x           charactervector, factorvector, or SummarizedExperiment
#' @param entry_var   kegg entry fvar
#' @param pathway_var kegg pathway fvar
#' @return character vector
#' @examples
#' if (requireNamespace('RCurl', quietly = TRUE)){
#'    x <- c("C07326", "C04742", "C18218", "C18218", NA_character_,
#'               NA_character_, "", "")
#'     kegg_entry_to_pathways(x)
#' }
#' @references http://www.kegg.jp/kegg/rest/keggapi.html
#' @noRd
kegg_entry_to_pathways <- function(x){
# Assert
    if (!requireNamespace('RCurl', quietly = TRUE)){
        stop("BiocManager::install('RCurl'). Then re-run.") }
# Satisfy check
    Entry <- Pathway <- . <- NULL
    x %<>% as.character()
# Map available values
    idx <- !is.na(x) & x!=''
    if (sum(idx)==0) return(character(length(x)))
    keggurl <- sprintf('http://rest.kegg.jp/link/pathway/%s',
                        paste0(unique(x[idx]), collapse = '+'))
    if (!RCurl::url.exists(keggurl)) return(character(length(x)))
    cachefile <- tempfile()
    download.file(keggurl, cachefile, quiet = TRUE)
    if (readLines(cachefile, n=1)=='') return(character(length(x)))
# Format and return
    fread(cachefile, header = FALSE, col.names = c("Entry", "Pathway"))      %>%
    extract(, Entry   := stri_replace_first_fixed(Entry, 'cpd:',  ''))       %>%
    extract(, Pathway := stri_replace_first_fixed(Pathway, 'path:', ''))     %>%
    extract(, list(Pathway = paste0(Pathway, collapse = ';')), by = 'Entry') %>%
    merge(  data.table(Entry = x), .,
            by = 'Entry', all.x = TRUE, sort = FALSE, ) %>%
    extract2('Pathway')
}


#===================================================================
#
#                    add_smiles
#                        pubchem_to_smiles
#
#===================================================================

#' Add smiles
#'
#' @param object character/factor vector with pubchem ids
#' @return character/factor vector
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' add_smiles(object[1:10, ])
#' @references https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial
#' @export
add_smiles <- function(object){
# Satify CHECK
    . <- NULL
# Assert
    assert_is_subset('PUBCHEM', fvars(object))
# Map to smiles
    PUBCHEMIDS <- fdata(object)$PUBCHEM %>% split_extract(1,';')
    SMILES <- rep(NA_character_, length(PUBCHEMIDS))
    idx <- !is.na(PUBCHEMIDS)
    SMILES[idx] <- PUBCHEMIDS[idx] %>%
            (function(object) split(object, ceiling(seq_along(object)/100))) %>%
            lapply(pubchem_to_smiles) %>%
            unlist() %>%
            unname()
    fdata(object)$SMILES <- SMILES
# Report
    npubchem <- sum(!is.na(fdata(object)$PUBCHEM))
    nsmiles  <- sum(!is.na(fdata(object)$SMILES))
    message('\t\tAdd SMILES: ', nrow(object), ' features -> ', 
            npubchem, ' map to PUBCHEM -> ', nsmiles, ' map to SMILES')
# Return
    object
}


#' Map a vector of pubchemids to (canonical) smiles
#'
#' @param x character/factor vector with pubchem ids
#' @return character/factor vector
#' @examples
#' x <- c(NA_character_, "10236635", "5283147", "91477", NA_character_)
#' pubchem_to_smiles(x)
#' @references https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial
#' @noRd
pubchem_to_smiles <- function(x){
# Satisfy CHECK
    . <- NULL
# Download pubchem smiles
    cachefile <- tempfile()
    resturl <- sprintf(
                    paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/',
                            'compound/cid/%s/property/CanonicalSMILES/CSV'),
                    paste0(unique(as.vector(na.exclude(x))), collapse = ','))
    download.file(resturl, cachefile, quiet = TRUE)
# Return
    data.table(CID = as.integer(as.character(x))) %>%
    merge(fread(cachefile), by = 'CID', all.x = TRUE, sort = FALSE) %>%
    extract2('CanonicalSMILES')
}


#==============================================================================
#
#                           .read_metabolon
#                            read_metabolon
#
#==============================================================================

#' @rdname read_metabolon
#' @export
.read_metabolon <- function(file, sheet = 'OrigScale',
    fidvar = '(COMP|COMP_ID)', sidvar = '(CLIENT_IDENTIFIER|Client ID)',
    sfile = NULL, by.x = 'sample_id', by.y = NULL, subgroupvar = 'Group'
){
# Assert
    assert_all_are_existing_files(file)
    . <- NULL
# Initial read
    sheet %<>% grep(excel_sheets(file), fixed = TRUE, value = TRUE)
    d_f <- read_excel(file, sheet, col_names = FALSE, .name_repair = 'minimal')
    fvar_rows <- which(!is.na(d_f %>% extract_dt_col(1))) %>% extract(1)
    svar_cols <- which(!is.na(d_f %>% extract_dt_row(1))) %>% extract(1)
    fvar_cols <- fdata_cols <- seq_len(svar_cols)
    svar_rows <- sdata_rows <- seq_len(fvar_rows)
    fvar_names <- extract_dt_row(d_f, fvar_rows) %>% extract(seq_len(svar_cols))
    svar_names <- extract_dt_col(d_f, svar_cols) %>%extract(seq_len(fvar_rows))
    fid_var <- fvar_names %>% extract(stri_detect_regex(., fid_var))
    sid_var <- svar_names %>% extract(stri_detect_regex(., sid_var))
    fid_rows  <- fdata_rows <- expr_rows <- (fvar_rows+1):nrow(d_f)
    sid_cols  <- sdata_cols <- expr_cols <- (svar_cols+1):ncol(d_f)
    fid_cols  <-  fvar_names %>% equals(fid_var) %>% which()
    sid_rows  <-  svar_names %>% is_in(sid_var) %>% which() %>% extract(1)
# Systematic read
    object <- read_rectangles(
                file,                       sheet      = sheet,
                fid_rows   = fid_rows,      fid_cols   = fid_cols,
                sid_rows   = sid_rows,      sid_cols   = sid_cols,
                expr_rows  = expr_rows,     expr_cols  = expr_cols,
                fvar_rows  = fvar_rows,     fvar_cols  = fvar_cols,
                svar_rows  = svar_rows,     svar_cols  = svar_cols,
                fdata_rows = fdata_rows,    fdata_cols = fdata_cols,
                sdata_rows = svar_rows,     sdata_cols = sdata_cols,
                transpose  = FALSE, verbose    = TRUE)
    assayNames(object)[1] <- paste0('metabolon')
# Update sdata/fdata                        Group   HMDB_ID -> HMDB_ID
    nsv <- length(svars((object)))
    nfv <- length(fvars((object)))
    svars(object)[nsv] %<>% stri_replace_first_regex('^([^ ]+)[ ]+([^ ]+)','$1')
    fvars(object)[nfv] %<>% stri_replace_first_regex('^([^ ]+)[ ]+([^ ]+)','$2')
    object %<>% merge_sfile(sfile = sfile, by.x = by, by.y = sfileby)
    if (is.null(subgroupvar)) subgroupvar <- 'Group'
    object %<>% add_subgroup(subgroupvar)
# Return
    object
}

#' Read metabolon
#' @param file            metabolon xlsx filepath
#' @param sheet           xls sheet number or name
#' @param fid_var         feature_id fvar
#' @param sid_var         sampleid svar
#' @param sfile           sample file
#' @param sfileby         sample file mergeby column
#' @param by              metabolon file mergeby column
#' @param subgroupvar     subgroup svar
#' @param fname_var       featurename fvar
#' @param impute          whether to impute
#' @param add_kegg_pathways  whether to add kegg pathways
#' @param add_smiles      whether to add smiles
#' @param pca             whether to pca
#' @param fit       fit model: NULL, 'limma', 'lm', 'lme', 'lmer', 'wilcoxon'
#' @param formula         designmat formula
#' @param block           block svar
#' @param coefs           NULL or character vector: model coefficients to test
#' @param contrastdefs    NULL or character vector: coefficient contrasts to test
#' @param verbose         whether to msg
#' @param plot            whether to plot
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' read_metabolon(file, pca = TRUE, fit = 'limma', block='SUB')
#' @export
read_metabolon <- function(file, sheet = 'OrigScale',
    fid_var = '(COMP|COMP_ID)', sid_var = '(CLIENT_IDENTIFIER|Client ID)',
    sfile = NULL, sfileby = NULL, by = NULL, subgroupvar = 'Group',
    fname_var = 'BIOCHEMICAL',
    impute  = FALSE, add_kegg_pathways = FALSE, add_smiles = FALSE,
    pca = FALSE, fit = NULL, formula = NULL, block = NULL, 
    coefs = NULL, contrastdefs = NULL, verbose = TRUE, plot = TRUE
){
# Read
    subgroup <- if (is.null(subgroupvar)) quo(NULL) else sym(subgroupvar)
    object <- .read_metabolon(
        file = file, sheet = sheet, fid_var = fid_var, sid_var = sid_var, 
        sfile = sfile, sfileby = sfileby, by = by, subgroupvar = subgroupvar)
# Prepare
    assert_is_subset(fname_var, fvars(object))
    fdata(object)$feature_name <- fdata(object)[[fname_var]]
    fdata(object) %<>% pull_columns(c('feature_id', 'feature_name'))
    object %<>% log2transform(verbose = TRUE)
    if (impute)             object %<>% impute_systematic_nondetects(
                                            subgroup = !!subgroup, plot = FALSE)
    if (add_kegg_pathways)  object %<>% add_kegg_pathways('KEGG', 'KEGGPATHWAY')
    if (add_smiles)         object %<>% add_smiles('SMILES', 'PUBCHEM')
# Analyze
    object %<>% analyze(pca = pca, fit = fit, subgroupvar = subgroupvar, 
                    formula = formula, block = block, 
                    coefs=coefs, contrastdefs = contrastdefs, 
                    verbose = verbose, plot = plot)
# Return
    object
}

