
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
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% add_kegg_pathways()
#' @references http://www.kegg.jp/kegg/rest/keggapi.html
#' @noRd
add_kegg_pathways <- function(
    object, entry_var = 'KEGG', pathway_var = 'KEGGPATHWAYS'
){
# Add KEGG Pathways
    fdata(object)[[pathway_var]] <- fdata(object)[[entry_var]]     %>%
                                    split_extract_fixed(';', 1) %>%
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
#' object <- read_metabolon(file, plot = FALSE)
#' add_smiles(object[1:10, ])
#' @references https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial
#' @export
add_smiles <- function(object){
# Satify CHECK
    . <- NULL
# Assert
    assert_is_subset('PUBCHEM', fvars(object))
# Map to smiles
    PUBCHEMIDS <- fdata(object)$PUBCHEM %>% split_extract_fixed(';', 1)
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
    fidvar = 'BIOCHEMICAL', # '(COMP|COMP_ID)', 
    sidvar = '(CLIENT_IDENTIFIER|Client ID)',
    sfile = NULL, by.x = 'sample_id', by.y = NULL, subgroupvar = 'Group', 
    verbose  = TRUE
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
    fidvar <- fvar_names %>% extract(stri_detect_regex(., fidvar))
    sidvar <- svar_names %>% extract(stri_detect_regex(., sidvar))
    fid_rows  <- fdata_rows <- expr_rows <- (fvar_rows+1):nrow(d_f)
    sid_cols  <- sdata_cols <- expr_cols <- (svar_cols+1):ncol(d_f)
    fid_cols  <-  fvar_names %>% equals(fidvar) %>% which()
    sid_rows  <-  svar_names %>% is_in(sidvar) %>% which() %>% extract(1)
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
        transpose  = FALSE,         verbose    = verbose)
    assayNames(object)[1] <- paste0('metabolon')
# Update sdata/fdata                        Group   HMDB_ID -> HMDB_ID
    nsv <- length(svars((object)))
    nfv <- length(fvars((object)))
    svars(object)[nsv] %<>% stri_replace_first_regex('^([^ ]+)[ ]+([^ ]+)','$1')
    fvars(object)[nfv] %<>% stri_replace_first_regex('^([^ ]+)[ ]+([^ ]+)','$2')
    object %<>% merge_sample_file(sfile = sfile, by.x = by.x, by.y = by.y)
    object %<>% add_subgroup(subgroupvar, verbose = verbose)
# Return
    object
}


#' Read metabolon xlsxfile
#' @param file          metabolon xlsx file
#' @param sheet         excel sheet (number or string)
#' @param fidvar        featureid var
#' @param sidvar        samplid var
#' @param sfile         sample file
#' @param by.x          `file`  mergeby column
#' @param by.y          `sfile` mergeby column
#' @param subgroupvar   subgroup var
#' @param fnamevar      featurename fvar
#' @param kegg_pathways TRUE or FALSE: add kegg pathways?
#' @param smiles        TRUE or FALSE: add smiles?
#' @param impute        TRUE or FALSE: impute group-specific NA values?
#' @param plot          TRUE or FALSE
#' @param label         fvar
#' @param pca           TRUE or FALSE
#' @param pls           TRUE or FALSE
#' @param fit           model engine: 'limma', 'lm', 'lme(r)', 'wilcoxon' or NULL
#' @param formula       model formula
#' @param block         model blockvar: string or NULL
#' @param coefs         model coefficients of interest:    character vector or NULL
#' @param contrasts     coefficient contrasts of interest: character vector or NULL
#' @param palette       NULL or colorvector
#' @param verbose       TRUE or FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('atkin18.metabolon.xlsx')
#' read_metabolon(file, plot = TRUE, block = 'Subject')
#' @export
read_metabolon <- function(file, sheet = 'OrigScale',
    fidvar = 'BIOCHEMICAL', # '(COMP|COMP_ID)', 
    sidvar = '(CLIENT_IDENTIFIER|Client ID)',
    sfile = NULL, by.x = 'sample_id', by.y = NULL, subgroupvar = 'Group',
    fnamevar = 'BIOCHEMICAL', kegg_pathways = FALSE, smiles = FALSE,
    impute  = TRUE, plot = FALSE, pca = plot, pls = plot, label = 'feature_id',
    fit = if (plot) 'limma' else NULL, formula = ~ subgroup, block = NULL, 
    coefs = NULL, contrasts = NULL, palette = NULL, verbose = TRUE
){
# Read
    object <- .read_metabolon(
        file    = file,    sheet       = sheet, 
        fidvar  = fidvar,  sidvar      = sidvar, 
        sfile   = sfile,   by.x        = by.x, 
        by.y    = by.y,    subgroupvar = subgroupvar, 
        verbose = verbose)
# Prepare
    assert_is_subset(fnamevar, fvars(object))
    fdata(object)$feature_name <- fdata(object)[[fnamevar]]
    fdata(object) %<>% pull_columns(c('feature_id', 'feature_name'))
    object %<>% log2transform(verbose = verbose)
    if ({{impute}})     object %<>% impute(plot = FALSE)
    if (kegg_pathways)  object %<>% add_kegg_pathways('KEGG', 'KEGGPATHWAY')
    if (smiles)         object %<>% add_smiles('SMILES')
# Analyze
    object %<>% analyze(
        pca        = pca,           pls        = pls,
        fit        = fit,           formula    = formula,
        block      = block,         coefs      = coefs,
        contrasts  = contrasts,     plot       = plot,
        label      = label,
        palette    = palette,       verbose    = verbose)
# Return
    object
}

