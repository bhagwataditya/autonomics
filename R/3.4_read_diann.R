
#' proteingroup to isoforms
#' @param x proteingroups string vector
#' @param unique whether to remove duplicates
#' @return string vector
#' @examples 
#'  (x <- c('Q96JP5;Q96JP5-2', 'Q96JP5', 'Q96JP5-2;P86791'))
#'  pg_to_isoforms(x)
#'  pg_to_canonical(x)
#'  pg_to_isoforms( x, unique = FALSE)
#'  pg_to_canonical(x, unique = FALSE)
#' # .pg_to_isoforms(x[1])   # unexported dot functions
#' # .pg_to_canonical(x[1])  # operate on scalars
pg_to_canonical <- function(x, unique = TRUE){
    assert_is_character(x)
    unname(vapply(x, .pg_to_canonical, character(1), unique = unique))
}

.pg_to_canonical <- function(x, unique = TRUE){
    z <- unlist(stri_split_fixed(x, ';'))
    z %<>% split_extract_fixed('-', 1)
    if ({{unique}}) z %<>% unique()
    paste0(z, collapse = ';')
}

#' rdname pg_to_canonical
#' @export
pg_to_isoforms <- function(x, unique = TRUE){
    assert_is_character(x)
    unname(vapply(x, .pg_to_isoforms, character(1), unique = unique))
}

.pg_to_isoforms <- function(x, unique = TRUE){
    z <- unlist(stri_split_fixed(x, ';'))
    z %<>% split_extract_fixed('-', 2)
    z[z=='NA'] <- '0'
        # Sometimes the canonical isoform is NOT isoform-1 !
        # https://www.uniprot.org/uniprotkb/Q9H0P0/entry#sequences
    if ({{unique}}) z %<>% unique()
    paste0(z, collapse = ',')
}

#' Extract common substring
#' @param a first string
#' @param b second string
#' @return  string
#' @examples
#' # Sequences
#'   a <- "heart-specific Fatty acid binding protein"
#'   b <- "Fatty acid binding protein, isoform 3"
#'   extract_common_substr(a, b)
#'
#'   a <- "Small nuclear ribonucleoprotein-associated proteins B and B'"
#'   b <- "Small nuclear ribonucleoprotein-associated protein N"
#'   extract_common_substr(a, b)
#' @references https://stackoverflow.com/questions/28261825
#' @noRd
extract_common_substr <- function(a, b){

    tt <- drop(attr(adist(a, b, counts=TRUE), "trafos"))

    # Nothing in common
    if (!stri_detect_regex(tt, 'M+')) return('')

    # Something in common
    aa  <-  stri_sub(tt, stri_locate_all_regex(tt, '[DM]+')[[1]]) %>%
            paste0(collapse = '') %>% trimws()
            # paste is required because multiple substrings can be found
    #bb <- tt %>%
    # stri_sub(stri_locate_all_regex(tt, '[IM]+')[[1]]) %>%
    # paste0(collapse = '')

    stri_sub(a, stri_locate_all_regex(aa, 'M+')[[1]]) %>%
    paste0(collapse = '') %>%
    trimws()

    # different  = c(a %>%
    #    stri_sub(stri_locate_all_regex(aa, 'D+')[[1]]) %>%
    #    trimws(),
    # b %>%
    #    stri_sub(stri_locate_all_regex(bb, 'I+')[[1]])) %>%
    #    trimws())

}

#' Commonify strings
#' @param x character vector
#' @examples
#' # NO DIFFERENCES
#'    x <- c( 'Retrotransposon Gag-like protein 8B',
#'            'Retrotransposon Gag-like protein 8B')
#'    commonify_strings(x)
#' # TAILS DIFFER
#'    x <- c( 'Histone H2B type 1-K',
#'            'Histone H2B type 1-C/E/F/G/I')
#'    commonify_strings(x)
#'    x <- c("Small nuclear ribonucleoprotein-associated proteins B and B'",
#'           "Small nuclear ribonucleoprotein-associated protein N")
#'    commonify_strings(x)
#' # MORE COMPLEX DIFFERENCES
#'    x <- c( 'Fatty acid binding protein, isoform 3',
#'            'Fatty acid binding protein',
#'            'heart-specific Fatty acid binding protein',
#'            'heart-specific Fatty acid binding protein, isoform 3')
#'    commonify_strings(x)
#' # NOTHING IN COMMON
#'    x <- c('ABC1', 'DEF2')
#'    commonify_strings(x)
#' @noRd
commonify_strings <- function(x){
    . <- NULL
    common <- Reduce(extract_common_substr, x)
    alternate  <- if (common==''){  x
                } else {            stri_replace_first_fixed(x, common, '') %>%
                                    stri_replace_first_fixed(', ', '')      %>%
                                    trimws()
                }
    if (all(alternate == '')) return(common)

    alternate                          %>%
    unique()                           %>%
    (function(s){s[s==''] <- '.'; s})  %>%
    sort()                             %>%
    #magrittr::extract(.!='')          %>%
    paste0(collapse='|')             %>%
    paste0('(', ., ')')              %>%
    paste0(common, .)
}

commonify_collapsed_strings <- function(x, sep = ';'){
    commonify_strings(unlist(stri_split_fixed(x, sep)))
}

#' Forge PG descriptions
#' @param uniprot string vector (without duplicates)
#' @param protein NULL or string vector (with same length as uniprot)
#' @param fastadt data.table
#' @return string vector
#' @examples
#' # Without file 
#'     uniprot <- c('Q96JP5;Q96JP5-2', 'O75822', 'Q96AC1;Q96AC1-3;Q9BQL6')
#'     protein <- c('ZFP91', 'EIF3J', 'FERM2;FERM1')
#'     forge_pg_descriptions(uniprot, protein)
#' # With file
#'     file <- download_data('uniprot_hsa_20140515.fasta')
#'     fastadt <- read_fastahdrs(file)
#'     forge_pg_descriptions(uniprot, fastadt = fastadt)
#' @export
forge_pg_descriptions <- function(
    uniprot, Protein.Name = NULL, fastadt = NULL
){
    assert_is_character(uniprot)
    assert_has_no_duplicates(uniprot)
    
    if (is.null(fastadt)){
        assert_is_character(Protein.Name)
        pgdt <- data.table(uniprot = uniprot, Protein.Name = Protein.Name)
        pgdt$Protein.Name %<>% stri_replace_all_regex('_[A-Z]+', '') # rm '_HUMAN'
        pgdt[, isoform := pg_to_isoforms(uniprot) ]
        pgdt[, feature_id := Protein.Name]
        #pgdt[stri_detect_fixed(Protein.Name, ';'), feature_id := commonify_collapsed_strings(feature_id, ';'), by = 'feature_id']
        pgdt[isoform!='0', feature_id := paste0(feature_id, '(', isoform, ')') ]
        PG0 <- uniprot
        pgdt %<>% extract(PG0, on = 'uniprot')
        pgdt$feature_id
    } else {
        pgdt <- data.table(proId = uniprot, uniprot = uniprot)
        pgdt %<>% separate_rows(uniprot, sep = ';') %>% data.table()
        fastadt %<>% extract(, c('uniprot', 'protein', 'canonical', 'isoform'))
        pgdt %<>% merge(fastadt, by = 'uniprot', sort = FALSE)
        # pgdt %<>% drop_inferior(verbose = FALSE)
            # DIA-NN has a non-razor approach
            # P63151, Q66LE6, P63151;Q66LE6 are three different proteingroups !
            # curation is incompatible with this setup
        pgdt %<>% extract(order(proId, protein, isoform))
        pgdt %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = c('proId', 'canonical')) # Collapse PG isoforms
        pgdt[, feature_id := protein]
        pgdt[, isoform := stri_replace_all_fixed(isoform, ';', ',')]
        pgdt[isoform!='0', feature_id := paste0(feature_id, '(', isoform, ')')]
        pgdt[, isoform := NULL]
        pgdt %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = c('proId'))              # Collapse PG paralogs
        pgdt %<>% extract(uniprot, on = 'proId')
        pgdt$feature_id
    }
}

#' diann precursor quantity
#' @export
PRECURSOR_QUANTITY <- 'Precursor.Quantity'

#' @rdname read_diann
#' @export
.read_diann_precursors <- function(
    file, 
    precursor_quantity = PRECURSOR_QUANTITY, 
    fastadt  = NULL, 
    organism = NULL
){
# Assert
    assert_all_are_existing_files(file)
    assert_is_subset(precursor_quantity, c('Precursor.Quantity', 'Precursor.Normalised'))
# Read
    anncols <- c('Run', 'Genes', 'Protein.Names', 'Protein.Group', 'Precursor.Id', 
                 'Q.Value', 'PG.Q.Value', 'Global.PG.Q.Value')
    numcols <- c(precursor_quantity, 'PG.Quantity', 'PG.MaxLFQ')
    cols <- c(anncols, numcols)
    dt <- fread(file, select = cols)
    for (col in numcols){
        dt[[col]] %<>% stri_replace_first_fixed(',', '.') # 1977.16 but 1,35E+11
        dt[[col]] %<>% as.numeric()
    }
    setnames(dt, c('Genes', 'Protein.Names', 'Protein.Group'), 
                 c('gene',  'protein',       'uniprot'))
# Filter/Annotate
    # dt %<>% extract(Lib.Q.Value <= 0.05)
    # dt %<>% extract(Lib.PG.Q.Value <= 0.01)
    pgdt <- unique(dt[, .(uniprot, protein)])
    pgdt[, feature_name := forge_pg_descriptions(uniprot, protein, fastadt)]
    pgdt[, organism := split_extract_fixed(protein, '_', 2)]
    pgdt[, protein := NULL]
    dt %<>% .merge(pgdt, by = 'uniprot')
    dt %<>% extract(order(feature_name, Run, -get(precursor_quantity)))
# Summarize
    dt[, Precursor.No := seq_len(.N), by = c('uniprot', 'Run')]
    dt[, PG.Top1 :=     rev(sort(get(precursor_quantity)))[1],                  by = c('uniprot', 'Run')]
    dt[, PG.Top3 := sum(rev(sort(get(precursor_quantity)))[1:3], na.rm = TRUE), by = c('uniprot', 'Run')]
    dt[, PG.Sum  := sum(         get(precursor_quantity),        na.rm = TRUE), by = c('uniprot', 'Run')]
    cols <- c('Run', 'gene', 'feature_name', 'organism',  'protein', 'uniprot', 
              'Precursor.No', 'Precursor.Id', 
              'PG.Quantity', 'PG.Top1', 'PG.Top3', 'PG.Sum', 'PG.MaxLFQ', precursor_quantity, 
              'Q.Value', 'PG.Q.Value', 'Global.PG.Q.Value')
    dt %<>% pull_columns(cols)
# Filter
    if (!is.null(organism)){    idx <- dt$organism == organism
                                dt %<>% extract(idx)             }
# Return
    dt 
}

#' @rdname read_diann
#' @export
.read_diann_proteingroups <- function(
    file, precursor_quantity = PRECURSOR_QUANTITY, fastadt = NULL){
    dt <- .read_diann_precursors(file, precursor_quantity = precursor_quantity, 
                                 fastadt = fastadt)
    cols <- c('uniprot', 'feature_name', 'protein', 'First.Protein.Description', 
              'gene', 'Run', 'PG.Quantity', 
              'PG.Top1', 'PG.Top3', 'PG.Sum', 'PG.MaxLFQ', 'Precursor.No')
    dt %<>% extract(, cols, with = FALSE )
    dt[, N.Precursor := .N, by = c('uniprot', 'Run')]
    dt[, Precursor.No := NULL]
    dt %<>% unique()
    dt
}


#' Contaminants URL
#' @examples 
#' CONTAMINANTSURL
#' @export
CONTAMINANTSURL <- paste0(
    'http://lotus1.gwdg.de/mpg/mmbc/maxquant_input.nsf', 
    '/7994124a4298328fc125748d0048fee2/$FILE/',
    'contaminants.fasta')


#' Downloads contaminants
#' @param url contaminants file url (string)
#' @param overwrite TRUE or FALSE: overwrite existiung download?
#' @return filename (string)
#' @examples
#' download_contaminants('www.doesntexist.de', overwrite = TRUE) # msg
#' # download_contaminants()                  # download first time
#' # download_contaminants(overwrite = TRUE)  # download each  time
#' @export
download_contaminants <-  function(url = CONTAMINANTSURL, overwrite = FALSE){
    destdir <- file.path(R_user_dir("autonomics", "cache"), "maxquant")
    dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
    destfile <- paste0(destdir, '/contaminants.fasta')
    if (overwrite | !file.exists(destfile)){
        tryCatch(
            download.file(url, destfile, mode = 'wb'),
            error = function(e){
                 message('Automatic download failed: ', url, 
                         '\nDownload manually into\n', 
                         destfile)
                 destfile <<- NULL
            }
        )
    }
    return(destfile)
}

#' Read contaminants
#' @param file contaminant file
#' @return data.table
#' @examples
#' file <- download_contaminants()
#' dt <- read_contaminants(file)
#' @export
read_contaminants <-  function(file = download_contaminants()){
    if (!requireNamespace('Biostrings', quietly = TRUE)){
        stop("BiocManager::install('Biostrings'). Then re-run.") }
    assert_all_are_existing_files(file)
    assert_is_identical_to_true(substr(file, nchar(file)-4, nchar(file)) == 'fasta')
    y <- Biostrings::readAAStringSet(file)
    y %<>% names()
    y %<>% split_extract_fixed(' ', 1)
    y
}

#' Rm contaminants
#'
#' Rm contaminants from DIA-NN SumExp
#' @param object         SummarizedExperiment
#' @param contaminants   uniprots (character vector)
#' @param verbose        TRUE or FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('dilution.report.tsv')
#' object <- read_diann(file)
#' # object %<>% rm_diann_contaminants()
#' @export
rm_diann_contaminants <- function(
    object, contaminants = read_contaminants(), verbose = TRUE
){
    fdt0 <- fdt(object)
    fdt0[, uniprot := feature_id ]
    fdt0 %<>% separate_rows(uniprot, sep = ';') %>% data.table()
    fdt0[, contaminant := FALSE]
    fdt0[uniprot %in% contaminants, contaminant := TRUE]
    fdt0[, uniprot := NULL]
    fdt0 %<>% extract(, .(contaminant = any(contaminant)), by = 'feature_id')
    object %<>% merge_fdata(fdt0, by.x = 'feature_id', by.y = 'feature_id')
    object %<>% filter_features(!contaminant, verbose = verbose)
}

#' Read diann
#'
#' @param file               'report.tsv' file
#' @param fastadt             NULL or data.table
#' @param quantity           'PG.MaxLFQ', 'PG.Quantity', 'PG.Top1', 'PG.Top3', or 'PG.Sum'
#' @param precursor_quantity 'Precursor.Quantity' or 'Precursor.Normalized'
#' @param simplify_snames     TRUE/FALSE : simplify (drop common parts in) samplenames ?
#' @param contaminants        string vector: contaminant uniprots
#' @param impute              TRUE/FALSE : impute group-specific NA values?
#' @param plot                TRUE/FALSE
#' @param pca                 TRUE/FALSE : compute and plot pca ?
#' @param fit                 model engine: 'limma', 'lm', 'lme(r)', 'wilcoxon' or NULL
#' @param formula             model formula
#' @param block               model blockvar: string or NULL
#' @param coefs               model coefficients          of interest: string vector or NULL
#' @param contrasts           model coefficient contrasts of interest: string vector or NULL
#' @param feature_id          string: feature for summary plot
#' @param sample_id           string: sample  for summary plot
#' @param palette             color palette: named string vector
#' @param verbose             TRUE/FALSE
#' @return  data.table / SummarizedExperiment
#' @examples 
#' # Read & Analyze
#'    file <- download_data('dilution.report.tsv')
#'    object <- read_diann(file)
#'    snames(object) %<>% split_extract_fixed('_', 3)
#'    object$subgroup <- as.numeric(object$sample_id)
#'    analyze(object)
#'     
#' # Read data.table (lower-level)
#'     file <- download_data('dilution.report.tsv')
#'    (PR   <- .read_diann_precursors(file))       # precursor    dt
#'    (PG   <- .read_diann_proteingroups(file))    # proteingroup dt
#'    
#' # Compare Summarizations
#'      PG[PG.Quantity==PG.Top1] # matches      : 26063 (85%) proteingroups
#'      PG[PG.Quantity!=PG.Top1] # doesnt match :  4534 (15%) proteingroups
#'      run <- 'IPT_HeLa_1_DIAstd_Slot1-40_1_9997'
#'      PR[uniprot=='Q96JP5;Q96JP5-2' & Run == run, 1:6] #    match:    8884 ==   8884
#'      PR[uniprot=='P36578'          & Run == run, 1:6] # no match:  650887 != 407978
#'      PR[PG.Quantity != PG.Top1][feature_name == unique(feature_name)[1]][Run == unique(Run)[1]][1:2, 1:6]
#'      PR[PG.Quantity != PG.Top1][feature_name == unique(feature_name)[2]][Run == unique(Run)[1]][1:2, 1:6]
#'      PR[PG.Quantity != PG.Top1][feature_name == unique(feature_name)[3]][Run == unique(Run)[1]][1:3, 1:6]
#' @export
read_diann <- function(
    file, fastadt = NULL, 
    quantity = 'PG.MaxLFQ', 
    precursor_quantity = PRECURSOR_QUANTITY, 
    simplify_snames = TRUE,
    contaminants = character(0), 
    impute = FALSE, plot = FALSE, 
    pca = plot, fit = if (plot) 'limma' else NULL, formula = NULL, block = NULL,
    coefs = NULL, contrasts = NULL, feature_id = NULL, sample_id = NULL, 
    palette = make_subgroup_palette(object), verbose = TRUE
){
# Assert
    assert_all_are_existing_files(file)
    assert_are_identical(names(fread(file, select = 1:3, nrows = 0)), 
                         c('File.Name', 'Run', 'Protein.Group'))
    if (!is.null(fastadt))  assert_is_data.table(fastadt)
    assert_is_subset(quantity, 
         c('PG.MaxLFQ', 'PG.Quantity', 'PG.Top1', 'PG.Top3', 'PG.Sum'))
    if (!is.null(contaminants))  assert_is_character(contaminants)
# SumExp
    # values
        dt <- .read_diann_proteingroups(file, precursor_quantity = precursor_quantity, fastadt = fastadt)
        precursors <- data.table::dcast(dt, uniprot ~ Run, value.var = 'N.Precursor')
        precursors %<>% dt2mat()
        dcast_quantity <- function(quantity) data.table::dcast(dt, uniprot ~ Run, value.var = quantity)
        assays0 <- Map(dcast_quantity, quantity)           # Map retains names - lapply doesnt
        assays0 %<>% lapply(dt2mat)
        assays0 %<>% lapply(zero_to_na, verbose = verbose) # lapply allows for multiple arguments - Map doesnt
        assays0 %<>% lapply(nan_to_na,  verbose = verbose)
        assays0 %<>% lapply(log2)
        assays0 %<>% c(list(precursors = precursors))
        object <- SummarizedExperiment(assays0)
        analysis(object)$nfeatures <- nrow(object)
    # fdt
        fdt(object)$feature_id <- rownames(object)                                                             # feature_id
        fdt0 <- unique(dt[, .(feature_id = uniprot, organism = protein)])                                      # organism
        fdt0 %<>% uncollapse(organism, sep = ';')
        fdt0[, organism := split_extract_fixed(organism, '_', 2)]
        fdt0 <- fdt0[, .(organism = paste_unique(organism, collapse = ';')), by = 'feature_id']
        object %<>% merge_fdata(fdt0, by.x = 'feature_id', by.y = 'feature_id')
        fdt0 <- unique(dt[, .(feature_id = uniprot, feature_name, protein, gene, First.Protein.Description)])  # gene,protein,etc 
        object %<>% merge_fdata(fdt0, by.x = 'feature_id', by.y = 'feature_id')
        fdt(object) %<>% pull_columns(c('feature_id', 'gene', 'feature_name', 'protein', 'organism'))
    # sdt
        snames(object) <- colnames(object)
        if (simplify_snames)  snames(object) %<>% simplify_snames()
        object$subgroup  <- infer_subgroup( object$sample_id)
# Filter. Impute. Analyze
    if (length(contaminants)>0){
        object %<>% rm_diann_contaminants(contaminants, verbose = verbose)
    }
    object %<>% rm_missing_in_all_samples(verbose = verbose)
    object %<>% extract(order(rowVars(values(.), na.rm = TRUE)), )
    object %<>% filter_exprs_replicated_in_some_subgroup(verbose = verbose)
    if ({{impute}})   object %<>% impute()
    assays <- c(quantity, 'precursors')
    for (assay in assays)  object %<>% add_assay_means(assay)
    object %<>% analyze(
        pca          = pca,           fit       = fit,
        formula      = formula,       block     = block,      
        coefs        = coefs,         contrasts = contrasts,   
        verbose      = verbose,       plot      = plot,
        feature_id   = feature_id,    sample_id = sample_id,  
        palette      = palette)
    object
}

has_one_level <- function(x) length(unique(x))==1
x <- paste0('pi_exp_', c('wt_r1', 'wt_r2', 'kd_r1', 'kd_r2'))
simplify_snames <- function(x){
    sep <- guess_sep(x)
    if (sep == 'NOSEP')  return(x)
    x <- data.table(x = x)
    x <- x[, tstrsplit(x, split = sep)]
    idx <- !vapply(x, has_one_level, logical(1))
    x <- x[, idx, with = FALSE]
    x <- Reduce(function(a,b) paste(a, b, sep = sep), x)
    x
}

# sampleids <- paste0('pi_exp_', c('wt_r1', 'wt_r2', 'kd_r1', 'kd_r2'))
infer_subgroup <- function(sampleids){
    sep <- guess_sep(sampleids)
    if (sep == 'NOSEP')  return(rep('group0', length(sampleids)))
    n <- nfactors(sampleids)
    subgroups <- sampleids %>% split_extract_fixed(sep, 1:(n-1))
    if (any(duplicated(subgroups)))  return(subgroups)
    return('group0')
}

