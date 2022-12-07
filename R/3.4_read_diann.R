
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


#' diann precursor quantity
#' @export
PRECURSOR_QUANTITY <- 'Precursor.Quantity'


cols <- function(file)  names(fread(file, nrow = 0))
col1 <- function(file)  cols(file)[1]
col2 <- function(file)  cols(file)[2]
col3 <- function(file)  cols(file)[3]


# x <- 'Q15149;Q15149-3;Q15149-4;Q15149-8'
# uniprot2isoforms(x)
uniprot2isoforms <- function(x){
    x %<>% stri_split_fixed(';') 
    x %<>% unlist() 
    x %<>% split_extract_fixed('-', 2)
    x[x=='NA'] <- '0'
    #x %<>% unique()  #doesnt work with the way diann organises its proteingroups
    x %<>% sort()
    paste0(x, collapse = ',')
}

#' @rdname read_diann
#' @export
.read_diann_precursors <- function(
    file, precursor_quantity = PRECURSOR_QUANTITY, verbose  = TRUE
){
# Assert
    assert_diann_file(file)
    assert_is_subset(precursor_quantity, c('Precursor.Quantity', 'Precursor.Normalised'))
# Read
    anncols <- c('Run', 'Genes', 'Protein.Names', 'Protein.Group', 'Precursor.Id',
                 'Q.Value', 'PG.Q.Value', 'Global.PG.Q.Value', 'Stripped.Sequence')
    numcols <- c(precursor_quantity, 'PG.Quantity', 'PG.MaxLFQ')
    cols <- c(anncols, numcols)
    dt <- fread(file, select = cols)                  # 1977.16 but 1,35E+11
    for (col in numcols){   dt[, (col) := stri_replace_first_fixed(get(col), ',', '.') ] 
                            dt[, (col) := as.numeric(get(col))  ] }
    setnames(dt, 'Genes',         'gene')
    setnames(dt, 'Protein.Names', 'protein')
    setnames(dt, 'Protein.Group', 'uniprot')
    setnames(dt, 'Precursor.Id',  'precursor')
    setnames(dt, 'Stripped.Sequence', 'sequence')
# Order precursors
    dt <- dt[, .SD[rev(order(get(precursor_quantity)))], by = c('uniprot', 'Run')]
    dt[, iprecursor := seq_len(.N),                      by = c('uniprot', 'Run')]
    dt[, nprecursor := length(unique(precursor)),        by = c('uniprot', 'Run')]
    dt[, npeptide   := length(unique(sequence)),         by = c('uniprot', 'Run')]
# Order proteingroups
    pgdt <- dt[, .(uniprot, Run, npeptide, nprecursor, PG.MaxLFQ)]
    pgdt %<>% unique()
    pgdt %<>% extract(, .(npeptide   = sum(npeptide), 
                          nprecursor = sum(nprecursor), 
                     log2.PG.MaxLFQ  = log2(sum(PG.MaxLFQ, na.rm = TRUE))), by = 'uniprot')
    pgdt %<>% extract(order(-npeptide, -nprecursor, -log2.PG.MaxLFQ))
    dt[, uniprot := factor(uniprot, pgdt$uniprot)]
    dt %<>% extract(order(uniprot, Run, iprecursor))
    dt[, uniprot := as.character(uniprot)]
# Intuify protein
    pgdt <- unique(dt[, .(uniprot, protein)])
    pgdt %<>% uncollapse(protein, sep = ';')                                     #     uncollapse
    pgdt[, organism := split_extract_fixed(protein, '_', 2)]                     #     drop organism
    pgdt[, protein  := split_extract_fixed(protein, '_', 1)]                     # 
    pgdt[, protein := commonify_strings(protein), by = c('uniprot', 'organism')] #     commonify  (within proteingroup/organism)
    pgdt %<>% recollapse(by = c('uniprot', 'organism'), sep = ';')               #     recollapse (within proteingroup/organism)
    pgdt[, protein := paste0(protein, '_', organism)]                            #     add organism
    pgdt %<>% recollapse(by = 'uniprot')                                         #     recollapse (within proteingroup)
# Add feature_id
    pgdt[, isoform := uniprot2isoforms(uniprot), by = 'uniprot']
    pgdt[, feature_id := paste0(protein, '-', isoform)]
    assert_has_no_duplicates(pgdt$feature_id)
    pgdt[, c('isoform') := NULL]
    # pgdt[, feature_name := forge_pg_descriptions(uniprot, protein, fastadt)]     #     add feature_name
    dt %<>% .merge(pgdt, by = 'uniprot')
    dt %<>% pull_columns(c('gene', 'protein', 'organism', 'feature_id', 'uniprot', 
                'Run', 'npeptide', 'nprecursor', 'iprecursor', 'precursor', 'sequence'))
# Summarize
    dt[, PG.Top1 :=     rev(sort(get(precursor_quantity)))[1],                  by = c('uniprot', 'Run')]
    dt[, PG.Top3 := sum(rev(sort(get(precursor_quantity)))[1:3], na.rm = TRUE), by = c('uniprot', 'Run')]
    dt[, PG.Sum  := sum(         get(precursor_quantity),        na.rm = TRUE), by = c('uniprot', 'Run')]
    dt 
}

#' @rdname read_diann
#' @export
.read_diann_proteingroups <- function(
    file, precursor_quantity = PRECURSOR_QUANTITY
){
    dt <- .read_diann_precursors(file, precursor_quantity = precursor_quantity)
    dt[, sequence := sequence[1], by = c('uniprot', 'Run')]
    cols <- c('gene', 'feature_id', 'protein', 'organism', 'uniprot', 'Run',
              'npeptide', 'nprecursor', 'sequence',
              'PG.Quantity', 'PG.Top1', 'PG.Top3', 'PG.Sum', 'PG.MaxLFQ', 
              'PG.Q.Value', 'Global.PG.Q.Value')
    dt %<>% extract(, cols, with = FALSE )
    dt %<>% unique()
    #dt[, .SD[.N>1], by = c('Run', 'feature_id')] # single row per run/protein - yes!
    dt
}



#' Read diann
#'
#' @param file               'report.tsv' file
#' @param precursor_quantity 'Precursor.Quantity' or 'Precursor.Normalized'
#' @param simplify_snames     TRUE/FALSE : simplify (drop common parts in) samplenames ?
#' @param contaminants        string vector: contaminant uniprots
#' @param impute              TRUE / FALSE : impute group-specific NA values?
#' @param plot                TRUE / FALSE
#' @param pca                 TRUE / FALSE : compute and plot pca ?
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
#' # Read
#'    file <- download_data('dilution.report.tsv')
#'    .read_diann_precursors(file)      #    precursors longdt
#'    .read_diann_proteingroups(file)   # proteingroups longdt
#'    fdt(read_diann(file))             # proteingroups sumexp
#' # Compare
#'     PR <- .read_diann_precursors(file)
#'     PG <- .read_diann_proteingroups(file)
#'     PG[PG.Quantity==PG.Top1] # matches      : 26063 (85%) proteingroups
#'     PG[PG.Quantity!=PG.Top1] # doesnt match :  4534 (15%) proteingroups
#'     run <- 'IPT_HeLa_1_DIAstd_Slot1-40_1_9997'
#'     PR[uniprot=='Q96JP5;Q96JP5-2' & Run == run, 1:6] #    match:    8884 ==   8884
#'     PR[uniprot=='P36578'          & Run == run, 1:6] # no match:  650887 != 407978
#'     PR[PG.Quantity != PG.Top1][feature_name == unique(feature_name)[1]][Run == unique(Run)[1]][1:2, 1:6]
#'     PR[PG.Quantity != PG.Top1][feature_name == unique(feature_name)[2]][Run == unique(Run)[1]][1:2, 1:6]
#'     PR[PG.Quantity != PG.Top1][feature_name == unique(feature_name)[3]][Run == unique(Run)[1]][1:3, 1:6]
#' @export
read_diann <- function(
    file, 
    precursor_quantity = PRECURSOR_QUANTITY, 
    simplify_snames = TRUE,
    contaminants = character(0), 
    impute = FALSE, plot = FALSE, 
    pca = plot, fit = if (plot) 'limma' else NULL, formula = NULL, block = NULL,
    coefs = NULL, contrasts = NULL, feature_id = NULL, sample_id = NULL, 
    palette = NULL, verbose = TRUE
){
# Assert
    if (!is.null(contaminants))  assert_is_character(contaminants)
# SumExp
    dt <- .read_diann_proteingroups(file, precursor_quantity = precursor_quantity)
    object <- SummarizedExperiment(list(
        log2.MaxLFQ      = dcast_diann(dt, 'PG.MaxLFQ',   fill = NA, log2 = TRUE),
        log2.Quantity    = dcast_diann(dt, 'PG.Quantity', fill = NA, log2 = TRUE),
        log2.Top1        = dcast_diann(dt, 'PG.Top1',     fill = NA, log2 = TRUE),
        log2.Top3        = dcast_diann(dt, 'PG.Top3',     fill = NA, log2 = TRUE),
        log2.Sum         = dcast_diann(dt, 'PG.Sum',      fill = NA             ),
        npeptide         = dcast_diann(dt, 'npeptide',    fill = 0              ),
        nprecursor       = dcast_diann(dt, 'nprecursor',  fill = 0              ), 
        sequence         = dcast_diann(dt, 'sequence',    fill = '')))
    sdt(object)$sample_id  <- snames(object)
    fdt(object)$feature_id <- fnames(object)
    analysis(object)$nfeatures <- nrow(object)
# fdt
    cols <- c('PG.MaxLFQ', 'PG.Quantity', 'PG.Top1', 'PG.Top3', 'PG.Sum', 
              'PG.Q.Value', 'sequence', 'Run', 'npeptide', 'nprecursor')
    dt[, (cols) := NULL]
    dt %<>% unique()
    assert_are_identical(nrow(dt), nrow(object)) # if not more fields need to be NULLed
    object %<>% merge_fdt(dt)
    fdt(object)$npeptide    <- rowMeans(assays(object)$npeptide)                  %>% round(digits = 1)
    fdt(object)$nprecursor  <- rowMeans(assays(object)$nprecursor)                %>% round(digits = 1)
    fdt(object)$log2.MaxLFQ <- rowMeans(assays(object)$log2.MaxLFQ, na.rm = TRUE) %>% round(digits = 1)
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
    object %<>% analyze(
        pca          = pca,           fit       = fit,
        formula      = formula,       block     = block,      
        coefs        = coefs,         contrasts = contrasts,   
        verbose      = verbose,       plot      = plot,
        feature_id   = feature_id,    sample_id = sample_id,  
        palette      = palette)
    object
}

# file <- download_data('dilution.report.tsv')
# dt <- .read_diann_proteingroups(file, sequence = TRUE)
# dcast_diann(dt, quantity = 'PG.MaxLFQ',  fill = NA, log2 = TRUE )[1:3, 1:3]
# dcast_diann(dt, quantity = 'nprecursor', fill = 0               )[1:3, 1:3]
# dcast_diann(dt, quantity = 'npeptide',   fill = 0               )[1:3, 1:3]
# dcast_diann(dt, quantity = 'sequence',   fill = ''              )[1:3, 1:3]
dcast_diann <- function(dt, quantity, fill, log2 = FALSE){
    mat <- data.table::dcast(dt, feature_id ~ Run, value.var = quantity, fill = fill)
    mat %<>% dt2mat()
    mat %<>% extract(unique(dt$feature_id), ) # preserve original order
    if (is.na(fill))  mat %<>% zero_to_na() %>% nan_to_na()
    if (log2)         mat %<>% log2()
    mat
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

