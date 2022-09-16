
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
#' @param Protein.Group string vector (without duplicates)
#' @param Protein.Names NULL or string vector (with same length as Protein.Group)
#' @param fastadt data.table
#' @return string vector
#' @examples
#' # Without file 
#'     Protein.Group <- c('Q96JP5;Q96JP5-2', 'O75822', 'Q96AC1;Q96AC1-3;Q9BQL6')
#'     Protein.Names <- c('ZFP91', 'EIF3J', 'FERM2;FERM1')
#'     forge_pg_descriptions(Protein.Group, Protein.Names)
#' # With file
#'     file <- download_data('uniprot_hsa_20140515.fasta')
#'     fastadt <- read_fastahdrs(file)
#'     forge_pg_descriptions(Protein.Group, fastadt = fastadt)
#' @export
forge_pg_descriptions <- function(
    Protein.Group, Protein.Name = NULL, fastadt = NULL
){
    assert_is_character(Protein.Group)
    assert_has_no_duplicates(Protein.Group)
    
    if (is.null(fastadt)){
        assert_is_character(Protein.Name)
        pgdt <- data.table(Protein.Group = Protein.Group, Protein.Name = Protein.Name)
        pgdt$Protein.Name %<>% stri_replace_all_regex('_[A-Z]+', '')
        pgdt[, isoform := pg_to_isoforms(Protein.Group) ]
        pgdt[, feature_id := Protein.Name]
        #pgdt[stri_detect_fixed(Protein.Name, ';'), feature_id := commonify_collapsed_strings(feature_id, ';'), by = 'feature_id']
        pgdt[isoform!='0', feature_id := paste0(feature_id, '(', isoform, ')') ]
        PG0 <- Protein.Group
        pgdt %<>% extract(PG0, on = 'Protein.Group')
        pgdt$feature_id
    } else {
        pgdt <- data.table(proId = Protein.Group, uniprot = Protein.Group)
        pgdt %<>% separate_rows(uniprot, sep = ';') %>% data.table()
        fastadt %<>% extract(, c('uniprot', 'entry', 'canonical', 'isoform'))
        pgdt %<>% merge(fastadt, by = 'uniprot', sort = FALSE)
        # pgdt %<>% drop_inferior(verbose = FALSE)
            # DIA-NN has a non-razor approach
            # P63151, Q66LE6, P63151;Q66LE6 are three different proteingroups !
            # curation is incompatible with this setup
        pgdt %<>% extract(order(proId, entry, isoform))
        pgdt %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = c('proId', 'canonical')) # Collapse PG isoforms
        pgdt[, feature_id := entry]
        pgdt[, isoform := stri_replace_all_fixed(isoform, ';', ',')]
        pgdt[isoform!='0', feature_id := paste0(feature_id, '(', isoform, ')')]
        pgdt[, isoform := NULL]
        pgdt %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = c('proId'))              # Collapse PG paralogs
        pgdt %<>% extract(Protein.Group, on = 'proId')
        pgdt$feature_id
    }
}


#' @rdname read_diann
#' @export
.read_diann_precursors <- function(file, fastadt = NULL){
# Read
    anncols <- c('Genes', 'Protein.Names', 'Protein.Group', 
                 'First.Protein.Description', 'Precursor.Id', 'Run')
    numcols <- c('PG.MaxLFQ', 'Precursor.Quantity', 'PG.Quantity')
    cols <- c(anncols, numcols)
    dt <- fread(file, select = cols)
    for (col in numcols){
        dt[[col]] %<>% stri_replace_first_fixed(',', '.') # 1977.16 but 1,35E+11
        dt[[col]] %<>% as.numeric()
    }
# Filter/Annotate
    # dt %<>% extract(Lib.Q.Value <= 0.05)
    # dt %<>% extract(Lib.PG.Q.Value <= 0.01)
    pgdt <- unique(dt[, .(Protein.Group, Protein.Names)])
    pgdt[, feature_name := forge_pg_descriptions(Protein.Group, Protein.Names, fastadt)]
    pgdt[, Protein.Names := NULL]
    dt %<>% .merge(pgdt, by = 'Protein.Group')
    dt %<>% extract(order(feature_name, Run, -Precursor.Quantity))
# Summarize/Return
    dt[, Precursor.No := seq_len(.N), by = c('Protein.Group', 'Run')]
    dt[, PG.Top1 :=     rev(sort(Precursor.Quantity))[1],                  by = c('Protein.Group', 'Run')]
    dt[, PG.Top3 := sum(rev(sort(Precursor.Quantity))[1:3], na.rm = TRUE), by = c('Protein.Group', 'Run')]
    dt[, PG.Sum  := sum(         Precursor.Quantity,        na.rm = TRUE), by = c('Protein.Group', 'Run')]
    cols <- c('Protein.Group', 'feature_name', 'First.Protein.Description', 
              'Genes', 'Run', 'Precursor.No', 'Precursor.Id',  'Precursor.Quantity', 
              'PG.Quantity', 'PG.Top1', 'PG.Top3', 'PG.Sum' , 'PG.MaxLFQ')
    dt %<>% extract(, cols, with = FALSE)
    dt
}

#' @rdname read_diann
#' @export
.read_diann_proteingroups <- function(file, fastadt = NULL){
    dt <- .read_diann_precursors(file, fastadt = fastadt)
    cols <- c('Protein.Group', 'feature_name', 'First.Protein.Description', 
              'Genes', 'Run', 'PG.Quantity', 
              'PG.Top1', 'PG.Top3', 'PG.Sum', 'PG.MaxLFQ', 'Precursor.No')
    dt %<>% extract(, cols, with = FALSE )
    dt[, N.Precursor := .N, by = c('Protein.Group', 'Run')]
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
    assert_all_are_existing_files(file)
    assert_is_identical_to_true(substr(file, nchar(file)-4, nchar(file)) == 'fasta')
    dt <- seqinr::read.fasta(file)
    dt %<>% vapply(attr, character(1), 'Annot') %>% unname()
    dt %<>% split_extract_fixed(' ', 1)
    dt %<>% substr(2, nchar(.))
    dt
}

#' Rm contaminants
#'
#' Rm contaminants from DIA-NN SumExp
#' @param object         SummarizedExperiment
#' @param contaminants   uniprots (character vector)
#' @param verbose        TRUE or FALSE
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('szymanski22.report.tsv')
#' object <- read_diann(file)
#' # object %<>% rm_diann_contaminants()
#' @export
rm_diann_contaminants <- function(
    object, contaminants = read_contaminants(), verbose = TRUE
){
    fdt0 <- fdt(object)
    fdt0[, uniprot := feature_id ]
    fdt0 %<>% separate_rows(uniprot, sep = ';') %>% data.table()
    fdt0[, `Potential contaminant`:= FALSE]
    fdt0[uniprot %in% contaminants, `Potential contaminant` := TRUE]
    fdt0[, uniprot := NULL]
    fdt0 %<>% extract(, .(`Potential contaminant` = any(`Potential contaminant`)), by = 'feature_id')
    object %<>% merge_fdata(fdt0, by.x = 'feature_id', by.y = 'feature_id')
    object %<>% filter_features(!`Potential contaminant`, verbose = verbose)
}

#' Read diann
#'
#' @param file         'report.tsv' file
#' @param fastadt       NULL or data.table
#' @param quantity     'PG.MaxLFQ', 'PG.Quantity', 'PG.Top1', 'PG.Top3', or 'PG.Sum'
#' @param contaminants  character vector with uniprot ids
#' @param plot          whether to plot
#' @param pca           whether to pca
#' @param fit           model fit engine: 'limma', 'lm', 'lmer', 'lme'
#' @param formula       model formula
#' @param block         block var (sdt)
#' @param coefs         character: coefficients to test
#' @param contrasts  character: coefficient contrasts to test
#' @param feature_id    string: summary plot feature
#' @param sample_id     string: summary plot sample
#' @param palette       character: color palette
#' @param verbose       whether to msg
#' @return  data.table / SummarizedExperiment
#' @examples 
#' # Read & Analyze
#'    file <- download_data('szymanski22.report.tsv')
#'    object <- read_diann(file)
#'    snames(object) %<>% split_extract_fixed('_', 3)
#'    object$subgroup <- as.numeric(object$sample_id)
#'    analyze(object)
#'     
#' # Read data.table (lower-level)
#'     file <- download_data('szymanski22.report.tsv')
#'    (PR   <- .read_diann_precursors(file))       # precursor    dt
#'    (PG   <- .read_diann_proteingroups(file))    # proteingroup dt
#' # Compare Summarizations
#'      PG[PG.Quantity==PG.Top1] # matches      : 26063 (85%) proteingroups
#'      PG[PG.Quantity!=PG.Top1] # doesnt match :  4534 (15%) proteingroups
#'      run <- 'IPT_HeLa_1_DIAstd_Slot1-40_1_9997'
#'      PR[Protein.Group=='Q96JP5;Q96JP5-2' & Run == run, 1:6] #    match:    8884 ==   8884
#'      PR[Protein.Group=='P36578'          & Run == run, 1:6] # no match:  650887 != 407978
#'      PR[PG.Quantity != PG.Top1][feature_name == unique(feature_name)[1]][Run == unique(Run)[1]][1:2, 1:6]
#'      PR[PG.Quantity != PG.Top1][feature_name == unique(feature_name)[2]][Run == unique(Run)[1]][1:2, 1:6]
#'      PR[PG.Quantity != PG.Top1][feature_name == unique(feature_name)[3]][Run == unique(Run)[1]][1:3, 1:6]
#' @export
read_diann <- function(
    file, fastadt = NULL, quantity = 'PG.MaxLFQ', contaminants = character(0), 
    impute = FALSE, plot = FALSE, 
    pca = plot, fit = if (plot) 'limma' else NULL, formula = NULL, block = NULL,
    coefs = NULL, contrasts = NULL, feature_id = NULL, sample_id = NULL, 
    palette = make_subgroup_palette(object), verbose = TRUE
){
# Assert
    assert_all_are_existing_files(file)
    if (!is.null(fastadt))  assert_is_data.table(fastadt)
    assert_is_subset(quantity, 
         c('PG.MaxLFQ', 'PG.Quantity', 'PG.Top1', 'PG.Top3', 'PG.Sum'))
    if (!is.null(contaminants))  assert_is_character(contaminants)
# SumExp
    # values
        dt <- .read_diann_proteingroups(file, fastadt = fastadt)
        precursors <- data.table::dcast(dt, Protein.Group ~ Run, value.var = 'N.Precursor')
        precursors %<>% dt2mat()
        dcast_quantity <- function(quantity) data.table::dcast(dt, Protein.Group ~ Run, value.var = quantity)
        assays0 <- Map(dcast_quantity, quantity)           # Map retains names - lapply doesnt
        assays0 %<>% lapply(dt2mat)
        assays0 %<>% lapply(zero_to_na, verbose = verbose) # lapply allows for multiple arguments - Map doesnt
        assays0 %<>% lapply(nan_to_na,  verbose = verbose)
        assays0 %<>% lapply(log2)
        assays0 %<>% c(list(precursors = precursors))
        object <- SummarizedExperiment(assays0)
        analysis(object)$nfeatures <- nrow(object)
    # fdt
        fdt(object)$feature_id <- rownames(object)
        fdt0 <- dt[, .(Protein.Group, feature_name, Genes, First.Protein.Description)]
        fdt0 %<>% unique()
        object %<>% merge_fdata(fdt0, by.x = 'feature_id', by.y = 'Protein.Group')
    # sdt
        snames(object) <- colnames(object)
        snames(object) %<>% drop_common()
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
          pca = pca,            fit = fit,         formula = formula, 
        block = block,        coefs = coefs,  contrasts = contrasts,
      verbose = verbose,      plot  = plot,     feature_id = feature_id,
    sample_id = sample_id,  palette = palette )
    object
}

has_one_level <- function(x) length(unique(x))==1
x <- paste0('pi_exp_', c('wt_r1', 'wt_r2', 'kd_r1', 'kd_r2'))
drop_common   <- function(x){
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

