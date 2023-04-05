

.read_fragpipe_fdt <- function(file){
# Read
    fdt0 <- fread(file, select = c('Protein', 'Indistinguishable Proteins', 'Protein Existence', 'Gene'))
    names(fdt0) <- c('fastahdr', 'proteingroup', 'existence', 'gene')
    fdt0[proteingroup != '', proteingroup := paste0(fastahdr, ', ', proteingroup)] 
    fdt0[proteingroup == '', proteingroup := fastahdr]
    fdt0$existence %<>% substr(1,1) %>% as.integer()
# Uncollapse    
    fdt1 <- uncollapse(fdt0, proteingroup, sep = ', ')
    fdt1[, protein      := proteingroup %>% split_extract_fixed('|', 3)]
    fdt1[, organism     := protein      %>% split_extract_fixed('_', 2)]
    fdt1[, uniprot      := proteingroup %>% split_extract_fixed('|', 2)]
    fdt1[, contaminant  := '']; fdt1[substr(proteingroup, 1,6) == 'contam', contaminant := '+']
    fdt1[, reviewed     := proteingroup %>% substr(1,2) %>% equals('sp')     %>% as.integer()]
    fdt1[, proteingroup := NULL ]
# Filter    
    fdt1 <- fdt1[, .SD[contaminant ==1 |  reviewed==max(reviewed) ], by = 'fastahdr']
    fdt1 <- fdt1[, .SD[contaminant ==1 | existence==min(existence)], by = 'fastahdr']
    fdt1[, reviewed := NULL]
    fdt1[, existence := NULL]
# Collapse
    fdt1 %<>% extract(order(protein))  # ensures proper order: VLCA VLCB  ...
    fdt1[, gene     := paste_unique(gene,       collapse = ';'), by = 'fastahdr']
    fdt1[, uniprot  := paste_unique(uniprot,    collapse = ';'), by = 'fastahdr']                 # ZNG1A_HUMAN     ZNG1B_HUMAN     ZNG1B_MOUSE
    fdt1[, protein  := protein %>% split_extract_fixed('_', 1)]                                   # ZNG1A           ZNG1B           ZNG1B
    fdt1[, protein  := commonify_strings(protein),               by = c('fastahdr', 'organism')]  # ZNG1(A|B)
    fdt1[, protein  := protein %>% paste0('_', organism) ]                                        # ZNG1(A|B)_HUMAN                 ZNG1B_MOUSE
    fdt1[, protein  := protein %>% paste_unique(collapse = ';'), by = 'fastahdr']                 # ZNG1(A|B)_HUMAN,ZNG1B_MOUSE
    fdt1[, organism := paste_unique(organism,   collapse = ';'), by = 'fastahdr']
    fdt1 %<>% unique()
# feature_id
    fdt1 %<>% extract(fdt0$fastahdr, on = 'fastahdr')
    fdt1[, feature_id := protein]
    fdt1 %<>% pull_columns(c('feature_id', 'protein', 'organism', 'gene'))
    fdt1
    # fdt1[gene %in% gene[duplicated(gene)]][order(gene)]
    # fdt1[uniprot %>% stri_detect_fixed(';')]
}

.read_fragpipe_mat <- function(file, pattern){
    select <- cols(file) %>% extract(stri_detect_regex(., pattern))
    select %<>% c('Protein', .)
    mat <- fread(file, select = select, integer64 = 'numeric') %>% un_int64()
    mat %<>% dt2mat()
    colnames(mat) %<>% stri_replace_first_regex(pattern, '')
    mat
}

.stri_any_regex <- function(str, pattern) any(stri_detect_regex(str, pattern))


#' Does any string have a regex
#' @param str      string vector
#' @param pattern  string
#' @return TRUE or FALSE
#' @examples 
#' str <- c('s1 Spectral Count', 's1 Unique Spectral Count')
#' patterns <- c('Spectral Count', '(?<!Unique) Spectral Count', 'Intensity')
#' stri_detect_regex(str, pattern = patterns[1])
#' stri_detect_regex(str, pattern = patterns[2])
#' stri_detect_regex(str, pattern = patterns[3])
#' stri_any_regex(   str, pattern = patterns)
#' @export
stri_any_regex <- function(str, pattern){
    vapply(pattern, function(pat) .stri_any_regex(str, pat), logical(1))
}



#' Read fragpipe
#' @param dir           directory with 'combined_protein.tsv'
#' @param file         'combined_protein.tsv' (full path)
#' @param contaminants  whether to include contaminants
#' @param verbose       whether to msg
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('multiorganism.combined_protein.tsv')
#' object <- read_fragpipe(file = file)
#' assayNames(object)
#' fdt(object)
#' @export
read_fragpipe <- function(
    dir = getwd(), 
    file = if (is_file(dir)) dir else file.path(dir, 'combined_protein.tsv'), 
    contaminants = FALSE,
    verbose = TRUE
){
# assays
    assert_fragpipe_tsv(file)
    pattern <- c(maxlfq         = ' MaxLFQ Intensity', 
                 intensity      = '(?<!MaxLFQ) Intensity', 
                 spectralcounts = '(?<!(Combined|Unique|Total)) Spectral Count',
                 uniquecounts   = '(?<!Combined) Unique Spectral Count', 
                 totalcounts    = '(?<!Combined) Total Spectral Count')
    pattern %<>% extract(stri_any_regex(cols(file), .))
    object <- mapply(.read_fragpipe_mat, pattern = pattern, MoreArgs = list(file = file), SIMPLIFY = FALSE)
    object %<>% SummarizedExperiment()
    intensity_assays <- names(pattern)
    intensity_assays %<>% intersect(c('maxlfq', 'intensity'))
    for (ass in intensity_assays){  
        assays(object)[[ass]] %<>% zero_to_na()
        object %<>% log2transform(assay = ass, verbose = TRUE)
    }
# sdt
    sdt(object) <- data.table(sample_id = colnames(object), subgroup = 'group0')
    object$subgroup <- infer_subgroup( object$sample_id)
# fdt/sdt
    fdt0 <- .read_fragpipe_fdt(file)
    assert_all_are_true(fdt0$fastahdr == rownames(object))
    rownames(object) <- fdt0$feature_id
    fdt(object) <- fdt0
    if (!contaminants)  object %<>% filter_features(contaminant == '', verbose = verbose)
    for (assay in assayNames(object))  object %<>% add_assay_means(assay)
# return
    object 
}
