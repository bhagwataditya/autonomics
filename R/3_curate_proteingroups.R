uncollapse <- function(dt, ..., sep = ';'){
    dt %>% 
    separate_rows(..., sep = sep) %>%
    data.table()
}

nastring_to_nachar <- function(x){ x[x=='NA'] <- NA_character_;  x }

nastring_to_1 <- function(x){ x[x=='NA'] <- 1; x}

extract_reviewed <- function(`Fasta headers`){
    `Fasta headers` %>% 
    substr(2, 3)    %>%
    equals('sp')    %>%
    as.integer()
}

extract_entryname <- function(`Fasta headers`){
    `Fasta headers` %>% 
    split_extract_fixed(' ', 1) %>% 
    split_extract_fixed('|', 3) %>% 
    split_extract_fixed('_', 1)
}

extract_genename <- function(`Fasta headers`){
    `Fasta headers` %>% 
    split_extract_fixed('GN=', 2)      %>% 
    split_extract_fixed(' ', 1)}

extract_uniprot <- function(`Fasta headers`){
    `Fasta headers` %>%
    split_extract_fixed(' ', 1) %>% 
    split_extract_fixed('|', 2) 
}

extract_canonical <- function(`Fasta headers`){
    `Fasta headers`   %>%
    extract_uniprot() %>%
    split_extract_fixed('-', 1) 
}

extract_isoform <- function(`Fasta headers`){
    `Fasta headers`   %>%
    extract_uniprot() %>%
    split_extract_fixed('-', 2) %>%
    nastring_to_1() %>%
    as.integer()
}

extract_proteinname <- function(`Fasta headers`){
    `Fasta headers` %>%
        split_extract_regex('_[A-Z]+ ', 2) %>% 
        split_extract_regex(' OS=', 1)
}


extract_fragment <- function(`Fasta headers`){
    `Fasta headers`                %>%
        extract_proteinname()          %>%
        stri_detect_fixed(  'ragment') %>%
        as.integer()
}

extract_existence <- function(`Fasta headers`){
    `Fasta headers`      %>%
    split_extract_fixed('PE=', 2) %>%
    split_extract_fixed(' ', 1) %>%
    nastring_to_nachar() %>%
    as.integer()
}
parse_fastahdrs <- function(`Fasta headers`){
    dt <- data.table(                                           # Reviewed
        Reviewed        = extract_reviewed(   `Fasta headers`), #   0 tr
        `Entry names`   = extract_entryname(  `Fasta headers`), #   1 sp
        `Gene names`    = extract_genename(   `Fasta headers`),
        Uniprot         = extract_uniprot(    `Fasta headers`), # Existence 
        Canonical       = extract_canonical(  `Fasta headers`), #   1 protein 
        Isoform         = extract_isoform(    `Fasta headers`), #   1 transcript
        `Protein names` = extract_proteinname(`Fasta headers`), #   2 homolog
        Fragment        = extract_fragment(   `Fasta headers`), #   3 prediction
        Existence       = extract_existence(  `Fasta headers`)) #   4 uncertain
    dt[, Existence := unique(.SD)[, Existence[!is.na(Existence)]], by = 'Canonical']
    # `unique`: for phosphosites the Fasta headers are sometimes
    #  replicated (when protein has multiple phosphosites)
    #  This duplication needs to be eliminated before proceeding.
    dt[]
}


# tests/testthat/test_3_curate_proteingroups.R
#' Curate annotations
#' @param dt data.table
#' @return data.table
#' @examples
#' profile <- download_data('billing19.proteingroups.txt')
#' fosfile <- download_data('billing19.phosphosites.txt')
#' fastafile <- download_data('uniprot_hsa_20140515.fasta')
#' prodt <- .read_pro(profile)
#' fosdt <- .read_fos(fosfile)
#' curate_proteingroups(prodt)[, 1:8]
#' curate_proteingroups(fosdt)[, 1:8]
#' curate_proteingroups(prodt, fastafile = fastafile)[, 1:8]
#' curate_proteingroups(fosdt, fastafile = fastafile)[, 1:8]
#' @export
curate_proteingroups <- function(dt, fastafile = NULL, verbose = TRUE){
    # Split
    if (verbose)  message('\tCurate proteingroups')
    if (verbose)  message('\t\tRm Contaminants/Reverse')
    idx <- dt$Reverse == '+' | dt$`Potential contaminant` == '+' 
    conrevdt <- dt[idx]
    dt %<>% extract(!idx)
    # Curate using MaxQuant
    if (verbose)  message('\t\tMaxQuant Fastahdrs')
    idcol <- if ('fosId' %in% names(dt)) 'fosId' else 'proId'
    anndt <- dt[, c(idcol, 'Fasta headers'), with = FALSE]
    if (verbose)  message('\t\t\tUncollapse, Drop truncated, Parse')
    anndt %<>% uncollapse(`Fasta headers`, sep = ';')
    anndt %<>% extract( stri_count_fixed(`Fasta headers`, '|') == 2)
    anndt %<>% cbind(parse_fastahdrs(anndt$`Fasta headers`))
    anndt$`Fasta headers` <- NULL
    anndt %<>% drop_inferior()
    
    if (verbose)  message('\t\t\tCollapse')
    anndt %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = idcol)
    # Curate using fastafile
    if (!is.null(fastafile)){
        if (verbose)  message('\t\tFastafile Fastahdrs')
        if (verbose)  message('\t\t\tUncollapse Uniprot accessions')
        anndt2 <- dt[, c(idcol, 'Uniprot'), with = FALSE]
        anndt2 %<>% uncollapse(`Uniprot`, sep = ';')
        
        if (verbose)  message('\t\t\tRead/Parse Fastafile ', fastafile)
        fastahdrs <- seqinr::read.fasta(fastafile)
        fastahdrs %<>% vapply(attr, character(1), 'Annot') %>% unname()
        fastadt <- parse_fastahdrs(fastahdrs)
        anndt2 %<>% merge(fastadt, by = 'Uniprot', sort = FALSE)
        anndt2 %<>% drop_inferior()
        
        if (verbose)  message('\t\t\tCollapse')
        anndt2 %<>% extract(, lapply(.SD, paste_unique, collapse = ';'), by = idcol)
        
        missingInFasta  <- setdiff(  dt[[idcol]], anndt2[[idcol]])
        if (verbose)  message('\t\t\t', length(missingInFasta), 
                              ' proteins missing in fastafile - revert to MaxQuant headers for these')
        anndt %<>% extract(missingInFasta,  on = idcol)
        anndt %<>% rbind(anndt2)
        anndt %<>% extract(order(as.integer(get(idcol))))
    }
    # Merge back annotations
    if (verbose)  message('\t\tBring back Contaminants/Reverse')
    setnames(anndt, 'Uniprot', 'Curated')
    dt %<>% merge(anndt, by = idcol, sort = FALSE)
    newcols <- setdiff(names(dt), names(conrevdt))   
    conrevdt[, (newcols) := '']
    dt %<>% rbind(conrevdt)
    dt %<>% pull_columns(names(anndt))
    dt$`Fasta headers` <- NULL
    dt[order(as.integer(get(idcol)))]
}

drop_inferior <- function(anndt, verbose = TRUE){
    idcol <- if ('fosId' %in% names(anndt)) 'fosId' else 'proId'
    
    if (verbose)  message('\t\t\tDrop proteins with NA `Entry names` / `Protein names`')
    anndt <- anndt[`Entry names` != 'NA']
    anndt <- anndt[`Protein names` != 'NA']
    
    if (verbose)  message('\t\t\tWithin ', idcol, ': drop trembl    in favour of swissprot')
    anndt <- anndt[, .SD[ Reviewed == max(Reviewed) ], by = idcol]
    
    if (verbose)  message('\t\t\t       ','     ','  drop fragments in favour of full proteins')
    anndt <- anndt[, .SD[ Fragment == min(Fragment) ], by = idcol]
    
    if (verbose)  message('\t\t\t       ','     ','  drop worse existences', 
                          ' (1=protein, 2=transcript, 3=homolog, 4=prediction, 5=uncertain)')
    if (!any(is.na(anndt$Existence))){
        anndt <- anndt[, .SD[Existence == min(Existence)], by = idcol] }
    anndt[, c('Reviewed', 'Fragment', 'Existence') := NULL]
    
    if (verbose)  message("\t\t\tDrop 'Isoform x of ' in `Protein names`")
    anndt$`Protein names` %<>% stri_replace_first_regex('Isoform [0-9A-Z]+ of ', '')
    
    if (verbose)  message("\t\t\tDrop '(Fragment)'    in `Protein names`")
    anndt$`Protein names` %<>% stri_replace_first_fixed(' (Fragment)', '')
    anndt[order(as.integer(get(idcol)), Canonical, Isoform)]
}
