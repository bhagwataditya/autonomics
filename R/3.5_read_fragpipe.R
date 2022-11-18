

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
    fdt1[, contaminant  := proteingroup %>% substr(1,6) %>% equals('contam') %>% as.integer()]
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
    # fdt1[gene %in% gene[duplicated(gene)]][order(gene)]
    # fdt1[uniprot %>% stri_detect_fixed(';')]
}



#' Read fragpipe
#' @param dir directory with 'combined_protein.tsv'
#' @examples 
#' file <- '../combined_protein.tsv'
read_fragpipe <- function(
    dir = getwd(), 
    file = if (is_file(dir)) dir else file.path(dir, 'combined_protein.tsv')
){
# Assert
    assert_all_are_existing_files(file)
# Read
    fdt0 <- .read_fragpipe_fdt(file)
    
}