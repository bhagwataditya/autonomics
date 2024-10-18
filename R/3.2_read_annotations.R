
#------------------------------------------------------
#
#      parse_uniprot_hdrs
#          extract_reviewed
#          extract_protein
#          extract_gene
#          extract_uniprot
#          extract_canonical
#          extract_isoform
#          extract_description
#          extract_fragment
#          extract_existence
#          FASTAFIELDS
#
#-------------------------------------------------------

extract_reviewed <- function(fastahdrs){
    fastahdrs                    %>%   # seqinr::read.fasta          gives '>sp|tr' fastahdrs
    split_extract_fixed('|',1)         %>%   # Biostrings::readAAStringSet gives  'sp|tr' fastahdrs
    substr(nchar(.)-1, nchar(.))       %>%   # this function now works with both scenarios
    equals('sp')                       %>%
    as.integer()
}

extract_protein <- function(fastahdrs){
    y <- fastahdrs %>% split_extract_fixed(' ', 1)
    idx <- stri_detect_fixed(y, '|')         # >IL2-FAT1  (645 aa)
    y[idx] %<>% split_extract_fixed('|', 3)  # >sp|P22234|PUR6_HUMAN
    #split_extract_fixed('_', 1)             # drop _HUMAN: not robust in multi-organism db!
    y
}

extract_gene <- function(fastahdrs){
    fastahdrs %>% 
    split_extract_fixed('GN=', 2)      %>% 
    split_extract_fixed(' ', 1)}

extract_uniprot <- function(fastahdrs){
    fastahdrs                    %>%
    split_extract_fixed(' ', 1)        %>% 
    split_extract_fixed('|', 2) 
}

extract_canonical <- function(fastahdrs){
    fastahdrs                    %>%
    extract_uniprot()                  %>%
    split_extract_fixed('-', 1) 
}

extract_isoform <- function(fastahdrs){
    fastahdrs                    %>%
    extract_uniprot()                  %>%
    split_extract_fixed('-', 2)        %>%
    nastring_to_0()                    %>%
    as.integer()
}

extract_description <- function(fastahdrs){
    fastahdrs                    %>%
    split_extract_regex('_[A-Z]+ ', 2) %>% 
    split_extract_regex(' OS=', 1)
}

extract_fragment <- function(fastahdrs){
    fastahdrs                    %>%
    extract_description()              %>%
    stri_detect_fixed(  'ragment')     %>%
    as.integer()
}

nastring_to_nachar <- function(x){ x[x=='NA'] <- NA_character_;  x }

extract_existence <- function(fastahdrs){
    fastahdrs                   %>%
    split_extract_fixed('PE=', 2)     %>%
    split_extract_fixed(' ', 1)       %>%
    nastring_to_nachar()              %>%
    as.integer()
}

FASTAFIELDS <- c('reviewed', 'protein', 'gene', 'fragment', 'existence', 'organism') 
    # description is too long, not so useful and also isoform specific . All other fields arent!

# UNIPROTHDR <- '(sp|tr)[|][A-Z0-9]+([-][0-9]+)?[|][^=]*[=]([A-Za-z]+ [A-Za-z]+) GN=([A-Z]+) PE=([0-5]) SV=([0-9]+)'

parse_uniprot_hdrs <- function(fastahdrs, fastafields = FASTAFIELDS){
    reviewed <- protein <- gene <- fragment <- existence <- organism <- dbid <- uniprot <- NULL
    dt <- data.table(dbid = extract_uniprot(fastahdrs) %>% split_extract_fixed('-', 1), fastahdrs = fastahdrs)
    dt %<>% extract(!duplicated(dbid))
    dt[,  uniprot := dbid ]                                                                 #   0 tr
    if ('reviewed'    %in% fastafields)  dt[, reviewed    := extract_reviewed( fastahdrs)]  #   1 sp
    if ('protein'     %in% fastafields)  dt[, protein     := extract_protein(  fastahdrs)]  # existence
    if ('gene'        %in% fastafields)  dt[, gene        := extract_gene(     fastahdrs)]  #   1 protein
    if ('fragment'    %in% fastafields)  dt[, fragment    := extract_fragment( fastahdrs)]  #   4 prediction
    if ('existence'   %in% fastafields)  dt[, existence   := extract_existence(fastahdrs)]  #   5 uncertain
    if ('organism'    %in% fastafields)  dt[, organism    := split_extract_fixed(protein, '_', 2)]
    if ('existence'   %in% fastafields)  dt[, existence   := unique(.SD)[, existence[!is.na(existence)]], by = 'uniprot'] 
        # `unique`: for phosphosites the fastahdrs are
        #  replicated when protein has multiple phosphosites
        #  This duplication needs to be eliminated before proceeding.
    dt[, fastahdrs := NULL]
    dt[]
}


#------------------------------------------------------
#
#      read_maxquant_hdrs
#      read_uniprotdt
#
#------------------------------------------------------



#' Read fasta hdrs
#' @param fastafile    string (or charactervector)
#' @param fastahdrs    character vector
#' @param fastafields  charactervector : which fastahdr fields to extract ?
#' @param force        whether to overwrite existing file
#' @param verbose      bool
#' @examples
#' # uniprot hdrs
#'      fastafile <- system.file('extdata/uniprot_hsa_20140515.fasta', package = 'autonomics')
#'      read_uniprotdt(fastafile)
#'      
#' # maxquant hdrs
#'     file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#'     dt <- .read_maxquant_proteingroups(file)
#'     parse_maxquant_hdrs(dt$`Fasta headers`)
#'
#'     profile <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
#'     fosfile <- system.file('extdata/billing19.phosphosites.txt',  package = 'autonomics' )
#'     prodt <- .read_maxquant_proteingroups(profile)
#'     fosdt <- .read_maxquant_phosphosites(fosfile, profile)
#'     parse_maxquant_hdrs(prodt$`Fasta headers`)
#'     parse_maxquant_hdrs(fosdt$`Fasta headers`)
#'     
#' # contaminant hdrs
#'     read_contaminantdt()
#'
#' @return data.table(uniprot, protein, gene, uniprot, reviewed, existence)
#' @note existence values are always those of the canonical isoform
#'       (no isoform-level resolution for this field)
#' @export
read_uniprotdt <- function(
    fastafile, fastafields = FASTAFIELDS, verbose = TRUE
){
# Assert
    if (is.null(fastafile)) return(NULL)
    assert_all_are_existing_files(fastafile)
    if (!requireNamespace('Biostrings', quietly = TRUE)){
        stop("BiocManager::install('Biostrings'). Then re-run.") }
# Read
    if (verbose)   cmessage('%sfastahdrs       %s', spaces(14), fastafile)
    fastahdrs <- Biostrings::readAAStringSet(fastafile)
    fastahdrs %<>% names()
    fastahdrs %<>% parse_uniprot_hdrs(fastafields = fastafields)
    #fastahdrs[, truncated := 0]
    fastahdrs
}


#' @rdname read_uniprotdt
#' @export
parse_maxquant_hdrs <- function(fastahdrs){
    cmessage('%smaxquant fastahdrs   ', spaces(14)) # Write Read rather than Parse to align with contaminantdt
    fastahdrs %<>% stri_split_fixed(';') %>% unlist() %>% unique()
    fastahdrs %<>% extract(. != '' )
    fastahdrs %<>% extract(stri_count_fixed(., '|')==2 ) # minimum requirement: >sp|A0AV96-2|RBM47_HUMAN
    fastadt <- parse_uniprot_hdrs(fastahdrs)
    fastadt[]
  # fastadt[           , truncated := 1 ]
  # fastadt[ idx==TRUE , truncated := 0 ]
}



#------------------------------------------------------------------------------------------
#
#    save_contaminant_hdrs
#        parse_uniprot_contaminants
#            annotate_uniprots_rest  :: .annotate_uniprot_rest ::  parse_ensp_contaminants
#            annotate_ensps_rest     :: .annotate_ensp_rest    ::    ens2org
#        parse_refseq_contaminants                             ::  taxon2org
#        parse_hinv_contaminants
#        parse_strep_contaminants
#
#    read_contaminant_hdrs
#
#------------------------------------------------------------------------------------------

#' Annotation Maps
#' @examples
#'     TAXON_TO_ORGNAME['9606']
#'    ABBREV_TO_ORGNAME['HSA']
#'  REVIEWED_TO_NUMBER['reviewed']
#' EXISTENCE_TO_NUMBER['Evidence at protein level']
#' @export
TAXON_TO_ORGNAME <-  c( `1280` = 'STAAU', 
                        `1498` = 'HATHI', 
                        `9606` = 'HUMAN', 
                        `9823` = 'PIG',
                        `9913` = 'BOVIN', 
                       `10090` = 'MOUSE', 
                      `105402` = 'ZOASP')

#' @rdname TAXON_TO_ORGNAME
#' @export
ABBREV_TO_ORGNAME <- c(   HSA = 'HUMAN',     MMU = 'MOUSE',    BTA = 'BOVIN', SSCO = 'PIG')

#' @rdname TAXON_TO_ORGNAME
#' @export
REVIEWED_TO_NUMBER <- c(unreviewed = '0', reviewed = '1')

#' @rdname TAXON_TO_ORGNAME
#' @export
EXISTENCE_TO_NUMBER <- c(    `Evidence at protein level` = 1, 
                          `Evidence at transcript level` = 2, 
                                `Inferred from homology` = 3,
                                               Predicted = 4 )

#' @rdname taxon2org
#' @export
ens2org <- function(x){
   no <- stri_match_first_regex(x, '([A-Z]+)([0-9]+)') %>% extract(, 3)
  gtp <- stri_match_first_regex(x, '([A-Z]+)([0-9]+)') %>% extract(, 2) %>% substr(nchar(.), nchar(.))
  org <- stri_match_first_regex(x, '([A-Z]+)([0-9]+)') %>% extract(, 2) %>% substr(4, nchar(.)-1)
  org[org==''] <- 'HSA'
  unname(ABBREV_TO_ORGNAME[org])
}


#' taxon/ens to organism
#' @param x character vector
#' @examples
#' taxon2org( x = c('9606', '9913') )
#'   ens2org( x = c('ENSP00000377550', 'ENSBTAP00000038329') )
#' @return character vector
#' @export
taxon2org <- function(x){
    mappings <- c(`9606` = 'HUMAN', `10090` = 'MOUSE', `9913` = 'BOVIN')
    assert_is_subset(unique(x[!is.na(x)]), names(mappings))
    unname(TAXON_TO_ORGNAME[x])
}


UNIPROTCOLS <- c('accession', 'reviewed', 'id', 'gene_primary', 'protein_existence')

# .annotate_uniprot_rest('P00761')
.annotate_uniprot_rest <- function(x, columns = UNIPROTCOLS){
    
    # Direct access to a uniprot record. 
    # Advantage: to-the-point. 
    # Disadvantages: 1) return fields cannot be customized (I think)
    #                2) redundant ids are automatically turned into their new versions
    #                   e.g. Q32MB2 is turned into Q86Y46
    #                   which leads to confusinf results
    # url <- sprintf('https://rest.uniprot.org/uniprotkb/%s.tsv', canonical0)
    cols <- paste0(columns, collapse = ',')
    url <- 'https://rest.uniprot.org/uniprotkb/search?query='
    
    # The following has been commented out because though it makes the search more specific
    # It prevents ensps from being mapped
    # if (accessions)  url %<>% paste0('accession:')
    url %<>% paste0(x)
    url %<>% paste0('&format=tsv&fields=')
    url %<>% paste0(cols)
    dt <- fread(url, colClasses = 'character', showProgress = FALSE)
                   # colClasses prevents defaulting to logical when empty
    
    # Cbind original query id for two reasons
    # First reason is that the query accession can differ from the returned accession
    # When uniprot finds an accession to be redundant it returns the new (updated) one.
    # Second reason is that the returned data.table can be empty
    # That is undesirable for further downstream data.table::rbindlist
    cbind(dbid = x, dt)
}


#' Annotate uniprot/ensp
#' @param x        character vector
#' @param columns  character vector
#' @param verbose  TRUE or FALSE
#' @return
#' data.table(dbid, uniprot, reviewed, protein, gene, canonical, 
#'            isoform, fragment, existence, organism, full)
#' @examples
#' annotate_uniprot_rest( x = c('P00761', 'Q32MB2') )
#' annotate_uniprot_rest( x = c('ENSBTAP00000006074', 'ENSP00000377550') )
#' @export
annotate_uniprot_rest <- function( x, columns = UNIPROTCOLS, verbose = TRUE ){
    organism <- fragment <- existence <- reviewed <- protein <- NULL
    if (is.null(x))  return(NULL)
    if (verbose){            cmessage('\t\t\tAnnotate %d proteins through uniprot restapi', length(x))
        if (length(x) < 10)  cmessage('\t\t\t\t%s', paste0(x, collapse = ', '))
    }

    dt <- lapply(x, .annotate_uniprot_rest, columns = columns )
    dt %<>% rbindlist()
    setnames(dt, 
         c('Entry',     'Entry Name', 'Reviewed',  'Gene Names (primary)', 'Protein existence'), 
         c('uniprot',   'protein',    'reviewed',  'gene',                 'existence'))
    
    dt[                 ,  reviewed := REVIEWED_TO_NUMBER[reviewed]        ]
    dt[ is.na(reviewed) ,  reviewed := '0'                                  ]
    dt[                 ,  reviewed :=   as.integer(reviewed)               ]
    dt[                 , existence := EXISTENCE_TO_NUMBER[existence]       ]
    dt[ is.na(existence), existence := 5                                    ]
    dt[                 ,  fragment := 0                                    ]
    dt[                 ,  organism := split_extract_fixed(protein,'_', 2)  ]

    cols <- c( 'dbid', 'uniprot', 'reviewed', 'protein',  'gene', 'existence', 'fragment', 'organism' )
    dt %<>% pull_columns(cols)
    dt[]
}




parse_refseq_contaminants <- function(hdrs){
    # refseq fastahdrs lack proteinnames
    # REFSEQ:XP_986630 Tax_Id=10090 Gene_Symbol=Krt33b keratin complex 1, acidic, gene 3
    # So I tried annotating refseq ids using web interface
    # But nearly all ids are deprecated and difficult to map to current systems
    # So instead parse whatever can be parsed from fastahdrs
    # Make up proteinname by concatenating refseqid and org
      org <- hdrs %>% stri_match_first_regex('Tax_Id=([0-9]+)') %>% extract(, 2) %>% taxon2org()
      ids <- hdrs %>% split_extract_fixed(' ', 1) %>% split_extract_fixed(':', 2)
    genes <- hdrs %>% stri_match_first_regex('Gene_Symbol=([A-Za-z0-9]+)') %>% extract(, 2)

    data.table(
               dbid = ids,
           reviewed = 0,
            protein = paste0(ids, '_', org),
            uniprot = NA_character_ ,
               gene = genes, 
           fragment = 0,
          existence = 1,
           organism = org
    )
}


parse_hinv_contaminants <- function(hdrs){
    # hinv fastahdrs lack proteinnames
    # H-INV:HIT000016045 Tax_Id=9606 Gene_Symbol=- Similar to Keratin, type II cytoskeletal 8
    # So I looked into annotating hinv ids using web interface.
    # But H-InvDB seems an obsolete thing from the past.
    # Not worth investing time in. Its also only three ids
    # Parse (whatever can be) from fastahdrs instead
    # Make up a proteinname by concatenating hinvid and org
     gene <- NULL
      ids <- hdrs %>% split_extract_fixed(' ', 1) %>% split_extract_fixed(':', 2)
     orgs <- hdrs %>% split_extract_fixed(' ', 2) %>% split_extract_fixed('=', 2) %>% taxon2org()
    genes <- hdrs %>% split_extract_fixed(' ', 3) %>% split_extract_fixed('=', 2)
    y <- data.table( 
                dbid = paste0(ids), 
             uniprot = NA_character_, 
            reviewed = 0,
             protein = sprintf('%s_%s', ids, orgs),
                gene = genes, 
            fragment = 0, 
           existence = 0, 
            organism = orgs
    )
    y[gene=='-', gene := NA_character_]
    y[]
}


parse_strep_contaminants <- function(hdrs){
    # streptavidin fastahdr lacks everthing :D
    # Streptavidin (S.avidinii)
    # So just create manually
    data.table(
                 dbid = 'Streptavidin (S.avidinii)', 
              uniprot = 'SAVIDIN', 
             reviewed = 0,
              protein = 'SAVIDIN_STRAV',
                 gene = NA_character_, 
             fragment = 0, 
            existence = 1,
             organism = 'STRAV'
    )
}


#' Save contaminant hdrs
#' 
#' Save contaminant hdrs as data.table
#' 
#' @param confile contaminant confile
#' @param verbose TRUE or FALSE
#' @return data.table
#' @examples
#' save_contaminant_hdrs()
#' @noRd
save_contaminant_hdrs <- function(confile = download_contaminants(), verbose = TRUE){
# Assert
    if (  is.null(confile))  return(NULL)
    assert_are_identical(tools::file_ext(confile), 'fasta')
    tsvfile <- file.path(dirname(confile), 'contaminants.tsv')    # dont mv up - breaks when NULL
    if (file.exists(tsvfile))  return(fread(tsvfile))
    if (!requireNamespace('Biostrings', quietly = TRUE)){
        message("BiocManager::install('Biostrings'). Then re-run.") 
        return(NULL)  }
    protein <- organism <- dbid <- NULL
# Read
    fastahdrs <- Biostrings::readAAStringSet(confile)
    fastahdrs %<>% names()                                          # 245 contaminants
# Parse
    # uniprot fastahdrs lack proteinnames
    # Q32MB2 TREMBL:Q32MB2;Q86Y46 Tax_Id=9606 Gene_Symbol=KRT73 Keratin-73
    # Therefore annotate using uniprot restapi
    if (verbose)  cmessage('\t\tParse uniprot contaminants')
    idx <- stri_detect_regex(fastahdrs, '(SWISS-PROT|TREMBL)')
    x <- split_extract_fixed(fastahdrs[idx], ' ', 1) %>% split_extract_fixed('-', 1) %>% unique()
    uniprotdt <- annotate_uniprot_rest(x)
    uniprotdt[, dbid := paste0('CON__', dbid)]
    fastahdrs %<>% extract(!idx)                                    # 210 uniprots
    
    # ensp fastahdrs lack protein/gene names
    # ENSEMBL:ENSBTAP00000038329 (Bos taurus) 9 kDa protein
    # Therefore annotate using uniprot restapi
    if (verbose)  cmessage('\t\tParse ensembl contaminants')
    idx <- stri_detect_fixed(fastahdrs, 'ENSEMBL:')
    x <- fastahdrs[idx] %>% split_extract_fixed(' ', 1) %>%  split_extract_fixed(':', 2)
    ensembldt <-  annotate_uniprot_rest(x)
    ensembldt[ is.na(organism), organism := ens2org(dbid)]
    ensembldt[ is.na(protein) , protein := paste0(dbid, '_', organism) ]
    ensembldt[, dbid := paste0('CON__ENSEMBL:', dbid) ]
    fastahdrs %<>% extract(!idx)                                    #  25 ensps (8 tracable)

    if (verbose)  cmessage('\t\tParse  refseq contaminants')
    idx <- stri_detect_fixed(fastahdrs, 'REFSEQ:')
    x <- fastahdrs[idx]
    refseqdt <- parse_refseq_contaminants( x )
    refseqdt[, dbid := paste0('CON__REFSEQ:', dbid) ]
    fastahdrs %<>% extract(!idx)                                    #   6 refseqs

    if (verbose)  cmessage('\t\tParse   hinv contaminants')
    idx <- stri_detect_fixed(fastahdrs, 'H-INV:')
    x <- fastahdrs[idx] 
    hinvdt <- parse_hinv_contaminants( x )
    hinvdt[, dbid := paste0('CON__H-INV:', dbid)]
    fastahdrs %<>% extract(!idx)                                    #   3 hinvs

    if (verbose)  cmessage('\t\tParse strept contaminant')
    idx <- stri_detect_fixed(fastahdrs, 'treptavidin')
    x <- fastahdrs[idx] 
    strepdt <- parse_strep_contaminants(x)
    strepdt[, dbid := paste0('CON__', dbid)]
    fastahdrs %<>% extract(!idx)                                    #   1 streptavidin
    assert_is_empty(fastahdrs)
# Return
    contaminantdt <- rbind( uniprotdt, ensembldt, refseqdt, hinvdt, strepdt)
    fwrite(contaminantdt, tsvfile, sep = '\t')
}


#' @rdname read_uniprotdt
#' @export
read_contaminantdt <- function(force = FALSE, verbose = TRUE){
    file <- system.file('extdata/contaminants.tsv', package = 'autonomics')
    if (verbose)  cmessage('%scontamin fastahdrs%s%s', spaces(14), spaces(35-nchar('contamin fastahdrs')), file)
    fread(file)
}



#---------------------------------------------------------------------------------------
#
#       annotate_maxquant    :    .initialize
#                                 .uncollapse
#                                 .drop_rev
#                                 .annotate           :      ..merge_hdrdt
#                                 .drop_inferior
#                                 .add_featureid
#                                 .restore_rev
#                                 .recollapse
#                                 .join
#
#---------------------------------------------------------------------------------------

.initialize <- function(dt, idcol){
    uniprot <- NULL
    anndt <- dt[, .(id = get(idcol), dbid = uniprot, reverse = '')]
    setnames(anndt, 'id', idcol)
    anndt[]
}


#' Uncollapse/Recollapse
#' 
#' Uncollapse data.table cols
#' @param dt data.table
#' @param ... cols
#' @param sep string
#' @param by  string
#' @examples
#'(dt <- data.table::data.table(
#'           uniprot  = 'Q9BQL6;Q96AC1;Q96AC1-3', 
#'           protein  = 'FERM1_HUMAN;FERM2_HUMAN', 
#'           gene     = 'FERMT1;FERMT2'))
#'(dt %<>% uncollapse(protein, gene, sep = ';'))
#'(dt %>% recollapse(by = 'uniprot')) 
#' @export
uncollapse <- function(dt, ..., sep = ';'){
    dt %>% 
    separate_rows(..., sep = sep) %>%
    data.table()
}

#' @rdname uncollapse
#' @export
recollapse <- function(dt, by, sep = ';'){
    dt[, lapply(.SD, paste_unique, collapse = sep), by = by]
}

nastring_to_0 <- function(x){ 
    x[x=='NA'] <- 0; x
        # Sometimes the canonical isoform is NOT isoform-1 !
        # https://www.uniprot.org/uniprotkb/Q9H0P0/entry#sequences
}


.uncollapse <- function(anndt, idcol, verbose){
    
    isoform <- dbid <- NULL
    anndt %<>% uncollapse(dbid, sep = ';')
    
    nproteins <- length(unique(anndt$dbid))
    nfeatures <- length(unique(anndt[[idcol]]))
    if (verbose)      cmessage('%sAnnotate  Uncollapse%s%5d %ss into  %5d  proteins', spaces(4), spaces(5), nfeatures, idcol, nproteins )
        
    idx <- anndt[ , !stri_detect_fixed(dbid, 'H-INV') ]
    anndt[!idx, isoform := 0 ]
    anndt[ idx, isoform := split_extract_fixed(dbid, '-', 2) %>% nastring_to_0() %>% as.integer()]
    anndt[ idx,    dbid := split_extract_fixed(dbid, '-', 1)]
    anndt %<>% unique()
    anndt[]
}

.drop_rev <- function(anndt, idcol, verbose){
    reverse <- dbid <- mixedgroup <- NULL
    anndt[, reverse := '' ]
    anndt[ stri_detect_fixed(dbid, 'REV__') ,  reverse := '+'   ]
    anndt[, dbid := stri_replace_first_fixed(dbid, 'REV__', '') ]
    if (verbose){
        cmessage('%sDrop REV_%s%5d  proteins', spaces(14), spaces(25), length(unique(anndt$dbid)))
    }
    anndt[, mixedgroup := any(reverse == '+') & any(reverse == ''), by = idcol]  # Drop revs in mixed groups
    nmixedgroups <- anndt[mixedgroup == TRUE, length(unique(get(idcol)))]
    if (verbose & nmixedgroups > 0)  cmessage('%sDrop rev in %d mixed groups', spaces(14), nmixedgroups)
    anndt <- anndt[ mixedgroup == FALSE | (mixedgroup == TRUE & reverse == '') ,  ]
    anndt[, mixedgroup := NULL]
    anndt[]
}


..merge_hdrdt <- function(anndt, hdrdt, idcol, verbose, first = FALSE){
    uniprot <- dbid <- NULL
    if (is.null(hdrdt)){
        if (verbose){
            if (first)  cmessage('%sAnnotate%s%5d using %s', spaces(14), spaces(33), 0, get_name_in_parent(hdrdt)) 
            else        cmessage(          '%s%5d using %s',             spaces(55), 0, get_name_in_parent(hdrdt))
        }
        return(anndt)
    }
    intersect <- anndt[ is.na(uniprot), intersect(dbid, hdrdt$dbid ) ]
    anndt0 <- anndt[!dbid %in% intersect ]
    anndt1 <- anndt[ dbid %in% intersect, c(idcol, 'dbid', 'reverse', 'isoform'), with = FALSE]
    anndt1 %<>% merge(hdrdt, by = 'dbid', sort = FALSE, all.x = TRUE)
    anndt <- rbind(anndt0, anndt1, fill = TRUE)
    if (verbose)  cmessage('%s%5d using %s', 
                                spaces(55), length(unique(anndt1$dbid)), get_name_in_parent(hdrdt))
    anndt[]
}




.annotate <- function(anndt, uniprothdrs, contaminanthdrs, maxquanthdrs, restapi, idcol, verbose){
    uniprot <- dbid <- reviewed <- fragment <- existence <- NULL
    anndt %<>% cbind( uniprot = NA_character_,  reviewed = NA_integer_,    protein = NA_character_,
                         gene = NA_character_,  fragment = NA_integer_,  existence = NA_integer_,   
                     organism = NA_character_)
    anndt %<>% ..merge_hdrdt(uniprothdrs,     idcol = idcol, verbose = verbose, first = TRUE)
    anndt %<>% ..merge_hdrdt(contaminanthdrs, idcol = idcol, verbose = verbose)
    anndt %<>% ..merge_hdrdt(maxquanthdrs,    idcol = idcol, verbose = verbose)
    
    uniprots <- if (restapi)                anndt[is.na(uniprot), dbid]       else NULL
     ngroups <- if (restapi)  length(unique(anndt[is.na(uniprot)][[idcol]]))  else 0
    uniprotrestapi <- annotate_uniprot_rest(uniprots, verbose = FALSE) # 945 proteins .. 13.44 - 
    anndt %<>% ..merge_hdrdt(uniprotrestapi, idcol = idcol, verbose = verbose)
    
    anndt[ is.na(reviewed),   reviewed := 0 ]
    anndt[ is.na(fragment),   fragment := 1 ]
    anndt[ is.na(existence), existence := 5 ]
    anndt[]
}

spaces <- function(n)  paste0(rep(' ', n), collapse = '')

.curate <- function(anndt, idcol, verbose = TRUE){
# Assert
    uniprot <- existence <- fragment <- reviewed <- NULL
    anndt %<>% copy()
    if (verbose)   cmessage('%sFilter%s%5d  proteins', spaces(14), spaces(28), length(unique(anndt$dbid)))
# Prefer swissprot (over trembl)
    n0 <- nrow(anndt)
    anndt <- anndt[, .SD[ reviewed == max(reviewed) ], by = idcol]
    if ( nrow(anndt) < n0 & verbose )    cmessage('%swithin %s%s%5d  proteins  swissprot > trembl', 
                                                  spaces(14), idcol, spaces(22), length(unique(anndt$dbid))) 
# Prefer full proteins (over fragments)
    n0 <- nrow(anndt)
    anndt <- anndt[, .SD[ fragment == min(fragment) ], by = idcol]
    if ( nrow(anndt) < n0 & verbose )    cmessage('%s%5d  proteins  fullseq   > fragment', 
                                                  spaces(48), length(unique(anndt$dbid)))
# Prefer better existences (over worse)
    n0 <- nrow(anndt)
    anndt <- anndt[, .SD[existence == min(existence)], by = idcol]
    if ( nrow(anndt) < n0 & verbose )    cmessage('%s%5d  proteins  %s', spaces(48), 
                                                            length(unique(anndt$dbid)), 
                           'protein   > transcript > homolog > prediction > uncertain')
# Order
    anndt[, c('reviewed', 'fragment', 'existence') := NULL]
    idnumeric <- is_numeric_string(anndt[[idcol]][1])
    if (idnumeric){ anndt %<>% extract(order(as.integer(get(idcol)), uniprot))
    } else {        anndt %<>% extract(order(           get(idcol) , uniprot)) }
# Return
    anndt[]
}

.restore_rev <- function(anndt, idcol, verbose){
    reverse <- uniprot <- protein <- gene <- NULL
    if (verbose)  cmessage('%sAdd REV__%s%5d  proteins', spaces(14), spaces(25), length(unique(anndt$dbid)))
    anndt[ reverse == '+',    uniprot := paste0('REV__', uniprot   ) ]
    anndt[ reverse == '+',    protein := paste0('REV__', protein   ) ]
    anndt[ reverse == '+',       gene := paste0('REV__', gene      ) ]
    anndt[, c('reverse', 'dbid') := NULL ]
    anndt
}


.recollapse <- function(anndt, idcol, verbose){
    uniprot <- isoform <- NULL
    anndt %<>% extract(order(uniprot, isoform)) # we want the isoforms to be pasted in order!
    anndt %<>% recollapse(by = idcol, sep = ';')
    anndt %<>% pull_columns(c(idcol, 'uniprot', 'isoform', 'protein'))
    if (verbose)  cmessage('%sRecollapse%s%5d %ss', spaces(14), spaces(5), 
                           length(unique(anndt[[idcol]])), idcol)
    anndt
}


#' Clean Merge
#' @param dt1 data.table
#' @param dt2 data.table
#' @param by string
#' @examples
#' require(data.table)
#' dt1 <- data.table(feature_id = c('PG1', 'PG2'), gene    = c('G1', 'G2'))
#' dt2 <- data.table(feature_id = c('PG1', 'PG2'), protein = c('P1', 'P2'))
#' dt1 %<>% .merge(dt2, by = 'feature_id')
#' dt1
#' @export
.merge <- function(dt1, dt2, by){
    dt1 %<>% copy()
    dt2 %<>% copy()
    commoncols <- intersect(names(dt1), names(dt2))
    commoncols %<>% setdiff(by)
    if (length(commoncols) > 0)  dt1[, (commoncols) := NULL]
    dt1 %<>% merge(dt2, by = by, sort = FALSE, all.x = TRUE)
    dt1 %<>% pull_columns(names(dt2))
    dt1[]
}


.join <- function(dt, anndt, idcol){
    uniprot <- Original <- NULL
    setnames(dt, 'uniprot', 'Original')
    dt %<>% .merge(anndt, by = idcol)
    dt[is.na(uniprot), uniprot := Original]
    dt[, Original := NULL]
    dt %<>% extract(order(as.integer(get(idcol))))
    dt
}

.name <- function(dt, idcol){
# Initialize
    `Amino acid`                <- contaminant <- isoform <- protein <- NULL
    `Positions within proteins` <- reverse     <- NULL
# Add
    dt %<>% copy()
    dt[, contaminant := as.character(contaminant)]   # one dataset had NA_logical values
    dt[is.na(contaminant), contaminant := '']        # which made all features being filtered out in the next lines
    dt[, feature_id := protein]
  # contaminants come from multiple origin
  # so with new generalized approach this doesnt work anymore
  # if (length(unique(dt$organism)) == 1)  dt1$feature_id %<>% stri_replace_all_regex('_[^;]+', '')
    dt[isoform!=0, feature_id := paste0(feature_id, '(', stri_replace_all_fixed(isoform, ';', ''), ')')]
    if (idcol=='fosId'){
        dt[, feature_id := paste0(feature_id, '-', `Amino acid`) ]
        dt[, feature_id := paste0(feature_id, split_extract_fixed(`Positions within proteins`, ';', 1)) ]
    }
    dt[contaminant=='+', feature_id := paste0('CON__', feature_id)]
    dup <- dt[duplicated(feature_id), feature_id]
    dt[ feature_id %in% dup, feature_id := paste0(feature_id, '_', get(idcol)) ]
    dt %<>% pull_columns(c(idcol, 'feature_id'))
    assert_all_are_non_missing_nor_empty_character(dt$feature_id)
    dt[]
}


#' Annotate maxquant
#' 
#' Annotate maxquant data.table
#' 
#' Uncollapse, annotate, curate, recollapse, name
#' @param dt              `data.table` : output of `read_maxquant_(proteingroups|phosphosites)`
#' @param uniprothdrs     `data.table` : output of `read_uniprotdt`
#' @param contaminanthdrs `data.table` : output of `read_uniprotdt`
#' @param maxquanthdrs    `data.table` : output of `read_uniprotdt`
#' @param restapi         `logical(1)` : use uniprot restapi to complete missing annotations ?
#' @param verbose         `logical(1)` : message ?
#' @return \code{data.table}
#' @examples
#' # Fukuda 2020: contaminants + maxquanthdrs
#' #-----------------------------------------
#'           file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#'             dt <- .read_maxquant_proteingroups(file)
#'             dt[, 1:2]
#'      uniprothdrs <- NULL
#'  contaminanthdrs <- read_contaminantdt()
#'     maxquanthdrs <- parse_maxquant_hdrs(dt$`Fasta headers`); dt$`Fasta headers` <- NULL
#'           dt %<>% annotate_maxquant(uniprothdrs, contaminanthdrs, maxquanthdrs)
#'           dt[                 , 1:9]
#'           dt[    reverse== '+', 1:9]
#'           dt[contaminant== '+', 1:9]
#'                                               
#' # Billing 2019: uniprothdrs + contaminants + maxquanthdrs
#' #--------------------------------------------------------
#' profile <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
#' fosfile <- system.file('extdata/billing19.phosphosites.txt',  package = 'autonomics')
#'  upfile <- system.file('extdata/uniprot_hsa_20140515.fasta',  package = 'autonomics')
#' prodt <- .read_maxquant_proteingroups(profile);         prodt[, 1:2]
#' fosdt <- .read_maxquant_phosphosites(fosfile, profile); fosdt[, 1:3]
#'     uniprothdrs <- read_uniprotdt(upfile)
#' contaminanthdrs <- read_contaminantdt()
#'    maxquanthdrs <- parse_maxquant_hdrs(prodt$`Fasta headers`)
#' annotate_maxquant(prodt, uniprothdrs, contaminanthdrs, maxquanthdrs)[, 1:8]
#' annotate_maxquant(fosdt, uniprothdrs, contaminanthdrs, maxquanthdrs)[, 1:8]
#' @md
#' @export
annotate_maxquant <- function(
    dt, uniprothdrs, contaminanthdrs, maxquanthdrs, restapi = FALSE, verbose = TRUE
){
# Assert
    dt %<>% copy()
    if (!is.null(uniprothdrs))  assert_fastadt(uniprothdrs)
    idcol <- if ('fosId' %in% names(dt)) 'fosId' else 'proId'
# Annotate
    anndt  <-  .initialize( dt, idcol )
    anndt %<>% .uncollapse(     idcol,   verbose = verbose )
    anndt %<>% .drop_rev(       idcol,   verbose = verbose )   # drops REV__ and drops rev in mixed groups
    anndt %<>% .annotate(            uniprothdrs = uniprothdrs,
                                 contaminanthdrs = contaminanthdrs,
                                    maxquanthdrs = maxquanthdrs, 
                                         restapi = restapi, 
                                           idcol = idcol, 
                                         verbose = verbose)
    anndt %<>% .curate(         idcol,   verbose = verbose )
    anndt %<>% .restore_rev(    idcol,   verbose = verbose)
    anndt %<>% .recollapse(     idcol,   verbose = verbose)
       dt %<>% .join(    anndt, idcol )
       dt %<>% .name(           idcol)
# Return
    dt[]
}



#---------------------------------------------------------------------------
#
#                       annotate_uniprot_ws
#
#---------------------------------------------------------------------------

paste_unique <- function(x, collapse) paste0(unique(x), collapse=collapse)

.annotate_uniprot_ws <- function(fdt, upws, columns='ENSEMBL', collapse=';'){
    # Map uniprot -> columns
        resdt <- data.table(UniProt.ws::select(upws, fdt$uniprot, columns))
    # Trim whitespace
        resdt <- resdt[, lapply(.SD, trimws),                          by = 'From']
    # Collapse per uniprot
        resdt <- resdt[, lapply(.SD, paste_unique, collapse=collapse), by = 'From']
    # Ensure original order
        merge.data.table(
            fdt, resdt, by.x = 'uniprot', by.y = 'From', all.x = TRUE, sort = FALSE)
}

#' Annotate uniprotids using UniProt.ws
#'
#' Annotate uniprotids in data.table or SummarizedExperiment
#'
#' data.table: column "uniprot" with uniprotid values (e.g. "Q5JTY5"). \cr
#'
#' SummarizedExperiment:  svar "feature_id" with collapsed uniprot values for 
#' a single protein ("A0A068F9P7") or a proteingroup ("A0A068F9P7;F1Q7I4").
#' On these the mapping is performed. 1:many mappings are collapsed and only then returned.
#' @param x          data.table/SummarizedExperiment with (s)var "uniprot" \cr
#'                   with single (data.table) or collapsed (SummarizedExperiment) uniprot values
#' @param upws       return value of Uniprot.ws::Uniprot.ws()
#' @param columns    subset of UniProt.ws::columns(up)
#' @param collapse   string: used in paste(collapse = .)
#' @param ...        used for S3 generic definition
#' @return SummarizedExperiment/data.table
#' @examples
#' # data.table
#'     x <- data.table::data.table(uniprot = c('A0A068F9P7', 'Q7ZVA2'))
#'     # upws <- UniProt.ws::UniProt.ws(taxId=7955)
#'     # annotate_uniprot_ws(x, upws)
#'
#' # SummarizedExperiment
#'     file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#'     x <- read_maxquant_proteingroups(file)
#'     x %<>% extract(1:10, )
#'     fdt(x)[1:3, ]
#'     # x %<>% annotate_uniprot_ws(upws)
#'     # fdt(x)[1:3, ]
#' @noRd
annotate_uniprot_ws <- function(x, ...)  UseMethod('annotate_uniprot_ws')


#' @noRd
annotate_uniprot_ws.data.table <- function(
    x, upws, columns=c('xref_ensembl'), collapse=';', ...
){
    # Assert valid inputs
        if (!requireNamespace('UniProt.ws', quietly = TRUE)){
            message("`BiocManager::install('UniProt.ws')`. Then re-run.")
            return(x)
        }
        assert_is_data.table(x)
        assert_is_subset('uniprot', names(x))
        assert_is_identical_to_true(all(!is.na(x$uniprot)))
        assert_is_all_of(upws, 'UniProt.ws')
        assert_is_character(columns)
        assert_is_subset(columns, UniProt.ws::columns(upws))
        assert_is_a_string(collapse)
        columns %<>% setdiff(names(x))
    # Chunk data (to avoid UniProt.ws warning)
        n <- ceiling(nrow(x)/99)
        chunks <- rep(seq_len(n), each=99)
        chunks %<>% extract(seq_len(nrow(x)))
        chunk <- NULL
        x[, chunk := chunks]
    # Call backend for each chunk
        x <- x[, .annotate_uniprot_ws(.SD, upws, columns = columns, collapse = collapse), by = 'chunk']
    # Return
        x$chunk <- NULL
        x[]
}


#' @noRd
annotate_uniprot_ws.SummarizedExperiment <- function(
    x, upws, columns = c('ENSEMBL'), collapse=';', ...
){
    # Assert
        assert_is_all_of(upws, 'UniProt.ws')
    # Split proteingroups into proteins and extract uniprot
        fdt <- data.table(fdata(x))
        fdt %<>% uncollapse(uniprot, sep=collapse)
        fdt[, uniprot := split_extract_fixed(uniprot, '-', 1)]
    # Annotate
        uniprot <- NULL
        # fdt %<>% extract(uniprot %in% UniProt.ws::keys(upws, keytype = 'UniProtKB'))  # latest version returns only small number of keys!
        fdt %<>% annotate_uniprot_ws.data.table(upws, columns = columns, collapse = collapse)
    # Collapse
        fdt <- fdt[, lapply(.SD, trimws),                          by = 'feature_id'] # trim whitespace
        fdt <- fdt[, lapply(.SD, paste_unique, collapse=collapse), by = 'feature_id'] # collapse
        x %<>% merge_fdt(fdt)
        x
}


#---------------------------------------------------------------------------
#
#                      add_psp
#
#---------------------------------------------------------------------------

#' Add psp
#' 
#' Add PhosphoSitePlus literature counts
#' 
#' Go to www.phosphosite.org                   \cr
#' Register and Login.                         \cr
#' Download Phosphorylation_site_dataset.gz'.  \cr
#' Save into: file.path(R_user_dir('autonomics','cache'),'phosphositeplus')
#' @param object     SummarizedExperiment
#' @param pspfile    phosphositeplus file
#' @return  SummarizedExperiment
#' @examples 
#' fosfile <- system.file('extdata/billing19.phosphosites.txt',  package = 'autonomics')
#' profile <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_phosphosites(fosfile = fosfile, profile = profile)
#' fdt(object)
#' object %<>% add_psp()
#' fdt(object)
#' @export
add_psp <- function(
    object, 
    pspfile = file.path(R_user_dir('autonomics', 'cache'), 
            'phosphositeplus', 'Phosphorylation_site_dataset.gz')
){
# Initialize
    AA <- ACC_ID <- `Amino acid` <- LT_LIT <- MOD_RSD <- MS_CST <- MS_LIT <- psp <- NULL
# Read
    assert_is_all_of(object, 'SummarizedExperiment')
    if (!file.exists(pspfile))  return(object)
    dt <- data.table::fread(pspfile)
    dt[is.na(LT_LIT), LT_LIT := 0]
    dt[is.na(MS_LIT), MS_LIT := 0]
    dt[is.na(MS_CST), MS_CST := 0]
    dt[, psp := LT_LIT + MS_LIT + MS_CST]
    dt$AA <- dt$MOD_RSD %>% substr(1, 1)
    dt$MOD_RSD %<>% substr(2, nchar(.)-2)
# Rm duplicates
    idx <- cduplicated(dt[, .(ACC_ID, MOD_RSD)])
    # dt[idx]  # Seems like a (single) psiteplus bug: two different rows for the same site
    #          # On the website it says 'under review' for this particular site
    #          # Rm these
    dt %<>% extract(!idx)
    dt %<>% extract(, .(ACC_ID, MOD_RSD, psp, AA))
    # Look at example
    # dt %<>% extract('ZZZ3',   on = 'PROTEIN')
    # dt %<>% extract(c('S606-p', 'S613-p'), on = 'MOD_RSD')
# Merge
    fdt(object)$uniprot1  <- split_extract_fixed(fdt(object)$Curated, ';', 1)
    fdt(object)$Position1 <- split_extract_fixed(fdt(object)$`Positions within proteins`, ';', 1)
    object %<>% merge_fdt(dt, by.x = c('uniprot1', 'Position1'), 
                                by.y = c('ACC_ID', 'MOD_RSD'))
    idx <- fdt(object)[, !is.na(psp) &  AA != `Amino acid`]
    fdt(object)$psp[idx] <- NA
    fdt(object)$AA <- NULL
    object
}



