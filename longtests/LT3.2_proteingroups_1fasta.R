#============================================================================
#
#                     read_fastahdrs
#                    parse_fastahdrs
#
#============================================================================

context('`read_fastahdrs/parse_fastahdrs`')
fastafile <- download_data('uniprot_hsa_20140515.fasta')
fastadt <- read_fastahdrs(fastafile)

test_that('`read_fastahdrs` reads all lines', {
    nrow(fastadt)             # 88 698
    x <- readChar(fastafile, file.info(fastafile)$size)
    x %<>% substr(2, nchar(.))
    x %<>% stri_split_regex('(\r)?(\n)[>]') %>% unlist()
    expect_identical(length(x), nrow(fastadt))
})

test_that('`read_fastahdrs` reads first protein',
    expect_equal(
        fastadt['P31946', on = 'uniprot'],
        data.table(
            uniprot     = 'P31946', 
            reviewed    = 1,
            protein     = '1433B_HUMAN', 
            gene        = 'YWHAB',
            canonical   = 'P31946', 
            isoform     = 0,
            fragment    = 0,
            existence   = 1, 
            organism    = 'HUMAN'))
)

test_that('`read_fastahdrs` reads intermediate swissprot protein',
    expect_equal(
        fastadt['Q9BUJ2-4', on = 'uniprot'],
        data.table(
            uniprot     = 'Q9BUJ2-4', 
            reviewed    = 1,
            protein     = 'HNRL1_HUMAN', 
            gene        = 'HNRNPUL1',
            canonical   = 'Q9BUJ2', 
            isoform     = 4,
            fragment    = 0,
            existence   = 1, 
            organism    = 'HUMAN'))
)

test_that('`read_fastahdrs` reads intermediate trembl protein',
    expect_equal(
        fastadt['G5E9N3', on = 'uniprot'],
        data.table(
            uniprot     = 'G5E9N3', 
            reviewed    = 0,
            protein     = 'G5E9N3_HUMAN', 
            gene        = 'RETSAT',
            canonical   = 'G5E9N3', 
            isoform     = 0,
            fragment    = 0,
            existence   = 4, 
            organism    = 'HUMAN'))
)

test_that('`read_fastahdrs` reads last trembl protein',
    expect_equal(
        fastadt['R4GMM2', on = 'uniprot'],
        data.table(
            uniprot     = 'R4GMM2', 
            reviewed    = 0,
            protein     = 'R4GMM2_HUMAN', 
            gene        = 'PARD6A',
            canonical   = 'R4GMM2', 
            isoform     = 0,
            fragment    = 0,
            existence   = 4, 
            organism    = 'HUMAN'))
)


#============================================================================
#
#       maxquant_curate
#          fasta_curate
#
#============================================================================

context('`curate_annotate_maxquant/curate_annotate_fastafile`')
proteinfile <- download_data('billing19.proteingroups.txt')
phosphofile <- download_data('billing19.phosphosites.txt')
pro <- .read_maxquant_proteingroups(proteinfile)
fos <- .read_maxquant_phosphosites(phosphofile, proteinfile, quantity = quantity)
pro1 <- curate_annotate_maxquant(pro)
fos1 <- curate_annotate_maxquant(fos)
pro2 <- curate_annotate_fastafile(pro, fastadt)
fos2 <- curate_annotate_fastafile(fos, fastadt)
anncols <- c('protein', 'gene','canonical', 'isoform', 'organism')

test_that('`curate` preserves rows',     {
    expect_equal( nrow(pro1), nrow(pro))  # pro maxquant
    expect_equal( nrow(pro2), nrow(pro))  #     fasta
    expect_equal( nrow(fos1), nrow(fos))  # fos maxquant
    expect_equal( nrow(fos2), nrow(fos))  #     fasta
})

test_that('`curate` preserves cols',     { 
    expect_equal( setdiff(names(pro),  names(pro1)), character(0))  # pro maxquant
    expect_equal( setdiff(names(pro),  names(pro2)), character(0))  #     fasta
    expect_equal( setdiff(names(fos),  names(fos1)), character(0))  # fos maxquant
    expect_equal( setdiff(names(fos),  names(fos2)), character(0))  #     fasta
})

test_that('`curate` adds anncols',       { 
    expect_setequal( setdiff(names(pro1), names(pro)), anncols)     # pro maxquant
    expect_setequal( setdiff(names(pro2), names(pro)), anncols)     #     fasta
    expect_setequal( setdiff(names(fos1), names(fos)), anncols)     # fos maxquant
    expect_setequal( setdiff(names(fos2), names(fos)), anncols)     #     fasta
})

test_that('`curate` preserves contents', {
    procols <- intersect(names(pro), names(pro1)) %>% setdiff(c('uniprot', 'fastahdrs'))
    foscols <- intersect(names(fos), names(fos1)) %>% setdiff(c('uniprot', 'fastahdrs'))
    expect_equal(pro1[, procols, with = FALSE],  pro[, procols, with = FALSE])   # pro maxquant
    expect_equal(pro2[, procols, with = FALSE],  pro[, procols, with = FALSE])   #     fastafile
    expect_equal(fos1[, foscols, with = FALSE],  fos[, foscols, with = FALSE])   # fos maxquant
    expect_equal(fos2[, foscols, with = FALSE],  fos[, foscols, with = FALSE])   #     fasta
})

test_that('Curated uniprots are a subset', {
    expect_true(
        all(is_collapsed_subset(
                  pro1[reverse=='' & contaminant == '']$uniprot, 
                  pro[ reverse=='' & contaminant == '']$uniprot)))
    expect_true(
        all(is_collapsed_subset( 
                  fos1[reverse=='' & contaminant == '']$uniprot, 
                  fos[ reverse=='' & contaminant == '']$uniprot)))
})


#============================================================================
#
#      add_feature_id
#
#============================================================================

context('`add_feature_id`')
proteinfile <- download_data('billing19.proteingroups.txt')
phosphofile <- download_data('billing19.phosphosites.txt')
pro1 <- .read_maxquant_proteingroups(proteinfile) %>% curate_annotate()
pro2 <- .read_maxquant_proteingroups(proteinfile) %>% curate_annotate(fastadt)
fos1 <- .read_maxquant_phosphosites( phosphofile, proteinfile, quantity = quantity) %>% curate_annotate()
fos2 <- .read_maxquant_phosphosites( phosphofile, proteinfile, quantity = quantity) %>% curate_annotate(fastadt)

pro1b <- add_feature_id(pro1)
pro2b <- add_feature_id(pro2)
fos1b <- add_feature_id(fos1)
fos2b <- add_feature_id(fos2)

test_that('`add_feature_id` preserves rows',     {
    expect_equal( nrow(pro1b), nrow(pro1))  # pro maxquant
    expect_equal( nrow(pro2b), nrow(pro2))  #     fasta
    expect_equal( nrow(fos1b), nrow(fos1))  # fos maxquant
    expect_equal( nrow(fos2b), nrow(fos2))  #     fasta
})

test_that('`add_feature_id` preserves cols',     { 
    expect_equal( setdiff(names(pro1),  names(pro1b)), character(0))    # pro maxquant
    expect_equal( setdiff(names(pro2),  names(pro2b)), character(0))    #     fasta
    expect_equal( setdiff(names(fos1),  names(fos1b)), character(0))    # fos maxquant
    expect_equal( setdiff(names(fos2),  names(fos2b)), character(0))    #     fasta
})

test_that('`add_feature_id` adds `feature_id`',       { 
    expect_setequal( setdiff(names(pro1b), names(pro1)), 'feature_id')  # pro maxquant
    expect_setequal( setdiff(names(pro2b), names(pro2)), 'feature_id')  #     fasta
    expect_setequal( setdiff(names(fos1b), names(fos1)), 'feature_id')  # fos maxquant
    expect_setequal( setdiff(names(fos2b), names(fos2)), 'feature_id')  #     fasta
})

test_that('feature_ids are unique',       {
    expect_true(has_no_duplicates(pro1b$feature_id))
    expect_true(has_no_duplicates(pro2b$feature_id))
    expect_true(has_no_duplicates(pro1b$feature_id))
    expect_true(has_no_duplicates(pro2b$feature_id))
})
