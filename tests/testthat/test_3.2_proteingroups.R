#============================================================================
#
#      .read_proteingroups
#      .read_phosphosites
#
#============================================================================

context('.read_proteingroups/.read_phosphosites')
proteinfile <- download_data('billing19.proteingroups.txt')
phosphofile <- download_data('billing19.phosphosites.txt')
pro <- .read_proteingroups(proteinfile, verbose = TRUE)
fos <- .read_phosphosites(phosphofile, proteinfile, verbose = TRUE)
prodt <- fread(proteinfile, colClasses = c(id='character'))
fosdt <- fread(phosphofile, colClasses = c(id='character'), integer64='numeric')
fosdt %<>% extract(fos$fosId, on = 'id')

test_that('id values match', {
    expect_identical(pro$proId, prodt$id)
    expect_identical(fos$fosId, fosdt$id)
})

test_that('Uniprot values match', {
    expect_identical(pro$Uniprot, prodt$`Majority protein IDs`)
    expect_identical(fos$Uniprot, fosdt$Proteins)
})

test_that('peptide counts match', {
    expect_identical(   pro$`Razor + unique peptides STD(L).E00(M).E01(H).R1`, 
                      prodt$`Razor + unique peptides STD(L).E00(M).E01(H).R1`)
    expect_identical(   fos$`Razor + unique peptides STD(L).E00(M).E01(H).R1`, 
                      fosdt$`Razor + unique peptides STD(L).E00(M).E01(H).R1`)
})

test_that('normalized ratios match', {
    expect_identical(pro$`Ratio H/M normalized STD(L).E02(M).E05(H).R8`, 
                   prodt$`Ratio H/M normalized STD(L).E02(M).E05(H).R8`   )
    expect_identical(fos$`Ratio H/M normalized STD(L).E02(M).E05(H).R8`, 
                   fosdt$`Ratio H/M normalized STD(L).E02(M).E05(H).R8___1`)
})


#============================================================================
#
#      drop_differing_uniprots
#
#============================================================================

context('`drop_differing_uniprots`')
proteinfile <- download_data('billing19.proteingroups.txt')
phosphofile <- download_data('billing19.phosphosites.txt')
prodt <- .read_proteingroups(proteinfile,              verbose = TRUE)
fosdt <- .read_phosphosites( phosphofile, proteinfile, verbose = TRUE)
fosdt1 <- drop_differing_uniprots(fosdt, prodt,        verbose = TRUE)

test_that('`drop_differing_uniprots` preserves colnames', {
    expect_setequal(names(fosdt), names(fosdt1))
})

test_that('`drop_differing_uniprots` preserves contents', {
    fosdt1 %<>% extract(, names(fosdt), with = FALSE)
    cols <- setdiff(names(fosdt), c('Uniprot', 'Positions within proteins'))
    expect_identical(fosdt[ , cols, with = FALSE], 
                     fosdt1[, cols, with = FALSE] )
})

test_that('`drop_differing_uniprots` uniprots are subset of original', {
    usplit <- function(x)  unlist(stri_split_fixed(x, ';'))
    is_string_subset <- function(x, y)  is_subset(usplit(x), usplit(y))
    expect_true(is_string_subset(fosdt1$Uniprot[   1], fosdt$Uniprot[   1]))
    expect_true(is_string_subset(fosdt1$Uniprot[  10], fosdt$Uniprot[  10]))
    expect_true(is_string_subset(fosdt1$Uniprot[ 100], fosdt$Uniprot[ 100]))
})

test_that('`drop_differing_uniprots` maintains integrity between `Positions within proteins` and `Uniprot`', {
    fosdt  %<>% extract( , c('fosId', 'Uniprot', 'Positions within proteins'), with = FALSE)
    fosdt1 %<>% extract( , c('fosId', 'Uniprot', 'Positions within proteins'), with = FALSE)
    fosdt  %<>% uncollapse(Uniprot, `Positions within proteins`)
    fosdt1 %<>% uncollapse(Uniprot, `Positions within proteins`)
    fosdt  %<>% merge(fosdt1, by = c('fosId', 'Uniprot'))
    expect_identical(fosdt$`Positions within proteins.x`, 
                     fosdt$`Positions within proteins.y`)
})


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
        fastadt['P31946', on = 'Uniprot'],
        data.table(
            Reviewed  = 1,
            Entry     = '1433B', 
            Gene      = 'YWHAB',
            Uniprot   = 'P31946', 
            Canonical = 'P31946', 
            Isoform   = 1,
            Protein   = '14-3-3 protein beta/alpha',
            Fragment  = 0,
            Existence = 1)))

test_that('`read_fastahdrs` reads intermediate swissprot protein',
    expect_equal(
        fastadt['Q9BUJ2-4', on = 'Uniprot'],
        data.table(
            Reviewed  = 1,
            Entry     = 'HNRL1', 
            Gene      = 'HNRNPUL1',
            Uniprot   = 'Q9BUJ2-4', 
            Canonical = 'Q9BUJ2', 
            Isoform   = 4,
            Protein   = 'Isoform 4 of Heterogeneous nuclear ribonucleoprotein U-like protein 1', 
            Fragment  = 0,
            Existence = 1)))

test_that('`read_fastahdrs` reads intermediate trembl protein',
    expect_equal(
        fastadt['G5E9N3', on = 'Uniprot'],
        data.table(
            Reviewed  = 0,
            Entry     = 'G5E9N3', 
            Gene      = 'RETSAT',
            Uniprot   = 'G5E9N3', 
            Canonical = 'G5E9N3', 
            Isoform   = 1,
            Protein   = 'All-trans-13,14-dihydroretinol saturase, isoform CRA_c', 
            Fragment  = 0,
            Existence = 4)))

test_that('`read_fastahdrs` reads last trembl protein',
    expect_equal(
        fastadt['R4GMM2', on = 'Uniprot'],
        data.table(
            Reviewed  = 0,
            Entry     = 'R4GMM2', 
            Gene      = 'PARD6A',
            Uniprot   = 'R4GMM2', 
            Canonical = 'R4GMM2', 
            Isoform   = 1,
            Protein   = 'Partitioning defective 6 homolog alpha',
            Fragment  = 0,
            Existence = 4)))


#============================================================================
#
#       maxquant_curate
#          fasta_curate
#
#============================================================================

context('`maxquant_curate/fasta_curate`')
proteinfile <- download_data('billing19.proteingroups.txt')
phosphofile <- download_data('billing19.phosphosites.txt')
pro <- .read_proteingroups(proteinfile, verbose = TRUE)
fos <- .read_phosphosites(phosphofile, proteinfile, verbose = TRUE)
pro1 <- curate_annotate_maxquant(prodt)
fos1 <- curate_annotate_maxquant(fosdt)
pro2 <- curate_annotate_fastafile(prodt, fastadt)
fos2 <- curate_annotate_fastafile(fosdt, fastadt)
anncols <- c("Entry", "Gene","Canonical", "Isoform", "Protein")

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
    procols <- intersect(names(pro), names(pro1)) %>% setdiff(c('Uniprot', 'Fasta headers'))
    foscols <- intersect(names(fos), names(fos1)) %>% setdiff(c('Uniprot', 'Fasta headers'))
    expect_equal(pro1[, procols, with = FALSE],  pro[, procols, with = FALSE])   # pro maxquant
    expect_equal(pro2[, procols, with = FALSE],  pro[, procols, with = FALSE])   #     fastafile
    expect_equal(fos1[, foscols, with = FALSE],  fos[, foscols, with = FALSE])   # fos maxquant
    expect_equal(fos2[, foscols, with = FALSE],  fos[, foscols, with = FALSE])   #     fasta
})

test_that('Curated uniprots are a subset', {
    assert_all_are_true(
        is_collapsed_subset, 
                  pro1[Reverse=='' & `Potential contaminant` == '']$Uniprot, 
                  pro[ Reverse=='' & `Potential contaminant` == '']$Uniprot)
    assert_all_are_true(
        is_collapsed_subset, 
                  fos1[Reverse=='' & `Potential contaminant` == '']$Uniprot, 
                  fos[ Reverse=='' & `Potential contaminant` == '']$Uniprot)
})


#============================================================================
#
#      add_feature_id
#
#============================================================================

context('`add_feature_id`')
fastafile   <- download_data('uniprot_hsa_20140515.fasta')
proteinfile <- download_data('billing19.proteingroups.txt')
phosphofile <- download_data('billing19.phosphosites.txt')
fastadt <- read_fastahdrs(fastafile)
pro1 <- .read_proteingroups(proteinfile) %>% curate_annotate()
pro2 <- .read_proteingroups(proteinfile) %>% curate_annotate(fastadt)
fos1 <- .read_phosphosites(phosphofile, proteinfile) %>% curate_annotate()
fos2 <- .read_phosphosites(phosphofile, proteinfile) %>% curate_annotate(fastadt)

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


#============================================================================
#
#     mqdt_to_mat
#
#============================================================================

context('mqdt_to_mat')

proteinfile <- download_data('billing19.proteingroups.txt')
phosphofile <- download_data('billing19.phosphosites.txt')
prodt <- .read_proteingroups(proteinfile = proteinfile)
prodt %<>% curate_annotate()
prodt %<>% add_feature_id()
quantity <- guess_maxquant_quantity(proteinfile)
pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
mat <- mqdt_to_mat(prodt, pattern = pattern)

dt <- fread(proteinfile, integer64 = 'numeric', colClasses = c(id = 'character'))
dt %<>% extract(, colnames(mat), with = FALSE)
dt %<>% data.matrix()
mat <- 2^mat - 1
rownames(mat) <- NULL
test_that('`mqdt_to_mat` preserves contents', expect_equal(dt, mat))


#============================================================================
#
#      dequantify
#
#============================================================================

context('dequantify')

# Ratios
    test_that('`dequantify` works for (normalized) ratios', {
        expect_identical( dequantify(                       # scalar
              'Ratio H/L WT(L).KD(M).OE(H).R1' ),
                        'WT(L).KD(M).OE(H).R1{H/L}' )
        expect_identical( dequantify(                       # vector
            c('Ratio H/L WT(L).KD(M).OE(H).R1',
              'Ratio M/L WT(L).KD(M).OE(H).R1')), 
            c('WT(L).KD(M).OE(H).R1{H/L}',
              'WT(L).KD(M).OE(H).R1{M/L}') )
        expect_identical( dequantify(                       # replicate (not run)
            c('Ratio H/L WT.R1(L).KD.R1(M).OE.R1(H)',
              'Ratio M/L WT.R1(L).KD.R1(M).OE.R1(H)')), 
                      c('WT.R1(L).KD.R1(M).OE.R1(H){H/L}',
                        'WT.R1(L).KD.R1(M).OE.R1(H){M/L}') )
        expect_identical( dequantify(                       # normalized ratios
            c('Ratio H/L normalized WT(L).KD(M).OE(H).R1',
              'Ratio M/L normalized WT(L).KD(M).OE(H).R1')), 
                           c('WT(L).KD(M).OE(H).R1{H/L}',
                             'WT(L).KD(M).OE(H).R1{M/L}') )
    })

# LFQ intensities
    test_that('`dequantify` works for (labeled) LFQ intensities', {
        expect_identical( dequantify(                           # vector
            c('LFQ intensity WT.R1', 'LFQ intensity KD.R1')),
                          c('WT.R1', 'KD.R1'))
        expect_identical( dequantify(                           # scalar
            'LFQ intensity WT.R1'),
                          'WT.R1' )
        expect_identical( dequantify(                           # labeled LFQs
            c('LFQ intensity L WT(L).KD(H).R1',
              'LFQ intensity H WT(L).KD(H).R1')), 
                            c('WT(L).KD(H).R1{L}', 
                              'WT(L).KD(H).R1{H}'))
    })

# Reporter intensities
    test_that('`dequantify` works for reporter intensities', {
        expect_identical( dequantify(                             # scalar
              'Reporter intensity 0 WT(0).KD(1).R1'),
                                   'WT(1).KD(2).R1{1}' )
        expect_identical( dequantify(                             # vector
            c('Reporter intensity 0 WT(0).KD(1).R1',
              'Reporter intensity 1 WT(0).KD(1).R1')),
                                 c('WT(1).KD(2).R1{1}', 
                                   'WT(1).KD(2).R1{2}'))
        expect_identical( dequantify(                             # 1-based
            c('Reporter intensity 1 WT(1).KD(2).R1',  
              'Reporter intensity 2 WT(1).KD(2).R1')),
                                 c('WT(1).KD(2).R1{1}', 
                                   'WT(1).KD(2).R1{2}'))
        expect_identical( dequantify(                             # label-based
            c('Reporter intensity 1 WT(126).KD(127).R1',
              'Reporter intensity 2 WT(126).KD(127).R1')),
                                 c('WT(1).KD(2).R1{1}', 
                                   'WT(1).KD(2).R1{2}'))
        expect_identical(dequantify(                              # corrected
            c('Reporter intensity corrected 1 WT(1).KD(2).R1',
              'Reporter intensity corrected 2 WT(1).KD(2).R1')),
           c('WT(1).KD(2).R1{1}', 
             'WT(1).KD(2).R1{2}'))
    })


#============================================================================
#
#      demultiplex
#
#============================================================================

context('demultiplex')


test_that('`demultiplex` works', {
    expect_identical( demultiplex(                                # scalar
        'WT(1).KD(2).R1{1}'), 
        'WT.R1')
    expect_identical( demultiplex(                                # uniplexed
        c('WT.R1', 'KD.R1')), 
        c('WT.R1', 'KD.R1'))
    expect_identical( demultiplex(
        c('WT(1).KD(2).R1{1}', 'WT(1).KD(2).R1{2}')),             # multiplexed
        c('WT.R1', 'KD.R1'))
    expect_identical( demultiplex(
        c('WT(L).KD(H).R1{L}', 'WT(L).KD(H).R1{H}')),             # labels
        c('WT.R1', 'KD.R1'))
    expect_identical( demultiplex(
        c('WT.R1(L).KD.R1(H){L}', 'WT.R1(L).KD.R1(H){H}')),       # replicates
        c('WT.R1', 'KD.R1'))
    expect_identical( demultiplex(
        c('WT(L).KD(M).OE(H).R1{H/L}', 'WT(L).KD(M).OE(H).R1{M/L}')), # ratios
        c('OE_WT.R1', 'KD_WT.R1') )
})


#============================================================================
#
#      read_phosphosites
#
#============================================================================

proteinfile <- download_data('billing19.proteingroups.txt')
phosphofile <- download_data('billing19.phosphosites.txt')
quantity <- guess_maxquant_quantity(proteinfile)
pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
fos <- read_phosphosites(phosphofile, proteinfile, plot = FALSE, impute = FALSE)

str(log2sites(fos))
dt <- fread(phosphofile, integer64 = 'numeric', colClasses = c(id = 'character'))
dt %<>% extract(fdt(fos)$fosId, on = 'id')
dt <- dt[, .SD, .SDcols = patterns(pattern, '___1')]
colnames(dt) %<>% stri_replace_first_fixed('___1', '')
colnames(dt) %<>% dequantify()
colnames(dt) %<>% demultiplex()

test_that('`read_phosphosites` preserves samples', 
    expect_equal(snames(fos), colnames(dt))
)

test_that('`read_phosphosites` preserves phosphosite values', {
    str(data.matrix(dt))
    str(log2sites(fos))
    expect_equal(snames(fos), colnames(dt))
})


#============================================================================
#
#      read_proteingroups
#
#============================================================================

proteinfile <- download_data('billing19.proteingroups.txt')
quantity <- guess_maxquant_quantity(proteinfile)
pattern <- MAXQUANT_PATTERNS_QUANTITY[[quantity]]
pro <- read_proteingroups(proteinfile, plot = FALSE, impute = FALSE)
mat <- log2proteins(pro)
rownames(mat) <- fdt(pro)$proId

test_that('`read_proteingroups(billing19)` returns SummarizedExperiment', {
    expect_s4_class(pro, 'SummarizedExperiment')
})

test_that('`read_proteingroups(billing19)` reads abundances', {
    dt <- fread(proteinfile)
    mat0 <- mqdt_to_mat(dt, pattern)
    rownames(mat0) <- dt$id
    colnames(mat0) %<>% dequantify()
    colnames(mat0) %<>% demultiplex()
    mat0 %<>% extract(rownames(mat), )
    expect_equal(mat, mat0)
})

test_that('`read_proteingroups(billing19)` reads uniprots', {
    fdt0 <- fread(proteinfile, select = c('id', 'Majority protein IDs'), colClasses = c(id = 'character'))
    fdt0 %<>% extract(fdt(pro)$proId, on = 'id')
    assert_all_are_true(is_collapsed_subset(fdt(pro)$Uniprot, fdt0$`Majority protein IDs`))
})


context('read_proteingroups: fukuda20')
metadata <- S4Vectors::metadata

test_that("read_proteingroups(file)", {
    file <- download_data('fukuda20.proteingroups.txt')
    object <- read_proteingroups(file, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'subgroup'    %in% svars(object))
})

test_that(  "read_proteingroups(file, pca=TRUE)", {
    file <- download_data('fukuda20.proteingroups.txt')
    object <- read_proteingroups(file, pca = TRUE, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(all(c('pca1', 'pca2') %in% svars(object)))
    expect_true(all(c('pca1', 'pca2') %in% fvars(object)))
    expect_true('pca' %in% names(metadata(object)))
})

test_that(  "read_proteingroups(file, fit='limma')", {
    file <- download_data('fukuda20.proteingroups.txt')
    object <- read_proteingroups(file, fit = 'limma', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0('fdr', FITSEP))))
})

# test_that(  "read_proteingroups(file, fit='lm')", {
#     file <- download_data('fukuda20.proteingroups.txt')
#     object <- read_proteingroups(file, fit = 'lm', plot = FALSE)
#     expect_s4_class(object, 'SummarizedExperiment')
#     expect_true('lm' %in% names(metadata(object)))
# })

test_that(  "read_proteingroups(file, fit='wilcoxon')", {
    file <- download_data('fukuda20.proteingroups.txt')
    object <- read_proteingroups(file, fit = 'wilcoxon', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), '~wilcoxon')))
    expect_true('wilcoxon' %in% names(metadata(object)))
})
