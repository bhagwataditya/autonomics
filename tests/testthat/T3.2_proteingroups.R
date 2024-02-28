require(testthat)


#============================================================================
#                                                                           #
             context('.read_(proteingroups|phosphosites)')                  #
#                                                                           #
#============================================================================


test_that(" .read_(proteingroups|phosphosites) ", {
    # read
        proteinfile <- download_data('billing19.proteingroups.txt')
        phosphofile <- download_data('billing19.phosphosites.txt')
        pro <- .read_maxquant_proteingroups(proteinfile, verbose = TRUE)
        fos <- .read_maxquant_phosphosites(phosphofile, proteinfile, verbose = TRUE)
        prodt <- fread(proteinfile, colClasses = c(id = 'character'))
        fosdt <- fread(phosphofile, colClasses = c(id = 'character'), integer64 = 'numeric')
        fosdt %<>% extract(fos$fosId, on = 'id')
    # ids
        expect_identical(pro$proId, prodt$id)
        expect_identical(fos$fosId, fosdt$id)
    # uniprots
        expect_identical(pro$uniprot, prodt$`Majority protein IDs`)
        expect_identical(fos$uniprot, fosdt$Proteins)
    # pecounts
        expect_identical(   pro$`Razor + unique peptides STD(L).E00(M).E01(H).R1`, 
                          prodt$`Razor + unique peptides STD(L).E00(M).E01(H).R1`)
        expect_identical(   fos$`Razor + unique peptides STD(L).E00(M).E01(H).R1`, 
                          fosdt$`Razor + unique peptides STD(L).E00(M).E01(H).R1`)
    # normalized ratios
        expect_identical(pro$`Ratio H/M normalized STD(L).E02(M).E05(H).R8`, 
                       prodt$`Ratio H/M normalized STD(L).E02(M).E05(H).R8`   )
        expect_identical(fos$`Ratio H/M normalized STD(L).E02(M).E05(H).R8`, 
                       fosdt$`Ratio H/M normalized STD(L).E02(M).E05(H).R8___1`)
})



#============================================================================
#                                                                           #
            context('drop_differing_uniprots')                              #
#                                                                           #
#============================================================================


test_that(" `drop_differing_uniprots` ", {
    # read
        proteinfile <- download_data('billing19.proteingroups.txt')
        phosphofile <- download_data('billing19.phosphosites.txt')
        prodt <- .read_maxquant_proteingroups(proteinfile)
        fosdt <- .read_maxquant_phosphosites( phosphofile, proteinfile)
        fosdt1 <- drop_differing_uniprots(fosdt, prodt, verbose = TRUE)
    # colnames
        expect_setequal(names(fosdt), names(fosdt1))
    # contents preserved
        fosdt1 %<>% extract(, names(fosdt), with = FALSE)
        cols <- setdiff(names(fosdt), c('uniprot', 'Positions within proteins'))
        expect_identical(fosdt[ , cols, with = FALSE], 
                         fosdt1[, cols, with = FALSE] )
    # uniprots are subset of original
        usplit <- function(x)  unlist(stri_split_fixed(x, ';'))
        is_string_subset <- function(x, y)  is_subset(usplit(x), usplit(y))
        expect_true(is_string_subset(fosdt1$uniprot[   1], fosdt$uniprot[   1]))
        expect_true(is_string_subset(fosdt1$uniprot[  10], fosdt$uniprot[  10]))
        expect_true(is_string_subset(fosdt1$uniprot[ 100], fosdt$uniprot[ 100]))
    # `Positions within proteins` and `uniprot` have same order
        fosdt  %<>% extract( , c('fosId', 'uniprot', 'Positions within proteins'), with = FALSE)
        fosdt1 %<>% extract( , c('fosId', 'uniprot', 'Positions within proteins'), with = FALSE)
        fosdt  %<>% uncollapse(uniprot, `Positions within proteins`)
        fosdt1 %<>% uncollapse(uniprot, `Positions within proteins`)
        fosdt  %<>% merge(fosdt1, by = c('fosId', 'uniprot'))
        expect_identical(fosdt$`Positions within proteins.x`, 
                         fosdt$`Positions within proteins.y`)
})



#============================================================================
#                                                                           #
            context('mqdt_to_mat')                                          #
#                                                                           #
#============================================================================


test_that(" `mqdt_to_mat` ", {
    # Read
        proteinfile <- download_data('billing19.proteingroups.txt')
        phosphofile <- download_data('billing19.phosphosites.txt')
        prodt <- .read_maxquant_proteingroups(proteinfile)
        uniprothdrs <- NULL
        contaminanthdrs <- read_contaminantdt()
        maxquanthdrs <- parse_maxquant_hdrs(prodt$`Fasta headers`);   prodt[, `Fasta headers` := NULL ]
        prodt %<>% annotate_maxquant(uniprothdrs = uniprothdrs, contaminanthdrs = contaminanthdrs, maxquanthdrs = maxquanthdrs, restapi = FALSE)
        quantity <- guess_maxquant_quantity(proteinfile)
        pattern <- MAXQUANT_PATTERNS[[quantity]]
        mat <- 2^mqdt_to_mat(prodt, pattern = pattern)
        
        prodt %<>% extract(, c('feature_id', colnames(mat)), with = FALSE)
        prodt %<>% dt2mat()
        prodt[is.nan(prodt)] <- NA
    # Preserves contents
        expect_equal(prodt, mat)
})



#============================================================================
#                                                                           #
               context('dequantify')                                        #
#                                                                           #
#============================================================================


# Ratios
    test_that(" dequantify : (normalized) ratios ", {
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
    test_that(" dequantify : (labeled) LFQ intensities ", {
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
    test_that(" dequantify : reporter intensities ", {
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
#                                                                           #
               context('demultiplex')                                       #
#                                                                           #
#============================================================================


test_that(" demultiplex ", {
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
#                                                                           #
             context('read_proteingroups : integer64')                      #
#                                                                           #
#============================================================================


test_that(" read_proteingroups: integer64 ", {
            file <- download_data('integer64.proteinGroups.txt')
            object <- .read_maxquant_proteingroups(file)
            expect_false( 'integer64' %in% sapply(object, class))
})


#============================================================================
#                                                                           #
             context('read_proteingroups : fukuda20')                       #
#                                                                           #
#============================================================================

test_that( " read_proteingroups: fukuda20 ", {
    
            file <- download_data('fukuda20.proteingroups.txt')
            object <- read_maxquant_proteingroups(file)
            expect_s4_class(object, 'SummarizedExperiment')
            expect_true( 'subgroup'    %in% svars(object))
})

    
test_that( " read_proteingroups: fukuda20, pca = TRUE ", {
    
            file <- download_data('fukuda20.proteingroups.txt')
            object <- read_maxquant_proteingroups(file, pca = TRUE)
            expect_s4_class(object, 'SummarizedExperiment')
            expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% svars(object)))
            expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% fvars(object)))
            expect_true('sample_id~pca' %in% names(metadata(object)))
})


test_that( " read_proteingroups: fukuda20, fit = 'limma' ", {
    
    file <- download_data('fukuda20.proteingroups.txt')
    object <- read_maxquant_proteingroups(file, fit = 'limma', plot = TRUE, label = NULL)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'limma'))))
    expect_true(any(stri_detect_fixed(fvars(object), paste0('fdr', FITSEP))))
})


test_that( " read_proteingroups: fukuda20, fit = 'lm' ", {
    
    file <- download_data('fukuda20.proteingroups.txt')
    object <- read_maxquant_proteingroups(file, fit = 'lm', plot = TRUE, label = NULL)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP,  'lm'))))
})


test_that( " read_proteingroups: fukuda20, fit = 'wilcoxon' ", {
    
    file <- download_data('fukuda20.proteingroups.txt')
    object <- read_maxquant_proteingroups(file, fit = 'wilcoxon')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'wilcoxon'))))
})



#============================================================================
#                                                                           #
             context(" read_proteingroups: billing16 ")                     #
#                                                                           #
#============================================================================


test_that(" read_proteingroups: billing16 ", {
        file <- download_data('billing16.proteingroups.txt')
        object <- read_maxquant_proteingroups(file)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true( 'subgroup'    %in% svars(object))
})


test_that(" read_proteingroups: billing16 , quantity = 'labeledintensity' ", {
        file <- download_data('billing16.proteingroups.txt')
        object <- read_maxquant_proteingroups(file, quantity = 'labeledintensity')
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('log2labeledintensity' %in% assayNames(object))
})

test_that(" read_proteingroups: billing16, invert = selection ", {
    # read_proteingroups(invert = TRUE)
        file <- download_data('billing16.proteingroups.txt')
        invert <- c('EM_E', 'BM_E', 'BM_EM')
        object <- read_maxquant_proteingroups(file, invert = invert)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_setequal(levels(object$subgroup), c('E_BM','E_EM','EM_BM'))
        expect_setequal(object$sample_id,  c("E_BM.R1",  "E_BM.R2",  "E_BM.R3",  "E_EM.R1", 
                                  "E_EM.R2", "E_EM.R3", "EM_BM.R1", "EM_BM.R2", "EM_BM.R3")) 
    # read_proteingroups %>% invert()
        object <- read_maxquant_proteingroups(file)
        object %<>% invert_subgroups(invert) 
        expect_setequal(levels(object$subgroup), c('E_BM','E_EM','EM_BM'))
        expect_setequal(object$sample_id,  c("E_BM.R1",  "E_BM.R2",  "E_BM.R3",  "E_EM.R1", 
                                  "E_EM.R2", "E_EM.R3", "EM_BM.R1", "EM_BM.R2", "EM_BM.R3")) 
})


if (requireNamespace('Biostrings', quietly = TRUE)){
test_that(" read_proteingroups: billing16, fastafile = uniprot_hsa_20140515.fasta", {
    # Read
        file <- download_data('billing16.proteingroups.txt')
        fastafile <- download_data('uniprot_hsa_20140515.fasta')
        invert <- c('EM_E', 'BM_E', 'BM_EM')
        object <- read_maxquant_proteingroups(file, invert = invert, fastafile = fastafile)
        dt <- .read_maxquant_proteingroups(file)
        dt %<>% extract(match(fdt(object)$proId, proId))
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('protein' %in% fvars(object))
        expect_true(any(nchar(fdt(object)$uniprot) < nchar(dt$uniprot)))
})}


test_that(" read_proteingroups: file, pca = TRUE ", {
    # Read
        file <- download_data('billing16.proteingroups.txt')
        invert <- c('EM_E', 'BM_E', 'EM_BM')
        object <- read_maxquant_proteingroups(file, invert = invert, pca = TRUE)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(         svars(object),  'pca1')))
        expect_true(any(stri_detect_fixed(         fvars(object),  'pca1')))
        expect_true(any(stri_detect_fixed(names(metadata(object)), 'pca')))
})


test_that(" read_proteingroups: billing16, impute = TRUE ", {
    # Read
        file <- download_data('billing16.proteingroups.txt')
        invert <- c('EM_E', 'BM_E', 'EM_BM')
        object <- read_maxquant_proteingroups(file, invert = invert, impute = TRUE)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('is_imputed' %in% SummarizedExperiment::assayNames(object))
})


test_that(" read_maxquant_proteingroups: file, fit = 'limma' ", {
    # Read
        file <- download_data('billing16.proteingroups.txt')
        invert <- c('EM_E', 'BM_E', 'EM_BM')
        formula <- ~0+subgroup
        object <- read_maxquant_proteingroups(file, invert = invert, 
                          fit = 'limma', formula = formula, coefs = c('BM_EM', 'E_EM', 'E_BM'))
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), 'limma')))
})


test_that(" read_maxquant_proteingroups: file, fit = 'lm' ", {
    # Read
        file <- download_data('billing16.proteingroups.txt')
        invert <- c('EM_E', 'BM_E', 'EM_BM')
        formula <- ~0 + subgroup
        object <- read_maxquant_proteingroups(file, invert = invert, fit = 'lm', formula = formula)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), 'lm')))
})


test_that(" read_maxquant_proteingroups: file, fit = 'wilcoxon' ", {
    # Read
        file <- download_data('billing16.proteingroups.txt')
        invert <- c('EM_E', 'BM_E', 'EM_BM')
        object <- read_maxquant_proteingroups(file, invert = invert, 
                                              fit = 'wilcoxon', formula = formula)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), 'wilcoxon')))
})




#============================================================================
#                                                                           #
             context(" read_proteingroups: billing19 ")                     #
#                                                                           #
#============================================================================

test_that(" read_proteingroups: billing19 ", {
    # SummarizedExperiment
        file <- download_data('billing19.proteingroups.txt')
        object <- read_maxquant_proteingroups(file = file)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('subgroup' %in% svars(object))
    # values
        dt <- fread(proteinfile)
        mat0 <- mqdt_to_mat(dt, pattern)
        rownames(mat0) <- dt$id
        colnames(mat0) %<>% dequantify()
        colnames(mat0) %<>% demultiplex()
        mat <- values(object)
        rownames(mat) <- fdt(object)$proId
        mat0 %<>% extract(rownames(mat), )
        expect_equal(mat, mat0)
    # uniprots
        fdt0 <- fread(proteinfile,   select = c('id', 'Majority protein IDs'), 
                                 colClasses = c( id = 'character'))
        fdt0 %<>% extract(fdt(object)$proId, on = 'id')
        expect_true(all(is_collapsed_subset(fdt(object)$uniprot, fdt0$`Majority protein IDs`)))
})


test_that(" read_proteingroups: billing19 , subgroups ", {
    # Read
        subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
        file <- download_data('billing19.proteingroups.txt')
        object <- read_maxquant_proteingroups(file, subgroups = subgroups)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(assertive.sets::are_set_equal( slevels(object, 'subgroup'), subgroups))
})


test_that(" read_proteingroups: billing19 , impute = TRUE ", {
    # Read
        subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
        file <- download_data('billing19.proteingroups.txt')
        object <- read_maxquant_proteingroups(file, subgroups = subgroups, impute = TRUE)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('is_imputed' %in% SummarizedExperiment::assayNames(object))
})


test_that(" read_proteingroups: billing19 , pca = TRUE ", {
    # Read
        subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
        file <- download_data('billing19.proteingroups.txt')
        object <- read_maxquant_proteingroups(file, subgroups = subgroups, pca = TRUE)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(         svars(object),  'pca1')))
        expect_true(any(stri_detect_fixed(         fvars(object),  'pca1')))
        expect_true(any(stri_detect_fixed(names(metadata(object)), 'pca' )))
})

test_that(" read_proteingroups: billing19 , fit = 'limma' ", {
    # Read
        subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
        file <- download_data('billing19.proteingroups.txt')
        object <- read_maxquant_proteingroups(file, subgroups = subgroups, fit = 'limma')
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), 'limma')))
})

test_that(" read_proteingroups: billing19 , fit = 'lm' ", {
    # Read
        subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
        file <- download_data('billing19.proteingroups.txt')
        object <- read_maxquant_proteingroups(file, subgroups = subgroups, fit = 'lm')
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), 'lm')))
})

test_that( " read_proteingroups: billing19 , fit = 'wilcoxon' ", {
    # Read
        subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
        file <- download_data('billing19.proteingroups.txt')
        object <- read_maxquant_proteingroups(file, subgroups = subgroups, fit = 'wilcoxon')
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), 'wilcoxon')))
})


#============================================================================
#                                                                           #
             context('read_phosphosites : billing19')                       #
#                                                                           #
#============================================================================



test_that(" read_phosphosites: billing19 ", {
    # sumexp
        subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
        proteinfile <- download_data('billing19.proteingroups.txt')
        phosphofile <- download_data('billing19.phosphosites.txt')
        object <- read_maxquant_phosphosites(
            phosphofile = phosphofile, proteinfile = proteinfile, subgroups = subgroups)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('log2proteins' %in% SummarizedExperiment::assayNames(object))
        expect_true('log2sites'    %in% SummarizedExperiment::assayNames(object))
        expect_true('log2diffs'    %in% SummarizedExperiment::assayNames(object))
    # snames
        dt <- fread(phosphofile, integer64 = 'numeric', colClasses = c(id = 'character'))
        dt %<>% extract(fdt(fos)$fosId, on = 'id')
        dt <- dt[, .SD, .SDcols = patterns(pattern, '___1')]
        colnames(dt) %<>% stri_replace_first_fixed('___1', '')
        colnames(dt) %<>% dequantify()
        colnames(dt) %<>% demultiplex()
        idx <- split_extract_fixed(colnames(dt), '.', 1) %in% subgroups
        dt %<>% extract(, idx, with = FALSE)
        expect_equal(snames(object), colnames(dt))
})



