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
#      read_proteingroups
#
#============================================================================

context('read_proteingroups: fukuda20')
metadata <- S4Vectors::metadata

test_that("read_proteingroups(file)", {
    file <- download_data('fukuda20.proteingroups.txt')
    object <- read_proteingroups(file, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'subgroup'    %in% svars(object))
})

test_that("read_proteingroups(file, subgroupvar = NULL) works", {
    file <- download_data('fukuda20.proteingroups.txt')
    object <- read_proteingroups(file, subgroupvar = NULL, plot = FALSE)
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
    expect_true('limma' %in% names(metadata(object)))
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
