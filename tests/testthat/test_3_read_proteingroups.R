#============================================================================
#
#      dequantify
#
#============================================================================

context('dequantify / demultiplex')

# Ratios
    test_that('Ratio', {
        
        expect_identical( dequantify(                           # scalar
              'Ratio H/L t0(L).t1(M).t2(H).R1' ),
                        't0(L).t1(M).t2(H).R1{H/L}' )
        
        expect_identical( dequantify(                           # vector
            c('Ratio H/L t0(L).t1(M).t2(H).R1',
              'Ratio M/L t0(L).t1(M).t2(H).R1')), 
                      c('t0(L).t1(M).t2(H).R1{H/L}',
                        't0(L).t1(M).t2(H).R1{M/L}') )
    })

# Normalized Ratios
    test_that('Ratio normalized', {
        
        expect_identical( dequantify(                           # scalar
              'Ratio H/L normalized t0(L).t1(M).t2(H).R1' ),
                                   't0(L).t1(M).t2(H).R1{H/L}' )
        
        expect_identical( dequantify(                           # vector
            c('Ratio H/L normalized t0(L).t1(M).t2(H).R1',
              'Ratio M/L normalized t0(L).t1(M).t2(H).R1')), 
                           c('t0(L).t1(M).t2(H).R1{H/L}',
                             't0(L).t1(M).t2(H).R1{M/L}') )
    })

# LFQ intensities
    test_that('LFQ intensity', {
        
        expect_identical( dequantify(                           # scalar
            'LFQ intensity t0.R1'),
                          't0.R1' )
        
        expect_identical( dequantify(                           # vector
            c('LFQ intensity t0.R1', 'LFQ intensity t1.R1')),
                          c('t0.R1', 't1.R1'))
    })

# LFQ intensties labeled
    test_that('LFQ intensity labeled', {
        
        expect_identical( dequantify(                           # scalar
            'LFQ intensity L t0(L).t1(H).R1'),
                            't0(L).t1(H).R1{L}' )
        
        expect_identical( dequantify(                           # vector
            c('LFQ intensity L t0(L).t1(H).R1',
              'LFQ intensity H t0(L).t1(H).R1')), 
                            c('t0(L).t1(H).R1{L}', 
                              't0(L).t1(H).R1{H}'))
    })

# Reporter intensities, default label, 0-based
    test_that('Reporter intensity, default label, 0-based', {
        
        expect_identical( dequantify(                           # scalar
              'Reporter intensity 0 t0(0).t1(1).R1'),
                                   't0(0).t1(1).R1{0}' )
        
        expect_identical( dequantify(                           # vector
            c('Reporter intensity 0 t0(0).t1(1).R1',
              'Reporter intensity 1 t0(0).t1(1).R1')),
                                 c('t0(0).t1(1).R1{0}', 
                                   't0(0).t1(1).R1{1}'))
    })

# Reporter intensities, default label, 1-based
    test_that('Reporter intensity, default label, 1-based', {
        
        expect_identical( dequantify(                            # scalar
              'Reporter intensity 1 t0(1).t1(2).R1'),
                                   't0(1).t1(2).R1{1}' )
        
        expect_identical( dequantify(                            # vector
            c('Reporter intensity 1 t0(1).t1(2).R1',  
              'Reporter intensity 2 t0(1).t1(2).R1')),
                                 c('t0(1).t1(2).R1{1}', 
                                   't0(1).t1(2).R1{2}'))
    })

# Reporter intensities, custom label, 0-based
    test_that('Reporter intensity, custom label, 0-based', {
        
        expect_identical( dequantify(                              # scalar
              'Reporter intensity 0 t0(126).t1(127).R1'),
                                   't0(126).t1(127).R1{126}' )
        
        expect_identical( dequantify(                              # vector
            c('Reporter intensity 0 t0(126).t1(127).R1',
              'Reporter intensity 1 t0(126).t1(127).R1')),
                                 c('t0(126).t1(127).R1{126}', 
                                   't0(126).t1(127).R1{127}'))
    })


# Reporter intensities, custom label, 1-based
    test_that('Reporter intensity, custom label, 1-based', {
        
        expect_identical( dequantify(                                # scalar
              'Reporter intensity 1 t0(126).t1(127).R1'),
                                   't0(126).t1(127).R1{126}' )
        
        expect_identical( dequantify(                                # vector
            c('Reporter intensity 1 t0(126).t1(127).R1', 
              'Reporter intensity 2 t0(126).t1(127).R1')),
                                 c('t0(126).t1(127).R1{126}', 
                                   't0(126).t1(127).R1{127}'))
    })


# Corrected reporter intensities, default label, 0-based
    test_that('Reporter intensity corrected, default label, 0-based)', {
        
        expect_identical(dequantify(                                  # scalar
              'Reporter intensity corrected 0 t0(0).t1(1).R1'),
                                             't0(0).t1(1).R1{0}' )
        
        expect_identical(dequantify(                                  # vector
            c('Reporter intensity corrected 0 t0(0).t1(1).R1',
              'Reporter intensity corrected 1 t0(0).t1(1).R1')),
           c('t0(0).t1(1).R1{0}', 
             't0(0).t1(1).R1{1}'))
    })

# Corrected reporter intensities, default label, 1-based
    test_that('Reporter intensity corrected, default label, 1-based)', {
        
        expect_identical( dequantify(                                 # scalar
              'Reporter intensity corrected 1 t0(1).t1(2).R1'),
                                             't0(1).t1(2).R1{1}' )
        
        expect_identical( dequantify(                                 # vector
            c('Reporter intensity corrected 1 t0(1).t1(2).R1',
              'Reporter intensity corrected 2 t0(1).t1(2).R1')), 
                                           c('t0(1).t1(2).R1{1}', 
                                             't0(1).t1(2).R1{2}'))
    })

# Corrected reporter intensities, custom label, 0-based
test_that('Reporter intensity corrected, custom label, 0-based', {

    expect_identical( dequantify(
          'Reporter intensity corrected 0 t0(126).t1(127).R1'),
                                         't0(126).t1(127).R1{126}')    # scalar
    
    expect_identical( dequantify(
        c('Reporter intensity corrected 0 t0(126).t1(127).R1',         # vector
          'Reporter intensity corrected 1 t0(126).t1(127).R1')), 
                                       c('t0(126).t1(127).R1{126}', 
                                         't0(126).t1(127).R1{127}'))
})
  
# Corrected reporter intensities, custom label, 1-based
test_that('Reporter intensity corrected, custom label, 1-based', {
    
    expect_identical( dequantify(
          'Reporter intensity corrected 1 t0(126).t1(127).R1'),       # scalar
                                         't0(126).t1(127).R1{126}' )
    
    expect_identical( dequantify(                                      # vector
        c('Reporter intensity corrected 1 t0(126).t1(127).R1', 
          'Reporter intensity corrected 2 t0(126).t1(127).R1')), 
                                       c('t0(126).t1(127).R1{126}', 
                                         't0(126).t1(127).R1{127}'))
})


#============================================================================
#
#      demultiplex
#
#============================================================================

context('demultiplex')


test_that('Ratio', {
    
    expect_identical( demultiplex('t0(L).t1(M).t2(H).R1{H/L}' ), 't2_t0.R1' )
    
    expect_identical( demultiplex(
        c('t0(L).t1(M).t2(H).R1{H/L}', 't0(L).t1(M).t2(H).R1{M/L}')), 
        c('t2_t0.R1', 't1_t0.R1') )
    
})

test_that('LFQ intensity', {
    
        expect_identical( demultiplex('t0.R1'),'t0.R1' )
        
        expect_identical( demultiplex(c('t0.R1', 't1.R1')), c('t0.R1', 't1.R1'))
        
        expect_identical( demultiplex('t0(L).t1(H).R1{L}'), 't0.R1' )
    
        expect_identical( demultiplex(
            c('t0(L).t1(H).R1{L}', 't0(L).t1(H).R1{H}')), 
            c('t0.R1',             't1.R1'))
})


test_that('Reporter intensity', {

    expect_identical( demultiplex('t0(126).t1(127).R1{126}'), 't0.R1')
    
    expect_identical( demultiplex(
        c('t0(126).t1(127).R1{126}', 't0(126).t1(127).R1{127}')), 
        c('t0.R1', 't1.R1'))
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
