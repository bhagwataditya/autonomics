require(testthat)
metadata <- S4Vectors::metadata

#============================================================================
#                                                                           #
             context(" read_metabolon ")                                    #
#                                                                           #
#============================================================================

test_that(  "read_metabolon", {
    # Read
        file <- download_data('atkin.metabolon.xlsx')
        object <- read_metabolon(file)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_false( 'Group'    %in% svars(object))
        expect_true('subgroup' %in% svars(object))
})

test_that(  "read_metabolon: subgroupvar = 'Time' ", {
    # Read
        file <- download_data('atkin.metabolon.xlsx')
        object <- read_metabolon(file, subgroupvar = 'Time')
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('subgroup' %in% svars(object))
        expect_false('Time'     %in% svars(object))
})

test_that(  "read_metabolon: pca = TRUE ", {
    # Read
        file <- download_data('atkin.metabolon.xlsx')
        object <- read_metabolon(file, pca = TRUE)
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% svars(object)))
        expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% fvars(object)))
        expect_true('sample_id~pca' %in% names(metadata(object)))
})

test_that(  "read_metabolon: fit = 'limma' ", {
    # Read
        file <- download_data('atkin.metabolon.xlsx')
        object <- read_metabolon(file, fit = 'limma', block = 'Subject')
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'limma'))))
})

test_that(  "read_metabolon: fit = 'lm' ", {
    # Read
        file <- download_data('atkin.metabolon.xlsx')
        object <- read_metabolon(file, fit = 'lm', block = 'Subject')
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lm'))))
})

test_that("read_metabolon: fit = 'lme' ", {
    # Read
        file <- download_data('atkin.metabolon.xlsx')
        object <- read_metabolon(file, block = 'Subject', fit = 'lme')
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lme'))))
})

test_that("read_metabolon: fit = 'lmer' ",{
    # Read
        file <- download_data('atkin.metabolon.xlsx')
        object <- suppressWarnings(read_metabolon(file, block = 'Subject', fit = 'lmer'))
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lmer'))))
})

test_that("read_metabolon: fit = 'wilcoxon' ", {
    # Read
        file <- download_data('atkin.metabolon.xlsx')
        object <- read_metabolon(file, block = 'Subject', fit = 'wilcoxon')
    # Test
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'wilcoxon'))))
})
