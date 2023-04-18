context('read_metabolon')
metadata <- S4Vectors::metadata

test_that(  "read_metabolon(file) works", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_false( 'Group'    %in% svars(object))
    expect_true('subgroup' %in% svars(object))
})

test_that(  "read_metabolon(file, subgroupvar = 'SET')", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, subgroupvar = 'SET')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('subgroup' %in% svars(object))
    expect_false('SET'     %in% svars(object))
})

test_that(  "read_metabolon(file, pca = TRUE)", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, pca = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% svars(object)))
    expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% fvars(object)))
    expect_true('sample_id~pca' %in% names(metadata(object)))
})

test_that(  "read_metabolon(file, fit = 'limma')", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, fit = 'limma', block = 'SUB', plot = TRUE, label = NULL)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'limma'))))
})

test_that(  "read_metabolon(file, fit = 'lm')", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, fit = 'lm', block = 'SUB', plot = TRUE, label = NULL)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lm'))))
})

test_that("read_metabolon(file, fit = 'lme')", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, block = 'SUB', fit = 'lme', plot = TRUE, label = NULL)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lme'))))
})

test_that("read_metabolon(file, fit = 'lmer')",{
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, block = 'SUB', fit = 'lmer', plot = TRUE, label = NULL)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lmer'))))
})

test_that("read_metabolon(file, fit = 'wilcoxon')", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, block = 'SUB', fit = 'wilcoxon', plot = TRUE, label = NULL)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'wilcoxon'))))
})
