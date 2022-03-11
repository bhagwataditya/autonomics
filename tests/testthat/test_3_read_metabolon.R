context('read_metabolon')
metadata <- S4Vectors::metadata

test_that(  "read_metabolon(file) works", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'Group'    %in% svars(object))
    expect_false('subgroup' %in% svars(object))
})

test_that(  "read_metabolon(file, subgroupvar='SET')", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, subgroupvar = 'SET', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('SET' %in% svars(object))
    expect_false('Group'    %in% svars(object))
    expect_false('subgroup' %in% svars(object))
})

test_that(  "read_metabolon(file, subgroupvar='SET', pca=TRUE)", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, subgroupvar = 'SET', pca=TRUE, plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(all(c('pca1', 'pca2') %in% svars(object)))
    expect_true(all(c('pca1', 'pca2') %in% fvars(object)))
    expect_true('pca' %in% names(metadata(object)))
})

test_that(  "read_metabolon(file, subgroupvar='SET', fit='limma')", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, subgroupvar = 'SET', fit='limma', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(
        stri_detect_fixed(fvars(object), paste0(FITSEP, 'limma'))))
})

test_that(  "read_metabolon(file, subgroupvar='SET', fit='lm')", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(
                file, subgroupvar = 'SET', impute = TRUE, fit='lm', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(
        stri_detect_fixed(fvars(object), paste0(FITSEP, 'lm'))))
})

test_that("read_metabolon(file, subgroupvar='SET', fit='limma', block='SUB')", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- suppressWarnings(read_metabolon(file, subgroupvar = 'SET', 
                  block = 'SUB', impute = TRUE, fit = 'limma', plot = FALSE))
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(
        stri_detect_fixed(fvars(object), paste0(FITSEP, 'limma'))))
    expect_true('dupcor' %in% names(metadata(object)))
})

# test_that("read_metabolon(file, subgroupvar='SET', fit='lme', block='SUB')", {
#     file <- download_data('atkin18.metabolon.xlsx')
#     sgvar <- 'SET'; block <- 'SUB'
#     object <- read_metabolon(file, subgroupvar = sgvar, block = block, 
#                             impute = TRUE, fit = 'lme', plot = FALSE)
#     expect_s4_class(object, 'SummarizedExperiment')
#     expect_true('lme' %in% names(metadata(object)))
# })

# test_that("read_metabolon(file, subgroupvar='SET', fit='lmer', block='SUB')",{
#     file <- download_data('atkin18.metabolon.xlsx')
#     sgvar <- 'SET'; block <- 'SUB'
#     object <- read_metabolon(file, subgroupvar = sgvar, block = block, 
#                             impute = TRUE, fit = 'lmer', plot = FALSE)
#     expect_s4_class(object, 'SummarizedExperiment')
#     expect_true('lmer' %in% names(metadata(object)))
# })

test_that("read_metabolon(file,subgroupvar='SET',fit='wilcoxon',block='SUB')", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, subgroupvar = 'SET', block = 'SUB', 
                            fit='wilcoxon', plot=FALSE, impute=TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(
        stri_detect_fixed(fvars(object), paste0(FITSEP, 'wilcoxon'))))
})
