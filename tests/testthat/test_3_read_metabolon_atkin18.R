context('read_metabolon')
metadata <- S4Vectors::metadata

test_that(  "read_metabolon(file) works", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'Group'    %in% svars(object))
    expect_false('subgroup' %in% svars(object))
})

test_that(  "read_metabolon(file, subgroupvar = NULL) works", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, subgroupvar = NULL, plot = FALSE)
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
    sgvar <- 'SET'
    object <- read_metabolon(file, subgroupvar = sgvar, fit='limma', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('limma' %in% names(metadata(object)))
})

test_that(  "read_metabolon(file, subgroupvar='SET', fit='lm')", {
    file <- download_data('atkin18.metabolon.xlsx')
    sgvar <- 'SET'
    object <- read_metabolon(
                file, subgroupvar = sgvar, impute = TRUE, fit='lm', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('lm' %in% names(metadata(object)))
})

test_that("read_metabolon(file, subgroupvar='SET', fit='limma', block='SUB')", {
    file <- download_data('atkin18.metabolon.xlsx')
    sgvar <- 'SET'; block <- 'SUB'
    object <- suppressWarnings(read_metabolon(file, subgroupvar = sgvar, 
                  block = block, impute = TRUE, fit = 'limma', plot = FALSE))
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('limma' %in% names(metadata(object)))
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
    sgvar <- 'SET'; block <- 'SUB'
    object <- read_metabolon(file, subgroupvar = sgvar, block = block, 
                            fit = 'wilcoxon', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('wilcoxon' %in% names(metadata(object)))
})
