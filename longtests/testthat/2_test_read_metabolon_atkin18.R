context('read_metabolon')       # late - early : 32 down, 75 up
metadata <- S4Vectors::metadata

test_that(  "read_metabolon(file) works", {
    file <- download_data('halama18.metabolon.xlsx')
    object <- read_metabolon(file, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'Group'    %in% svars(object))
    expect_false('subgroup' %in% svars(object))
})

test_that(  "read_metabolon(file, subgroupvar = NULL) works", {
    file <- download_data('halama18.metabolon.xlsx')
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
    formulastr <- formula2str(default_formula(
                                object, subgroupvar = sgvar, fit = 'limma'))
    expect_identical(names(dimnames(limma(object)))[2], formulastr)
})

test_that(  "read_metabolon(file, subgroupvar='SET', fit='lm')", {
    file <- download_data('atkin18.metabolon.xlsx')
    sgvar <- 'SET'
    object <- read_metabolon(
                file, subgroupvar = sgvar, impute = TRUE, fit='lm', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('lm' %in% names(metadata(object)))
    formulastr <- formula2str(default_formula(
                                object, subgroupvar = sgvar, fit = 'lm'))
    expect_identical(names(dimnames(metadata(object)$lm))[2], formulastr)
})

test_that("read_metabolon(file, subgroupvar='SET', fit='limma', block='SUB')", {
    file <- download_data('atkin18.metabolon.xlsx')
    sgvar <- 'SET'; block <- 'SUB'
    object <- read_metabolon(file, subgroupvar = sgvar, block = block, 
                            impute = TRUE, fit = 'limma', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('limma' %in% names(metadata(object)))
    formulastr <- formula2str(default_formula(
                                object, subgroupvar = sgvar, fit = 'limma'))
    expect_identical(names(dimnames(metadata(object)$limma))[2], formulastr)
    expect_true('dupcor' %in% names(metadata(object)))
})

test_that("read_metabolon(file, subgroupvar='SET', fit='lme', block='SUB')", {
    file <- download_data('atkin18.metabolon.xlsx')
    sgvar <- 'SET'; block <- 'SUB'
    object <- read_metabolon(file, subgroupvar = sgvar, block = block, 
                            impute = TRUE, fit = 'lme', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('lme' %in% names(metadata(object)))
    formulastr <- formula2str(default_formula(
                                object, subgroupvar = sgvar, fit = 'lme'))
    expect_identical(names(dimnames(metadata(object)$lme))[2], formulastr)
})

test_that("read_metabolon(file, subgroupvar='SET', fit='lmer', block='SUB')", {
    file <- download_data('atkin18.metabolon.xlsx')
    sgvar <- 'SET'; block <- 'SUB'
    object <- read_metabolon(file, subgroupvar = sgvar, block = block, 
                            impute = TRUE, fit = 'lmer', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('lmer' %in% names(metadata(object)))
    formulastr <- formula2str(default_formula(
                                object, subgroupvar = sgvar, fit = 'lmer'))
    expect_identical(names(dimnames(metadata(object)$lme))[2], formulastr)
})

test_that("read_metabolon(file,subgroupvar='SET',fit='wilcoxon',block='SUB')", {
    file <- download_data('atkin18.metabolon.xlsx')
    sgvar <- 'SET'; block <- 'SUB'
    object <- read_metabolon(file, subgroupvar = sgvar, block = block, 
                            fit = 'wilcoxon', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('wilcoxon' %in% names(metadata(object)))
    expect_identical(names(dimnames(metadata(object)$wilcoxon))[2], sgvar)
})
