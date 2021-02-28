context('read_metabolon')       # late - early : 32 down, 75 up
metadta <- S4Vectors::metadata

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

test_that(  "read_metabolon(file, subgroupvar = 'SET') works", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, subgroupvar = 'SET', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('SET' %in% svars(object))
    expect_false('Group'    %in% svars(object))
    expect_false('subgroup' %in% svars(object))
})

test_that(  "read_metabolon(file, subgroupvar = 'SET', pca=TRUE) works", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, subgroupvar = 'SET', pca=TRUE, plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(all(c('pca1', 'pca2') %in% svars(object)))
    expect_true(all(c('pca1', 'pca2') %in% fvars(object)))
    expect_true('pca' %in% names(metadata(object)))
})

test_that(  "read_metabolon(file, subgroupvar = 'SET', fit='limma') works", {
    file <- download_data('atkin18.metabolon.xlsx')
    sgvar <- 'SET'
    object <- read_metabolon(file, subgroupvar = sgvar, fit='limma', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('limma' %in% names(metadata(object)))
    formulastr <- formula2str(default_formula(
                                object, subgroupvar = sgvar, fit = 'limma'))
    expect_identical(names(dimnames(limma(object)))[2], formulastr)
    names(dimnames(limma(object)))
})
