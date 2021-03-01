context('read_somscan')       # late - early : 32 down, 75 up
metadata <- S4Vectors::metadata

test_that( "read_somascan(file) ", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'SampleGroup'    %in% svars(object))
    expect_false('subgroup' %in% svars(object))
})

test_that(  "read_somascan(file, subgroupvar = NULL) works", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, subgroupvar = NULL, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'SampleGroup'    %in% svars(object))
    expect_false('subgroup' %in% svars(object))
})

test_that(  "read_somascan(file, subgroupvar='SampleGroup')", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, subgroupvar = 'SampleGroup', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('SampleGroup'    %in% svars(object))
    expect_false('subgroup' %in% svars(object))
})

test_that(  "read_somascan(file, pca=TRUE)", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, pca=TRUE, plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(all(c('pca1', 'pca2') %in% svars(object)))
    expect_true(all(c('pca1', 'pca2') %in% fvars(object)))
    expect_true('pca' %in% names(metadata(object)))
})

test_that(  "read_somascan(file, fit='limma')", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, fit='limma', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('limma' %in% names(metadata(object)))
    formulastr <- formula2str(
        default_formula(object, subgroupvar = 'SampleGroup', fit = 'limma'))
    expect_identical(names(dimnames(metadata(object)$limma))[2], formulastr)
})

test_that(  "read_somascan(file, fit='lm')", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, fit='lm', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('lm' %in% names(metadata(object)))
    formulastr <- formula2str(default_formula(
                                object, subgroupvar = 'SampleGroup', fit = 'lm'))
    expect_identical(names(dimnames(metadata(object)$lm))[2], formulastr)
})

test_that("read_somascan(file, fit='limma', block='Subject_ID')", {
    file <- download_data('atkin18.somascan.adat')
    block <- 'Subject_ID'
    object <- read_somascan(file, block = block, fit = 'limma', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('limma' %in% names(metadata(object)))
    formulastr <- formula2str(default_formula(
                            object, subgroupvar = 'SampleGroup', fit = 'limma'))
    expect_identical(names(dimnames(metadata(object)$limma))[2], formulastr)
    expect_true('dupcor' %in% names(metadata(object)))
})

test_that("read_somascan(file, fit='lme', block='Subject_ID')", {
    file <- download_data('atkin18.somascan.adat')
    block <- 'Subject_ID'
    object <- read_somascan(file, block = block, fit = 'lme', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('lme' %in% names(metadata(object)))
    formulastr <- formula2str(default_formula(
                            object, subgroupvar = 'SampleGroup', fit = 'lme'))
    expect_identical(names(dimnames(metadata(object)$lme))[2], formulastr)
})

test_that("read_somascan(file,  fit='lmer', block='Subject_ID')", {
    file <- download_data('atkin18.somascan.adat')
    block <- 'Subject_ID'
    object <- read_somascan(file, block = block, fit = 'lmer', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('lmer' %in% names(metadata(object)))
    formulastr <- formula2str(default_formula(
                            object, subgroupvar = 'SampleGroup', fit = 'lmer'))
    expect_identical(names(dimnames(metadata(object)$lme))[2], formulastr)
})

test_that("read_somascan(file, fit='wilcoxon', block='Subject_ID')", {
    file <- download_data('atkin18.somascan.adat')
    block <- 'Subject_ID'
    object <- read_somascan(file, block = block, fit = 'wilcoxon', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('wilcoxon' %in% names(metadata(object)))
    expect_identical(names(dimnames(metadata(object)$wilcoxon))[2],
                    'SampleGroup')
})
