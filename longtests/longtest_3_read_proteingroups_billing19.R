require(testthat)
context('read_maxquant_proteingroups: billing19')
metadata <- S4Vectors::metadata

test_that("read_maxquant_proteingroups(file)", {
    file <- download_data('billing19.proteingroups.txt')
    object <- read_maxquant_proteingroups(file)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'subgroup'    %in% svars(object))
})

test_that("read_maxquant_proteingroups(file, subgroups)", {
    file <- download_data('billing19.proteingroups.txt')
    subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_maxquant_proteingroups(file, subgroups = subgroups)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(assertive.sets::are_set_equal( slevels(object, 'subgroup'), subgroups))
})

test_that("read_maxquant_phosphosites(file)", {
    file        <- download_data('billing19.phosphosites.txt')
    proteinfile <- download_data('billing19.proteingroups.txt')
    subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_maxquant_phosphosites(file, proteinfile = proteinfile, subgroups = subgroups)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('log2proteins' %in% SummarizedExperiment::assayNames(object))
    expect_true('log2sites'    %in% SummarizedExperiment::assayNames(object))
    expect_true('log2diffs'    %in% SummarizedExperiment::assayNames(object))
})

test_that("read_maxquant_proteingroups(file, pca = TRUE)", {
    file <- download_data('billing19.proteingroups.txt')
    subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_maxquant_proteingroups(file, subgroups = subgroups, pca = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(         svars(object),  'pca1')))
    expect_true(any(stri_detect_fixed(         fvars(object),  'pca1')))
    expect_true(any(stri_detect_fixed(names(metadata(object)), 'pca' )))
})

test_that("read_maxquant_proteingroups(file, impute = TRUE)", {
    file <- download_data('billing19.proteingroups.txt')
    select <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_maxquant_proteingroups(file, subgroups = subgroups, impute = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('is_imputed' %in% SummarizedExperiment::assayNames(object))
})

test_that("read_maxquant_proteingroups(file, fit = 'limma')", {
    file <- download_data('billing19.proteingroups.txt')
    subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_maxquant_proteingroups(file, subgroups = subgroups, impute = TRUE, fit = 'limma')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), 'limma')))
})

test_that("read_maxquant_proteingroups(file, fit='lm')", {
    file <- download_data('billing19.proteingroups.txt')
    subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_maxquant_proteingroups(file, subgroups = subgroups, impute = TRUE, fit = 'lm')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), 'lm')))
})

test_that("read_maxquant_proteingroups(file, fit = 'wilcoxon')", {
    file <- download_data('billing19.proteingroups.txt')
    subgroups <- sprintf('%s_STD', c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_maxquant_proteingroups(file, subgroups = subgroups, impute = TRUE, fit = 'wilcoxon')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), 'wilcoxon')))
})
