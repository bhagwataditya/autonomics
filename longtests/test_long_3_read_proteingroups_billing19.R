context('read_maxquant_proteingroups: billing19')
metadata <- S4Vectors::metadata

test_that("read_maxquant_proteingroups(file)", {
    file <- download_data('billing19.proteingroups.txt')
    object <- read_maxquant_proteingroups(file, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'subgroup'    %in% svars(object))
})

test_that("read_maxquant_proteingroups(file, subgroupvar = NULL)", {
    file <- download_data('billing19.proteingroups.txt')
    object <- read_maxquant_proteingroups(file, subgroupvar = NULL, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'subgroup'    %in% svars(object))
})

test_that("read_maxquant_proteingroups(file, select_subgroups)", {
    file <- download_data('billing19.proteingroups.txt')
    select_subgroups <- sprintf('%s_STD', 
                                c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_maxquant_proteingroups(
                    file, select_subgroups = select_subgroups, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(assertive::are_set_equal(
                    slevels(object, 'subgroup'), select_subgroups))
})

test_that("read_phosphosites(file)", {
    file        <- download_data('billing19.phosphosites.txt')
    proteinfile <- download_data('billing19.proteingroups.txt')
    select_subgroups <- sprintf('%s_STD', 
                                c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_phosphosites(file, proteinfile =proteinfile, 
                            select_subgroups = select_subgroups, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('proteingroups' %in% SummarizedExperiment::assayNames(object))
    expect_true('phosphosites'  %in% SummarizedExperiment::assayNames(object))
    expect_true('occupancies'   %in% SummarizedExperiment::assayNames(object))
})

test_that("read_maxquant_proteingroups(file, pca = TRUE)", {
    file <- download_data('billing19.proteingroups.txt')
    select_subgroups <- sprintf('%s_STD', 
                                c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_maxquant_proteingroups(file, select_subgroups = select_subgroups, 
                                pca = TRUE, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(all(c('pca1', 'pca2') %in% svars(object)))
    expect_true(all(c('pca1', 'pca2') %in% fvars(object)))
    expect_true('pca' %in% names(metadata(object)))
})

test_that("read_maxquant_proteingroups(file, impute = TRUE)", {
    file <- download_data('billing19.proteingroups.txt')
    select_subgroups <- sprintf('%s_STD', 
                                c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_maxquant_proteingroups(file, select_subgroups = select_subgroups, 
                                impute = TRUE, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('is_imputed' %in% SummarizedExperiment::assayNames(object))
})

test_that("read_maxquant_proteingroups(file, fit = 'limma')", {
    file <- download_data('billing19.proteingroups.txt')
    select_subgroups <- sprintf('%s_STD', 
                                c('E00','E01','E02','E05','E15','E30','M00'))
    object <- read_maxquant_proteingroups(file, select_subgroups = select_subgroups, 
                                impute = TRUE, fit = 'limma', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('limma' %in% names(metadata(object)))
})

# test_that("read_maxquant_proteingroups(file, fit='lm')", {
#     file <- download_data('billing19.proteingroups.txt')
#     select_subgroups <- sprintf('%s_STD', 
#                                 c('E00','E01','E02','E05','E15','E30','M00'))
#     object <- read_maxquant_proteingroups(file, select_subgroups = select_subgroups, 
#                                 impute = TRUE, fit = 'lm', plot = FALSE)
#     expect_s4_class(object, 'SummarizedExperiment')
#     expect_true('lm' %in% names(metadata(object)))
# })

# test_that("read_maxquant_proteingroups(file, fit='wilcoxon')", {
#     file <- download_data('billing19.proteingroups.txt')
#     select_subgroups <- sprintf('%s_STD', 
#                                 c('E00','E01','E02','E05','E15','E30','M00'))
#     object <- read_maxquant_proteingroups(file, select_subgroups = select_subgroups, 
#                                 impute = TRUE, fit = 'wilcoxon', plot = FALSE)
#     expect_s4_class(object, 'SummarizedExperiment')
#     expect_true('wilcoxon' %in% names(metadata(object)))
#     expect_identical(names(dimnames(metadata(object)$wilcoxon))[2], 'subgroup')
# })
