context('read_proteingroups: billing16')
metadata <- S4Vectors::metadata
require(magrittr)

test_that("read_proteingroups(file)", {
    file <- download_data('billing16.proteingroups.txt')
    object <- read_proteingroups(file, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'subgroup'    %in% svars(object))
})

test_that("read_proteingroups(file, subgroupvar = NULL)", {
    file <- download_data('billing16.proteingroups.txt')
    object <- read_proteingroups(file, subgroupvar = NULL, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'subgroup'    %in% svars(object))
})

test_that("read_proteingroups(file, quantity)", {
    file <- download_data('billing16.proteingroups.txt')
    object <- read_proteingroups(
                file, quantity = 'Intensity labeled', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_identical(metadata(object)$quantity, 'Intensity labeled')
})

test_that("read_proteingroups(file, invert_subgroups)", {
    file <- download_data('billing16.proteingroups.txt')
    invert_subgroups <- c('EM_E', 'BM_E', 'BM_EM')
    object <- read_proteingroups(
                    file, invert_subgroups = invert_subgroups, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_setequal(levels(object$subgroup), c('E_BM','E_EM','EM_BM'))
    expect_setequal(object$sample_id, 
                        c("E_BM.R1", "E_BM.R2", "E_BM.R3", "E_EM.R1", "E_EM.R2", 
                          "E_EM.R3", "EM_BM.R1", "EM_BM.R2", "EM_BM.R3")) 
    object <- read_proteingroups(file, plot = FALSE)
    object %<>% invert(subgroups = invert_subgroups) 
    expect_setequal(levels(object$subgroup), c('E_BM','E_EM','EM_BM'))
    expect_setequal(object$sample_id, 
                        c("E_BM.R1", "E_BM.R2", "E_BM.R3", "E_EM.R1", "E_EM.R2", 
                          "E_EM.R3", "EM_BM.R1", "EM_BM.R2", "EM_BM.R3")) 
})

if (requireNamespace('Biostrings', quietly = TRUE)){
    test_that("read_proteingroups(file, fastafile)", {
        file <- download_data('billing16.proteingroups.txt')
        fastafile <- download_data('uniprot_hsa_20140515.fasta')
        fastadt <- read_fastahdrs(fastafile)
        invert <- c('EM_E', 'BM_E', 'BM_EM')
        object0 <- read_proteingroups(file, invert = invert, plot=FALSE)
        object  <- read_proteingroups(file, invert = invert, fastadt = fastadt, plot = FALSE)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('canonical' %in% fvars(object))
        expect_true(any(nchar(fdata(object)$uniprot) <
                        nchar(fdata(object0)$uniprot)))
    })
}

test_that(  "read_proteingroups(file, pca=TRUE)", {
    file <- download_data('billing16.proteingroups.txt')
    invert_subgroups <- c('EM_E', 'BM_E', 'EM_BM')
    object <- read_proteingroups(file, invert_subgroups = invert_subgroups, 
                                pca = TRUE, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(all(c('pca1', 'pca2') %in% svars(object)))
    expect_true(all(c('pca1', 'pca2') %in% fvars(object)))
    expect_true('pca' %in% names(metadata(object)))
})

test_that(  "read_proteingroups(file, impute=TRUE)", {
    file <- download_data('billing16.proteingroups.txt')
    invert_subgroups <- c('EM_E', 'BM_E', 'EM_BM')
    object <- read_proteingroups(file, invert_subgroups = invert_subgroups, 
                                impute = TRUE, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('is_imputed' %in% SummarizedExperiment::assayNames(object))
})

test_that(  "read_proteingroups(file, fit='limma')", {
    file <- download_data('billing16.proteingroups.txt')
    invert_subgroups <- c('EM_E', 'BM_E', 'EM_BM')
    formula <- ~0+subgroup
    object <- read_proteingroups(file, invert_subgroups = invert_subgroups, 
                    fit = 'limma', formula = formula, 
                    contrasts = c('BM_EM', 'E_EM', 'E_BM'), plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('limma' %in% names(metadata(object)))
})

# test_that(  "read_proteingroups(file, fit='lm')", {
#     file <- download_data('billing16.proteingroups.txt')
#     invert_subgroups <- c('EM_E', 'BM_E', 'EM_BM')
#     formula <- ~0 + subgroup
#     object <- read_proteingroups(file, invert_subgroups = invert_subgroups, 
#                             fit = 'lm', formula = formula, plot = FALSE)
#     expect_s4_class(object, 'SummarizedExperiment')
#     expect_true('lm' %in% names(metadata(object)))
# })

test_that(  "read_proteingroups(file, fit='wilcoxon')", {
    file <- download_data('billing16.proteingroups.txt')
    invert_subgroups <- c('EM_E', 'BM_E', 'EM_BM')
    object <- read_proteingroups(file, invert_subgroups = invert_subgroups, 
                    fit = 'wilcoxon', subgroupvar = 'subgroup', plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('wilcoxon' %in% names(metadata(object)))
})
