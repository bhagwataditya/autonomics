require(testthat)
context('read_maxquant_proteingroups: billing16')
metadata <- S4Vectors::metadata

test_that("read_maxquant_proteingroups(file)", {
    file <- download_data('billing16.proteingroups.txt')
    object <- read_maxquant_proteingroups(file)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true( 'subgroup'    %in% svars(object))
})

test_that("read_maxquant_proteingroups(file, quantity)", {
    file <- download_data('billing16.proteingroups.txt')
    object <- read_maxquant_proteingroups(file, quantity = 'labeledintensity')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('log2labeledintensity' %in% assayNames(object))
})

test_that("read_maxquant_proteingroups(file, invert)", {
    file <- download_data('billing16.proteingroups.txt')
    invert <- c('EM_E', 'BM_E', 'BM_EM')
    object <- read_maxquant_proteingroups(file, invert = invert)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_setequal(levels(object$subgroup), c('E_BM','E_EM','EM_BM'))
    expect_setequal(object$sample_id, 
                        c("E_BM.R1", "E_BM.R2", "E_BM.R3", "E_EM.R1", "E_EM.R2", 
                          "E_EM.R3", "EM_BM.R1", "EM_BM.R2", "EM_BM.R3")) 
    object <- read_maxquant_proteingroups(file)
    object %<>% invert_subgroups(invert) 
    expect_setequal(levels(object$subgroup), c('E_BM','E_EM','EM_BM'))
    expect_setequal(object$sample_id, 
                        c("E_BM.R1", "E_BM.R2", "E_BM.R3", "E_EM.R1", "E_EM.R2", 
                          "E_EM.R3", "EM_BM.R1", "EM_BM.R2", "EM_BM.R3")) 
})

if (requireNamespace('Biostrings', quietly = TRUE)){
    test_that("read_maxquant_proteingroups(file, fastafile)", {
        file <- download_data('billing16.proteingroups.txt')
        fastafile <- download_data('uniprot_hsa_20140515.fasta')
        fastadt <- read_fastahdrs(fastafile)
        invert <- c('EM_E', 'BM_E', 'BM_EM')
        object0 <- read_maxquant_proteingroups(file, invert = invert)
        object  <- read_maxquant_proteingroups(file, invert = invert, fastadt = fastadt)
        expect_s4_class(object, 'SummarizedExperiment')
        expect_true('canonical' %in% fvars(object))
        expect_true(any(nchar(fdt(object)$uniprot) < nchar(fdt(object0)$uniprot)))
    })
}

test_that(  "read_maxquant_proteingroups(file, pca = TRUE)", {
    file <- download_data('billing16.proteingroups.txt')
    invert <- c('EM_E', 'BM_E', 'EM_BM')
    object <- read_maxquant_proteingroups(file, invert = invert, pca = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(         svars(object),  'pca1')))
    expect_true(any(stri_detect_fixed(         fvars(object),  'pca1')))
    expect_true(any(stri_detect_fixed(names(metadata(object)), 'pca')))
})

test_that(  "read_maxquant_proteingroups(file, impute = TRUE)", {
    file <- download_data('billing16.proteingroups.txt')
    invert <- c('EM_E', 'BM_E', 'EM_BM')
    object <- read_maxquant_proteingroups(file, invert = invert, impute = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('is_imputed' %in% SummarizedExperiment::assayNames(object))
})

test_that(  "read_maxquant_proteingroups(file, fit = 'limma')", {
    file <- download_data('billing16.proteingroups.txt')
    invert <- c('EM_E', 'BM_E', 'EM_BM')
    formula <- ~0+subgroup
    object <- read_maxquant_proteingroups(file, invert = invert, fit = 'limma', formula = formula, 
                    contrasts = c('BM_EM', 'E_EM', 'E_BM'))
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), 'limma')))
})

test_that(  "read_maxquant_proteingroups(file, fit = 'lm')", {
    file <- download_data('billing16.proteingroups.txt')
    invert <- c('EM_E', 'BM_E', 'EM_BM')
    formula <- ~0 + subgroup
    object <- read_maxquant_proteingroups(file, invert = invert, fit = 'lm', formula = formula)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), 'lm')))
})

test_that(  "read_maxquant_proteingroups(file, fit = 'wilcoxon')", {
    file <- download_data('billing16.proteingroups.txt')
    invert <- c('EM_E', 'BM_E', 'EM_BM')
    object <- read_maxquant_proteingroups(file, invert = invert, fit = 'wilcoxon', formula = formula)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), 'wilcoxon')))
})
