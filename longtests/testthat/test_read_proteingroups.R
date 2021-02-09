context('read_proteingroups / read_phosphosites')

# invert_ratios
    file <- download_data('billing16.proteingroups.txt')
    inv <- c('EM_E','BM_E','BM_EM')
    object <- read_proteingroups(
                file, invert_subgroups = inv, 
                contrastdefs = c('E_EM', 'E_BM', 'EM_BM'), 
                pca = FALSE, lmfit = FALSE, plot = FALSE)
    msg <- 'read_proteingroups(invert_ratios) works'
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
# quantity
    file <- download_data('billing16.proteingroups.txt')
    object <- read_proteingroups(
                file, quantity = 'Intensity labeled', 
                pca = FALSE, lmfit = FALSE, plot = FALSE)
    msg <- 'read_proteingroups(select) works'
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
# select_subgroups
    require(magrittr)
    file <-  download_data('billing19.proteingroups.txt')
    select_subgroups <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
    select_subgroups %<>% paste0('_STD')
    object <- read_proteingroups(
                file, select_subgroups = select_subgroups, 
                pca = FALSE, lmfit = FALSE, plot = FALSE)
    msg <- 'read_proteingroups(select) works'
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
# fastafile
    file <-  download_data('billing16.proteingroups.txt')
    fastafile <- download_data('uniprot_hsa_20140515.fasta')
    object <- read_proteingroups(
                file, fastafile=fastafile, 
                pca = FALSE, lmfit = FALSE, plot = FALSE)
    msg <- 'read_proteingroups(fastafile) works'
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
# phospho
    phosphofile <- download_data('billing19.phosphosites.txt')
    proteinfile <- download_data('billing19.proteingroups.txt')
    select_subgroups <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
    select_subgroups %<>% paste0('_STD')
    object <- read_phosphosites(
                phosphofile, proteinfile, select_subgroups = select_subgroups, 
                pca = FALSE, lmfit = FALSE, plot = FALSE)
    msg <- 'read_phosphosites() works'
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))

