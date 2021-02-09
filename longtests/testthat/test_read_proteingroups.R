context('invert')
    require(magrittr)
    file <- download_data('billing16.proteingroups.txt')
    object <- read_proteingroups(file, pca=FALSE, lmfit=FALSE, plot=FALSE)
    object %<>% invert(subgroups = c('EM_E','BM_E','BM_EM'))
    msg <- 'invert_ratios'
    test_that(msg, {
        expect_s4_class(object, 'SummarizedExperiment')
        expect_setequal(levels(object$subgroup), c('E_BM','E_EM','EM_BM'))
        expect_setequal(object$sample_id, 
                        c("E_BM.R1", "E_BM.R2", "E_BM.R3", "E_EM.R1", "E_EM.R2", 
                          "E_EM.R3", "EM_BM.R1", "EM_BM.R2", "EM_BM.R3")) })

context('read_proteingroups / read_phosphosites')
    # invert_ratios
        file <- download_data('billing16.proteingroups.txt')
        inv <- c('EM_E','BM_E','BM_EM')
        object <- read_proteingroups(
                    file, invert_subgroups = inv, 
                    contrastdefs = c('E_EM', 'E_BM', 'EM_BM'), 
                    pca = FALSE, lmfit = FALSE, plot = FALSE)
        msg <- 'read_proteingroups(invert_ratios)'
        test_that(msg, {
            expect_s4_class(object, 'SummarizedExperiment')
            expect_setequal(levels(object$subgroup), c('E_BM','E_EM','EM_BM'))
            expect_setequal(object$sample_id, 
                        c("E_BM.R1", "E_BM.R2", "E_BM.R3", "E_EM.R1", "E_EM.R2", 
                          "E_EM.R3", "EM_BM.R1", "EM_BM.R2", "EM_BM.R3")) })
    # quantity
        file <- download_data('billing16.proteingroups.txt')
        object <- read_proteingroups(
                    file, quantity = 'Intensity labeled', 
                    pca = FALSE, lmfit = FALSE, plot = FALSE)
        msg <- 'read_proteingroups(select)'
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    # select_subgroups
        require(magrittr)
        file <-  download_data('billing19.proteingroups.txt')
        select_subgroups <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
        select_subgroups %<>% paste0('_STD')
        object <- read_proteingroups(
                    file, select_subgroups = select_subgroups, 
                    pca = FALSE, lmfit = FALSE, plot = FALSE)
        msg <- 'read_proteingroups(select)'
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    # fastafile
        file <-  download_data('billing16.proteingroups.txt')
        fastafile <- download_data('uniprot_hsa_20140515.fasta')
        object <- read_proteingroups(
                    file, fastafile=fastafile, 
                    pca = FALSE, lmfit = FALSE, plot = FALSE)
        msg <- 'read_proteingroups(fastafile)'
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    # phospho
        phosphofile <- download_data('billing19.phosphosites.txt')
        proteinfile <- download_data('billing19.proteingroups.txt')
        select_subgroups <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
        select_subgroups %<>% paste0('_STD')
        object <- read_phosphosites(
                    phosphofile, proteinfile, 
                    select_subgroups = select_subgroups, 
                    pca = FALSE, lmfit = FALSE, plot = FALSE)
        msg <- 'read_phosphosites()'
        test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))

