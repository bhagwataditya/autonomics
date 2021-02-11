require(magrittr)
context('add_limma')
    
    msg <- 'add_limma("billing19.proteingroups")'
    file <- download_data('billing19.proteingroups.txt')
    select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
    select %<>% paste0('_STD')
    object <- read_proteingroups(file, select_subgroups = select, plot = FALSE)
    object %<>% add_limma(plot=FALSE)
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    
    msg <- 'add_limma("halama18.metabolon")'
    file <- download_data('halama18.metabolon.xlsx')
    object <- read_metabolon(file, plot = FALSE)
    object %<>% add_limma(plot=FALSE)
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    
    msg <- 'add_limma("billing19.rnacounts")'
    file <- download_data('billing19.rnacounts.txt')
    object <- read_rnaseq_counts(file, plot=FALSE)
    object %<>% add_limma(plot=FALSE)
    test_that(msg, expect_s4_class(object, 'SummarizedExperiment'))
    

context('extract_limma_summary')

# RNASEQCOUNTS
    # ~ 0 + subgroup | weights
        msg <- 'extract_limma_summary("billing19.rnacounts.txt")'
        file <- download_data('billing19.rnacounts.txt')
        object <- read_rnaseq_counts(file, limma=TRUE, plot=FALSE)
        limmadt <- extract_limma_summary(object)
        test_that(msg, expect_s3_class(limmadt, 'data.table'))

    # ~ 0 + subgroups
        msg <- 'extract_limma_summary("billing19.rnacounts.txt" no voom)'
        weights(object) <- NULL
        object %<>% add_limma(plot=FALSE)
        limmadt <- extract_limma_summary(object)
        test_that(msg, expect_s3_class(limmadt, 'data.table'))

# METABOLON
    # ~ 0 + subgroup
        msg <- 'extract_limma_summary("atkin18.metabolon.xlsx")'
        file <- download_data('atkin18.metabolon.xlsx')
        object <- read_metabolon(file, limma=TRUE, block='SUB', plot=FALSE)
        extract_limma_summary(object)
        test_that(msg, expect_s3_class(limmadt, 'data.table'))

   # ~ 0 + subgroup | block
        msg <- 'extract_limma_summary("atkin18.metabolon.xlsx, block")'
        object %<>% add_limma(plot=FALSE, block = 'SUB')
        limmadt <- extract_limma_summary(object)
        test_that(msg, expect_s3_class(limmadt, 'data.table'))

   # ~ 0 + subgroup + t2d | block
        msg <- 'extract_limma_summary("atkin18.metabolon.xlsx, block, T2D")'
        object %<>% add_limma(formula=~0+subgroup+T2D, block='SUB', plot=FALSE)
        limmadt <- extract_limma_summary(object)
        test_that(msg, expect_s3_class(limmadt, 'data.table'))

context('plot_contrastogram')
    # subgroup vector
    msg <- 'plot_contrastogram("billing19.proteingroups")'
    file <-  download_data('billing19.proteingroups.txt')
    select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
    select %<>% paste0('_STD')
    object <- read_proteingroups(
                file, select_subgroups = select, limma=TRUE, plot=FALSE)
    test_that(msg, expect_error(plot_contrastogram(object, curve=0.8), NA))
    
    # subgroup vector
    msg <- 'plot_contrastogram("fukuda20.proteingroups.txt")'
    file <-  download_data('fukuda20.proteingroups.txt')
    object <- read_proteingroups(file, limma=TRUE, plot=FALSE)
    test_that(msg, expect_error(plot_contrastogram(object), NA))
    
    # Ratios: self-contrasts
    msg <- 'plot_contrastogram("billing16.proteingroups.txt")'
    file <- download_data('billing16.proteingroups.txt')
    invert <- c('EM_E', 'BM_E', 'BM_EM')
    object <- read_proteingroups(
               file, invert_subgroups=invert, limma=TRUE, plot=FALSE)
    test_that(msg, expect_error(plot_contrastogram(object), NA))

        
context('plot_volcano')
    # proteingroup group ratios
    msg  <- 'plot_volcano("billing16.proteingroups.txt")'
    file <- download_data("billing16.proteingroups.txt")
    inv <- c('EM_E', 'BM_E', 'BM_EM')
    object <- read_proteingroups(
        file, invert_subgroups=inv, limma=TRUE, plot=FALSE)
    test_that(msg, expect_s3_class(plot_volcano(object), 'gg'))
    
    # metabolon intensities: complex design
    msg  <- 'plot_volcano("halama18.metabolon.xlsx")'
    file <- download_data('halama18.metabolon.xlsx')
    object <- read_metabolon(file, limma=TRUE, plot=FALSE)
    test_that(msg, expect_s3_class(plot_volcano(object, ntop=0), 'gg'))

    # proteingroup internalstandard ratios
    msg  <- 'plot_volcano("billing19.proteingroups.txt")'
    file <-  download_data('billing19.proteingroups.txt')
    select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
    select %<>% paste0('_STD')
    object <- read_proteingroups(
                 file, select_subgroups = select, limma=TRUE, plot=FALSE)
    test_that(msg, expect_s3_class(plot_volcano(object), 'gg'))

