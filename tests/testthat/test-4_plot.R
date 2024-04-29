require(magrittr)

#============================================================================
#                                                                           #
#            context(" plot_violins ")                                      #
#                                                                           #
#============================================================================


    file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
    object <- read_metabolon(file, fit = 'limma')
    object %<>% extract(, order(.$subgroup))
    control_features <- c('biotin','phosphate')
    fdata(object) %<>% cbind(control = .$feature_name %in% control_features)
    
    msg <- 'plot_violins: atkin.metabolon.xlsx'
    test_that(msg, expect_s3_class(
         plot_violins(object[1:12, ], x = 'feature_id', fill = 'feature_id'),
        'gg'))
    
    test_that(msg, expect_s3_class(
         plot_violins(object[, 1:12], x = 'sample_id',  fill = 'sample_id' ),
        'gg'))
    
    test_that(msg, expect_s3_class(
         plot_violins(object[, 1:12], x = 'sample_id',  fill = 'subgroup'),
        'gg'))
    
    test_that(msg, expect_s3_class(
         plot_violins(object[, 1:12], x = 'subgroup',   fill = 'subgroup',   group = 'sample_id'),
        'gg'))
    
    test_that(msg, expect_s3_class(
         plot_violins(object[1:4,  ], x = 'subgroup',   fill = 'subgroup',   facet = 'feature_id'),
        'gg'))


#============================================================================
#                                                                           #
#            context(" plot_volcano ")                                      #
#                                                                           #
#============================================================================
    
# test_that('plot_volcano: billing16.proteingroups ', {
#     file <- download_data("billing16.proteingroups.txt")
#     invert <- c('EM_E', 'BM_E', 'BM_EM')
#     object <- read_maxquant_proteingroups(file, invert = invert, fit = 'limma')
#     expect_s3_class(plot_volcano(object), 'gg')
# })

test_that(' plot_volcano: billing19.proteingroups ', {
    file <-  system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
    subgroups <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
    subgroups %<>% paste0('_STD')
    object <- read_maxquant_proteingroups(file, subgroups = subgroups, fit = 'limma')
    expect_s3_class(plot_volcano(object), 'gg')
})

test_that('plot_volcano: atkin.metabolon ', {
    file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
    object <- read_metabolon(file, fit = 'limma')
    expect_s3_class(plot_volcano(object), 'gg')
})

    