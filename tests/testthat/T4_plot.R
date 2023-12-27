require(testthat)
require(magrittr)

#============================================================================
#                                                                           #
             context(" plot_violins ")                                      #
#                                                                           #
#============================================================================


    file <- download_data('atkin.metabolon.xlsx')
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

