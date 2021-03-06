require(magrittr)
context('plot_violins')
    
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, fit='limma', plot = FALSE)
    object %<>% extract(, order(.$Group))
    control_features <- c('biotin','phosphate')
    fdata(object) %<>% cbind(control=.$feature_name %in% control_features)
    
    msg <- 'plot_violins("halama18.metabolon.xlsx")'
    test_that(msg, expect_s3_class(plot_violins(
        object[1:12, ], x=feature_id, fill=feature_id), 'gg'))

    test_that(msg, expect_s3_class(plot_violins(
        object[, 1:12], x=sample_id,  fill=sample_id ), 'gg'))

    test_that(msg, expect_s3_class(plot_violins(
        object[, 1:12], x=sample_id,  fill=Group), 'gg'))
    
    test_that(msg, expect_s3_class(plot_violins(
        object[, 1:12], x=Group,   fill=Group, group=sample_id), 'gg'))
    
    test_that(msg, expect_s3_class(plot_violins(
        object[1:4, ],  x=Group,   fill=Group, facet=feature_id), 'gg'))

