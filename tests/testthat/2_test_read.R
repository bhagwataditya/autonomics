context('read_metabolon')       # late - early : 32 down, 75 up

test_that(  "read_metabolon(file)", {
    file <- download_data('halama18.metabolon.xlsx')
    object <- read_metabolon(file, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('Group' %in% svars(object))
    expect_false('subgroup' %in% svars(object))
})

test_that(  "read_metabolon(file, subgroupvar = NULL)", {
    file <- download_data('halama18.metabolon.xlsx')
    object <- read_metabolon(file, subgroupvar = NULL, plot = FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('Group' %in% svars(object))
    expect_false('subgroup' %in% svars(object))
})

test_that(  "read_metabolon(file, subgroupvar = 'SET')", {
    file <- download_data('atkin18.metabolon.xlsx')
    object <- read_metabolon(file, subgroupvar = 'SET', plot=FALSE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('SET' %in% svars(object))
    expect_false('Group' %in% svars(object))
    expect_false('subgroup' %in% svars(object))
    object$Group
})


