context('read_somascan')       # late - early : 32 down, 75 up
metadata <- S4Vectors::metadata

test_that( "read_somascan(file) ", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('subgroup' %in% svars(object))
    expect_false('SampleGroup' %in% svars(object))
})

test_that(  "read_somascan(file, subgroupvar = 'Subject_ID')", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, subgroupvar = 'Subject_ID')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true('subgroup'    %in% svars(object))
    expect_false('Subject_ID' %in% svars(object))
    expect_true('SampleGroup' %in% svars(object))
})

test_that(  "read_somascan(file, pca = TRUE)", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, pca = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% svars(object)))
    expect_true(all(c('effect~sample_id~pca1', 'effect~sample_id~pca2') %in% fvars(object)))
    expect_true('sample_id~pca' %in% names(metadata(object)))
})

test_that(  "read_somascan(file, fit = 'limma')", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, fit = 'limma')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'limma'))))
})

test_that(  "read_somascan(file, fit = 'limma')", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, fit = 'limma')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'limma'))))
})


test_that(  "read_somascan(file, fit = 'lm')", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, fit = 'lm')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lm'))))
})

test_that("read_somascan(file, fit = 'limma', block = 'Subject_ID')", {
    file <- download_data('atkin18.somascan.adat')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'limma'))))
    # expect_true('dupcor' %in% names(metadata(object))) 
    # now only internal since data.table is returned
})

test_that(  "read_somascan(file, fit = 'limma', block = 'Subject_ID', plot = TRUE)", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, fit = 'limma', block = 'Subject_ID', plot = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'limma'))))
})

test_that("read_somascan(file, fit = 'lme', block = 'Subject_ID', plot = TRUE)", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, block = 'Subject_ID', fit = 'lme', plot = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lme'))))
})

test_that("read_somascan(file,  fit = 'lmer', block = 'Subject_ID')", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, block = 'Subject_ID', fit = 'lmer')
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'lmer'))))
})

test_that("read_somascan(file, fit = 'wilcoxon', block = 'Subject_ID', plot = TRUE)", {
    file <- download_data('atkin18.somascan.adat')
    object <- read_somascan(file, block = 'Subject_ID', fit = 'wilcoxon', plot = TRUE)
    expect_s4_class(object, 'SummarizedExperiment')
    expect_true(any(stri_detect_fixed(fvars(object), paste0(FITSEP, 'wilcoxon'))))
})
