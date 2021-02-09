context('download_omics_data')
for (dataset in AUTONOMICS_DATASETS){
    test_that(
        paste0(dataset, ' works'),
        expect_type(download_data(dataset), 'character')
    )
}

