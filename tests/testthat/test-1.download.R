for (dataset in AUTONOMICS_DATASETS){
    test_that( sprintf('download_data(%s)', dataset),
               expect_type(download_data(dataset), 'character'))
}

