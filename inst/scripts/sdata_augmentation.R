write_sdata_template <- function(object, file = tempfile()) {
    data.table::fwrite(sdata(object), file, sep = "\t", row.names = TRUE)
    invisible(file)
}

add_additonal_sdata <- function(object, file = NULL, conservative = TRUE) {
    original_sdata <- sdata((object))
    new_sdata      <-
        fread(file, sep = "\t", integer64 = "numeric", data.table = FALSE) %>%
        set_rownames(.[["V1"]]) %>%
        dplyr::select(-V1)
    assert_are_identical(nrow(original_sdata), nrow(new_sdata))
    if (conservative) {
        subset_new_sdata <- transfer_df_str(
            original_sdata, new_sdata[names(original_sdata)])
        assert_are_identical(original_sdata, subset_new_sdata)
        new_sdata[names(original_sdata)] <- subset_new_sdata
    }
    sdata(object) <- new_sdata
}

transfer_df_str <- function(origin, target) {
    assert_is_data.frame(origin)
    assert_is_data.frame(target)
    assert_are_identical(dim(origin), dim(target))

    mapply( FUN = as, target, sapply(origin, class), SIMPLIFY = FALSE) %>%
    as.data.frame() %>%
    set_rownames(rownames(target))
}
