
#' message dataframe
#'
#' message dataframe using sprintf syntax.
#' Use place holder '%s' for data.frame.
#'
#' @param format_string sprintf style format string
#' @param x data.frame
#' @return nothing returned
#' @examples
#' x <- data.frame(feature_id = c('F001', 'F002'), symbol = c('FEAT1', 'FEAT2'))
#' message_df('\t%s', x)
#'
#' x <- c(rep('PASS', 25), rep('FAIL', 25))
#' message_df(format_string = '%s', table(x))
#' @export
message_df <- function(format_string, x){
    format_string %>%
    sprintf(capture.output(print(x))) %>%
    paste0(collapse = '\n') %>%
    enc2utf8() %>%
    message()
}


#' Convenient (two way) duplicated
#' @param x vector
#' @return logical vector
#' @examples
#' require(magrittr)
#' c(1,2,3,4,5,2) %>% cduplicated()
#' @noRd
cduplicated <- function(x){
    duplicated(x) | duplicated(x, fromLast = TRUE)
}


#' Pull columns in a dataframe to the front
#' @param df         data.frame
#' @param first_cols character vector: columns to be pulled to the front
#' @param verbose    TRUE (default) or FALSE
#' @return dataframe with re-ordered columns
#' @examples
#' df <- data.frame(
#'    symbol = c('A1BG', 'A2M'),
#'    id     = c('1',    '2'),
#'    name   = c('alpha-1-B glycoprotein', 'alpha-2-macroglobulin'),
#'    type   = c('proteincoding', 'proteincoding'))
#' first_cols <- c('id', 'symbol', 'location', 'uniprot')
#' pull_columns(df, first_cols)
#' @export
pull_columns <- function(df, first_cols, verbose = TRUE){

    assert_is_data.frame(df)
    assert_is_character(first_cols)

    idx <- first_cols %in% names(df)
    if (any(!idx)){
        if (verbose) message(
            'pull_columns: ignore absent columns ',
            paste0(sprintf("'%s'", first_cols[!idx]), collapse = ', '))
        first_cols %<>% extract(idx)
    }

    if (is.data.table(df)){
        df[, c(first_cols, setdiff(names(df), first_cols)), with = FALSE]
    } else {
        df[, c(first_cols, setdiff(names(df), first_cols)), drop = FALSE]
    }
}



# note: earlier name was 'order_pres_factor'
#' Create factor with levels in order of appearance
#'
#' Creates a factor from a vector, where the levels are in (possibly reverse)
#' order of appearance in the vector, rather than being alphabetically sorted.
#' @param avector An atomic vector.
#' @param reverse FALSE (default) or TRUE: reverse order factor levels?
#' @return factor vector
#' @examples
#' factorify(c('x', 'z', 'a'))
#' @noRd
factorify <- function(avector, reverse = FALSE){
    myLevels <- unique(avector)
    if(reverse){myLevels <- rev(myLevels)}
    factor(avector, myLevels)
}
