utils::globalVariables('.')


#' Conveniently message dataframe
#'
#' Conveniently message dataframe using sprintf syntax.
#' Use place holder '%s' for data.frame.
#'
#' @param format_string sprintf style format string
#' @param x data.frame
#' @return NULL
#' @examples
#' x <- data.frame(feature_id = c('F001', 'F002'), symbol = c('FEAT1', 'FEAT2'))
#' cmessage_df('\t%s', x)
#'
#' x <- c(rep('PASS', 25), rep('FAIL', 25))
#' cmessage_df(format_string = '%s', table(x))
#' @noRd
cmessage_df <- function(format_string, x){
    format_string %>%
    sprintf(capture.output(print(x))) %>%
    paste0(collapse = '\n') %>%
    enc2utf8() %>%
    message()
}


#' Make vector components unique by appending spaces
#' @param x character or factor vector
#' @param method 'make.unique' or 'make.unique.spaces'
#' @param verbose TRUE (default) or FALSE
#' @return character vector
#' @seealso \code{\link[base]{make.unique}}
#' @examples
#' x <- c('A', 'B', 'C', 'A', 'D')
#' uniquify(x, 'make.unique')
#' uniquify(x, 'make.unique.spaces')
#' @noRd
uniquify <- function(x, method = 'make.unique.spaces', verbose = TRUE){
    idx <- cduplicated(x)
    if (any(idx)){
        uniquefun <- get(method)
        uniquex <- uniquefun(x)
        if (verbose){
            message(   '\t\tUniquify ( %s -> %s ) duplicates of',
                        x[idx][1], uniquefun(x[idx][c(1, 1)])[2])
            cmessage_df("\t\t\t  %s", table(x[idx]) %>% as.list() %>% unlist)
                # unlist(as.list(.)) is to prevent empty line
        }
    } else {
        uniquex <- x
    }
    uniquex
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
#' @noRd
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

  df %>% extract(, c(first_cols, setdiff(names(df), first_cols)), drop = FALSE)
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
