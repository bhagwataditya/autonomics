#' @importFrom AnnotationDbi mapIds
#' @importFrom affy          just.rma
#' @importFrom assertive     assert_all_are_existing_files   assert_is_a_bool
#' @importFrom assertive     assert_is_character  assert_is_a_number
#' @importFrom assertive     assert_is_a_string   assert_is_subset
#' @importFrom assertive     is_existing_file
#' @importFrom assertive     has_no_duplicates
#' @importFrom data.table    fread
#' @importFrom magrittr      %>%   %<>%   add   and   equals
#' @importFrom magrittr      extract   extract2   is_in
#' @importFrom magrittr      set_colnames   set_names
#' @importFrom purrr         detect_index
#' @importFrom readxl        read_excel
#' @importFrom Rsubread      featureCounts
#' @importFrom S4Vectors     DataFrame   metadata   metadata<-
#' @importFrom stringi       stri_count_fixed
#' @importFrom stringi       stri_detect_fixed   stri_detect_regex
#' @importFrom stringi       stri_extract_all_regex
#' @importFrom stringi       stri_extract_first_regex   stri_extract_first_words
#' @importFrom stringi       stri_extract_last_regex
#' @importFrom stringi       stri_split_fixed   stri_split_regex
#' @importFrom stringi       stri_replace_first_fixed   stri_replace_first_regex
#' @importFrom SummarizedExperiment  assays   assays<-   SummarizedExperiment
#' @importFrom SummarizedExperiment  rowData   rowData<-   colData   colData<-
#' @importFrom SummarizedExperiment  makeSummarizedExperimentFromExpressionSet
#' @importFrom utils         count.fields   getFromNamespace
NULL
