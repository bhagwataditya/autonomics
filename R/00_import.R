#' @importFrom AnnotationDbi mapIds
#' @importFrom affy          just.rma
#' @importFrom assertive     assert_all_are_existing_files   assert_is_a_bool
#' @importFrom assertive     assert_is_a_number   assert_is_a_string
#' @importFrom assertive     assert_is_character  assert_is_data.frame
#' @importFrom assertive     assert_is_numeric    assert_is_subset
#' @importFrom assertive     is_existing_file
#' @importFrom assertive     has_no_duplicates
#' @importFrom data.table    fread   fwrite   data.table   setkeyv
#' @importFrom magrittr      %>%   %<>%   add   and   equals
#' @importFrom magrittr      extract   extract2   is_in
#' @importFrom magrittr      set_colnames   set_names
#' @importFrom methods       as
#' @importFrom parallel      detectCores
#' @importFrom purrr         detect_index
#' @importFrom readxl        read_excel   excel_sheets
#' @importFrom Rsubread      featureCounts
#' @importFrom R.utils       gunzip
#' @importFrom S4Vectors     DataFrame   metadata   metadata<-
#' @import     stringi
#' @importFrom tools         file_ext   file_path_sans_ext
#' @importFrom SummarizedExperiment  assays   assays<-   SummarizedExperiment
#' @importFrom SummarizedExperiment  rowData   rowData<-   colData   colData<-
#' @importFrom SummarizedExperiment  makeSummarizedExperimentFromExpressionSet
#' @importFrom utils         capture.output   count.fields   download.file
#' @importFrom utils         getFromNamespace   installed.packages   unzip
NULL

