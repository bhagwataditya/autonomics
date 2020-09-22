#' @importFrom AnnotationDbi mapIds
#' @importFrom affy          just.rma
#' @importFrom assertive     assert_engine   assert_all_are_existing_files
#' @importFrom assertive     assert_all_are_in_closed_range
#' @importFrom assertive     assert_all_are_finite
#' @importFrom assertive     assert_all_are_greater_than
#' @importFrom assertive     assert_all_are_in_range
#' @importFrom assertive     assert_any_are_not_nan
#' @importFrom assertive     assert_all_are_whole_numbers
#' @importFrom assertive     assert_are_same_length
#' @importFrom assertive     has_names
#' @importFrom assertive     assert_is_all_of   assert_is_any_of
#' @importFrom assertive     assert_is_a_bool
#' @importFrom assertive     assert_is_a_number   assert_is_a_string
#' @importFrom assertive     assert_is_character  assert_is_data.frame
#' @importFrom assertive     assert_is_list
#' @importFrom assertive     assert_is_logical
#' @importFrom assertive     assert_is_non_empty
#' @importFrom assertive     assert_is_non_scalar   assert_is_not_null
#' @importFrom assertive     assert_is_numeric    assert_is_of_length
#' @importFrom assertive     assert_is_subset
#' @importFrom assertive     false   get_name_in_parent
#' @importFrom assertive     has_no_duplicates
#' @importFrom assertive     is_empty_character   is_existing_file
#' @importFrom colorspace    sequential_hcl
#' @importFrom data.table    data.table   fread   fwrite   setkeyv   set
#' @importFrom data.table    setnames    setorderv   tstrsplit   :=   .SD   .I
#' @importFrom edgeR         filterByExpr
#' @import     ggplot2
#' @importFrom graphics      lines    pie    title
#' @importFrom grDevices     hcl
#' @importFrom limma         duplicateCorrelation   lmFit   voom
#' @importFrom magrittr      %>%   %<>%   add   and   equals
#' @importFrom magrittr      divide_by   extract   extract2   is_in
#' @importFrom magrittr      set_colnames   set_names   set_rownames
#' @importFrom MASS          fitdistr
#' @importFrom matrixStats   rowAnys   colAnys   colWeightedMeans
#' @importFrom matrixStats   rowAlls   rowSds
#' @importFrom methods       as
#' @importFrom parallel      detectCores
#' @importFrom readxl        read_excel   excel_sheets
#' @importFrom rlang         enquo   eval_tidy   enexpr   expr_text   quo_name
#' @importFrom RCurl         url.exists
#' @importFrom Rsubread      featureCounts
#' @importFrom R.utils       gunzip
#' @importFrom S4Vectors     DataFrame   metadata   metadata<-
#' @importFrom seqinr        read.fasta
#' @importFrom stats         aggregate    approxfun   as.formula   lowess   sd
#' @importFrom stats         median   model.matrix   na.exclude   qnorm   rnorm
#' @import     stringi
#' @importFrom tidyr         separate_rows
#' @importFrom tools         file_ext   file_path_sans_ext
#' @importFrom SummarizedExperiment  assays   assays<-   SummarizedExperiment
#' @importFrom SummarizedExperiment  rowData   rowData<-   colData   colData<-
#' @importFrom SummarizedExperiment  makeSummarizedExperimentFromExpressionSet
#' @importFrom utils         adist   capture.output   count.fields
#' @importFrom utils         download.file    getFromNamespace
#' @importFrom utils         installed.packages    unzip
NULL

