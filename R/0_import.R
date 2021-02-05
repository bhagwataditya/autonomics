#' @importFrom abind       adrop
#' @importFrom AnnotationDbi mapIds
#' @importFrom assertive   assert_all_are_dirs   assert_all_are_existing_files
#' @importFrom assertive   assert_all_are_finite
#' @importFrom assertive   assert_all_are_greater_than
#' @importFrom assertive   assert_all_are_in_closed_range
#' @importFrom assertive   assert_all_are_in_range
#' @importFrom assertive   assert_all_are_less_than_or_equal_to
#' @importFrom assertive   assert_all_are_non_missing_nor_empty_character
#' @importFrom assertive   assert_all_are_not_matching_fixed
#' @importFrom assertive   assert_all_are_true
#' @importFrom assertive   assert_all_are_whole_numbers
#' @importFrom assertive   assert_any_are_not_nan
#' @importFrom assertive   assert_are_identical     assert_are_same_length
#' @importFrom assertive   assert_engine
#' @importFrom assertive   assert_has_names         assert_has_no_duplicates
#' @importFrom assertive   assert_is_all_of         assert_is_any_of
#' @importFrom assertive   assert_is_a_bool
#' @importFrom assertive   assert_is_a_number       assert_is_a_string
#' @importFrom assertive   assert_is_character      assert_is_data.frame
#' @importFrom assertive   assert_is_function       assert_is_identical_to_true
#' @importFrom assertive   assert_is_list
#' @importFrom assertive   assert_is_logical        assert_is_non_empty
#' @importFrom assertive   assert_is_non_scalar     assert_is_not_null
#' @importFrom assertive   assert_is_matrix         assert_is_numeric
#' @importFrom assertive   assert_is_of_length      assert_is_subset
#' @importFrom assertive   false                    get_name_in_parent
#' @importFrom assertive   has_names                has_no_duplicates
#' @importFrom assertive   is_empty_character   is_scalar   is_existing_file
#' @importFrom colorspace  sequential_hcl
#' @importFrom data.table  as.data.table   copy   data.table   dcast
#' @importFrom data.table  fread   fwrite
#' @importFrom data.table  is.data.table    .N   rbindlist   setkeyv   set
#' @importFrom data.table  setnames    setorderv   tstrsplit   :=   .SD   .I
#' @importFrom edgeR       filterByExpr
#' @import     ggplot2
#' @importFrom ggrepel     geom_text_repel
#' @importFrom graphics    lines    pie    title
#' @importFrom grDevices   hcl
#' @importFrom limma       contrasts.fit    duplicateCorrelation   eBayes
#' @importFrom limma       lmFit    makeContrasts   removeBatchEffect    voom
#' @importFrom grid        grid.layout  grid.draw   grid.newpage
#' @importFrom gridExtra   arrangeGrob   grid.arrange
#' @importFrom magrittr    %>%   %<>%   add   and   equals
#' @importFrom magrittr    divide_by   extract   extract2   is_in
#' @importFrom magrittr    multiply_by   subtract
#' @importFrom magrittr    set_colnames   set_names   set_rownames
#' @importFrom matrixStats   rowAnys   colAnys   colWeightedMeans
#' @importFrom matrixStats   rowAlls   rowSds
#' @importFrom methods       as    is
#' @importFrom MultiAssayExperiment  colData       colData<-
#' @importFrom MultiAssayExperiment  experiments   experiments<-
#' @importFrom MultiAssayExperiment  MultiAssayExperiment
#' @importFrom parallel      detectCores
#' @importFrom readxl        read_excel   excel_sheets
#' @importFrom rlang         as_string   enquo   eval_tidy   enexpr
#' @importFrom rlang         expr_text   quo_name   as_name   quo_is_null
#' @importFrom RCurl         url.exists
#' @importFrom Rsubread      featureCounts
#' @importFrom R.utils       gunzip
#' @importFrom S4Vectors     DataFrame   metadata   metadata<-
#' @importFrom seqinr        read.fasta
#' @importFrom stats         aggregate    approxfun   as.formula
#' @importFrom stats         dist   hclust   lowess   median   model.matrix
#' @importFrom stats         na.exclude   p.adjust qnorm   rnorm   sd
#' @import     stringi
#' @importFrom tidyr         separate_rows
#' @importFrom tools         file_ext   file_path_sans_ext
#' @importFrom SummarizedExperiment  assays  assays<-  assayNames  assayNames<-
#' @importFrom SummarizedExperiment  SummarizedExperiment
#' @importFrom SummarizedExperiment  rowData   rowData<-   colData   colData<-
#' @importFrom SummarizedExperiment  makeSummarizedExperimentFromExpressionSet
#' @importFrom utils         adist   capture.output   count.fields
#' @importFrom utils         download.file    getFromNamespace
#' @importFrom utils         installed.packages    unzip
NULL

utils::globalVariables(c('.'))
utils::globalVariables(c('subgroup', 'sample_id'))
utils::globalVariables(c('feature_id', 'feature_name'))
utils::globalVariables(c('pca1', 'pca2'))
utils::globalVariables('value')

