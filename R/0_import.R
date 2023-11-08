#' @importFrom abind                 adrop
#' @importFrom assertive.base        assert_all_are_true
#' @importFrom assertive.base        assert_any_are_true
#' @importFrom assertive.base        assert_are_identical
#' @importFrom assertive.base        assert_engine
#' @importFrom assertive.base        assert_is_identical_to_true
#' @importFrom assertive.base        false
#' @importFrom assertive.base        get_name_in_parent
#' @importFrom assertive.files       assert_all_are_dirs
#' @importFrom assertive.files       assert_all_are_existing_files
#' @importFrom assertive.files       is_existing_file
#' @importFrom assertive.numbers     assert_all_are_finite
#' @importFrom assertive.numbers     assert_all_are_greater_than
#' @importFrom assertive.numbers     assert_all_are_in_closed_range
#' @importFrom assertive.numbers     assert_all_are_in_range
#' @importFrom assertive.numbers     assert_all_are_less_than_or_equal_to
#' @importFrom assertive.numbers     assert_all_are_whole_numbers
#' @importFrom assertive.numbers     assert_any_are_not_nan
#' @importFrom assertive.sets        assert_is_subset
#' @importFrom assertive.sets        is_subset
#' @importFrom BiocGenerics cbind 
#' @importFrom colorspace  sequential_hcl
#' @importFrom data.table  as.data.table   copy  data.table  
#' @importFrom data.table  dcast  dcast.data.table
#' @importFrom data.table  fread   fwrite  is.data.table
#' @importFrom data.table  melt.data.table   .N   rbindlist   setkeyv   set
#' @importFrom data.table  setnames    setorderv   tstrsplit   :=   .SD   .I
#' @importFrom edgeR       filterByExpr
#' @import     ggplot2
#' @importFrom ggrepel     geom_text_repel
#' @importFrom graphics    lines    pie    title
#' @importFrom grDevices   hcl
#' @importFrom limma       contrasts.fit    duplicateCorrelation   eBayes
#' @importFrom limma       lmFit    makeContrasts   removeBatchEffect    
#' @importFrom limma       voom   vennDiagram
#' @importFrom graphics    layout
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
#' @importFrom R.utils       gunzip
#' @importFrom S4Vectors     DataFrame   metadata   metadata<-
#' @importFrom stats         aggregate    approxfun   as.formula   dist   hclust
#' @importFrom stats         lowess   lm  median   model.matrix   wilcox.test
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

