#' @importFrom abind                  adrop
#' @importFrom assertive.base         are_identical
#' @importFrom assertive.base         assert_all_are_true   assert_any_are_true
#' @importFrom assertive.base         assert_are_identical
#' @importFrom assertive.base         assert_engine
#' @importFrom assertive.base         assert_is_identical_to_false
#' @importFrom assertive.base         assert_is_identical_to_true
#' @importFrom assertive.base         false    get_name_in_parent
#' @importFrom assertive.files        assert_all_are_dirs   
#' @importFrom assertive.files        assert_all_are_existing_files
#' @importFrom assertive.files        is_existing_file
#' @importFrom assertive.numbers      assert_all_are_in_left_open_range
#' @importFrom assertive.numbers      assert_all_are_in_range
#' @importFrom assertive.numbers      assert_all_are_less_than_or_equal_to
#' @importFrom assertive.numbers      is_in_closed_range
#' @importFrom assertive.numbers      is_greater_than   is_greater_than_or_equal_to
#' @importFrom assertive.properties   assert_are_same_length
#' @importFrom assertive.properties   assert_has_names
#' @importFrom assertive.properties   assert_has_no_duplicates
#' @importFrom assertive.properties   assert_is_of_length
#' @importFrom assertive.properties   assert_is_not_null
#' @importFrom assertive.properties   assert_is_scalar
#' @importFrom assertive.properties   has_names   has_no_duplicates
#' @importFrom assertive.properties   is_empty    is_scalar  
#' @importFrom assertive.sets         assert_are_disjoint_sets
#' @importFrom assertive.sets         assert_is_subset   is_subset
#' @importFrom assertive.strings      assert_all_are_non_missing_nor_empty_character
#' @importFrom assertive.strings      assert_any_are_matching_regex
#' @importFrom assertive.strings      is_numeric_string
#' @importFrom assertive.types        assert_is_all_of
#' @importFrom assertive.types        assert_is_a_bool
#' @importFrom assertive.types        assert_is_a_number
#' @importFrom assertive.types        assert_is_a_string
#' @importFrom assertive.types        assert_is_character
#' @importFrom assertive.types        assert_is_data.frame
#' @importFrom assertive.types        assert_is_data.table
#' @importFrom assertive.types        assert_is_factor       assert_is_formula
#' @importFrom assertive.types        assert_is_function
#' @importFrom assertive.types        assert_is_list
#' @importFrom assertive.types        assert_is_matrix
#' @importFrom assertive.types        assert_is_numeric
#' @importFrom assertive.types        is_a_number   is_a_string   is_formula  
#' @importFrom BiocFileCache          BiocFileCache   bfcquery   bfcadd   bfcrpath
#' @importFrom BiocGenerics           cbind 
#' @importFrom codingMatrices         code_control  code_diff   code_deviation
#' @importFrom codingMatrices         code_helmert  contr.diff
#' @importFrom colorspace  sequential_hcl
#' @importFrom data.table  as.data.table   copy  data.table  
#' @importFrom data.table  dcast  dcast.data.table
#' @importFrom data.table  fread   fwrite  is.data.table
#' @importFrom data.table  melt.data.table   merge.data.table   
#' @importFrom data.table  .N   rbindlist   setkeyv   set
#' @importFrom data.table  setnames    setorderv   tstrsplit   :=   .SD   .I
#' @importFrom edgeR       filterByExpr
#' @import     ggplot2
#' @importFrom ggforce     facet_wrap_paginate
#' @importFrom ggrepel     geom_text_repel
#' @importFrom graphics    lines    pie    title
#' @importFrom grDevices   hcl    dev.off   pdf
#' @importFrom limma       contrasts.fit    duplicateCorrelation   eBayes
#' @importFrom limma       lmFit    makeContrasts   removeBatchEffect    
#' @importFrom limma       voom   vennDiagram
#' @importFrom graphics    abline   axis   box   image   layout   par
#' @importFrom grid        grid.layout  grid.draw   grid.newpage
#' @importFrom gridExtra   arrangeGrob   grid.arrange
#' @importFrom magrittr    %>%   %<>%   add   and   equals
#' @importFrom magrittr    divide_by   extract   extract2   is_in
#' @importFrom magrittr    multiply_by   subtract
#' @importFrom magrittr    set_colnames   set_names   set_rownames
#' @importFrom matrixStats   rowAnys   colAnys   colWeightedMeans
#' @importFrom matrixStats   rowAlls   rowSds   rowVars
#' @importFrom methods       as    is
#' @importFrom MultiAssayExperiment  colData       colData<-
#' @importFrom MultiAssayExperiment  experiments   experiments<-
#' @importFrom MultiAssayExperiment  MultiAssayExperiment
#' @importFrom parallel      detectCores
#' @importFrom RColorBrewer  brewer.pal
#' @importFrom readxl        read_excel   excel_sheets
#' @importFrom rlang         as_string   enquo   eval_tidy   enexpr
#' @importFrom rlang         expr_text   quo_name   as_name   quo_is_null
#' @importFrom R.utils       gunzip
#' @importFrom S4Vectors     DataFrame         metadata           metadata<-
#' @importFrom stats         aggregate         approxfun          as.formula   
#' @importFrom stats         as.dist           contrasts          contrasts<-
#' @importFrom stats         cor               dist
#' @importFrom stats         contr.sum         contr.treatment
#' @importFrom stats         hclust            IQR                lowess
#' @importFrom stats         lm                median             model.matrix
#' @importFrom stats         wilcox.test       na.exclude         p.adjust
#' @importFrom stats         pchisq            qnorm              quantile
#' @importFrom stats         rnorm             sd
#' @import     stringi
#' @importFrom tidyr         separate_rows
#' @importFrom tools         file_ext   file_path_sans_ext   R_user_dir
#' @importFrom SummarizedExperiment  assays  assays<-  assayNames  assayNames<-
#' @importFrom SummarizedExperiment  SummarizedExperiment
#' @importFrom SummarizedExperiment  rowData   rowData<-   colData   colData<-
#' @importFrom SummarizedExperiment  makeSummarizedExperimentFromExpressionSet
#' @importFrom utils         adist   capture.output   count.fields
#' @importFrom utils         download.file    getFromNamespace
#' @importFrom utils         installed.packages    unzip
NULL

utils::globalVariables('patterns')  # data.table constructs
utils::globalVariables('.')         # magrittr   constructs