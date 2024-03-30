#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom magrittr %<>% 
#' @export
magrittr::`%<>%`

#' @importFrom magrittr extract
#' @export
magrittr::extract

#' @importFrom data.table data.table
#' @export
data.table::data.table


#' @importFrom abind                  adrop
#' @importFrom BiocFileCache          BiocFileCache   bfcquery   bfcadd   bfcrpath
#' @importFrom BiocGenerics           cbind 
#' @importFrom codingMatrices         code_control  code_diff   code_deviation
#' @importFrom codingMatrices         code_helmert  contr.diff
#' @importFrom colorspace  sequential_hcl
#' @importFrom data.table  as.data.table   copy
#' @importFrom data.table  dcast  dcast.data.table
#' @importFrom data.table  fread   fwrite  is.data.table
#' @importFrom data.table  melt.data.table   merge.data.table   
#' @importFrom data.table  .N   rbindlist   setkeyv   set
#' @importFrom data.table  setnames    setorderv   tstrsplit   :=   .SD   .I
#' @importFrom edgeR       filterByExpr
#' @import     ggplot2
#' @importFrom ggforce     facet_wrap_paginate
#' @importFrom ggrepel     geom_text_repel    geom_label_repel
#' @importFrom graphics    lines    pie    title
#' @importFrom grDevices   hcl    dev.off   pdf
#' @importFrom limma       contrasts.fit    duplicateCorrelation   eBayes
#' @importFrom limma       lmFit    makeContrasts   removeBatchEffect    
#' @importFrom limma       voom   vennDiagram
#' @importFrom graphics    abline   axis   box   image   layout   par
#' @importFrom grid        grid.layout  grid.draw   grid.newpage
#' @importFrom gridExtra   arrangeGrob   grid.arrange
#' @importFrom magrittr    add   and   equals
#' @importFrom magrittr    divide_by   extract2   is_in
#' @importFrom magrittr    multiply_by   subtract
#' @importFrom magrittr    set_colnames   set_names   set_rownames
#' @importFrom matrixStats   rowAnys   colAnys   colWeightedMeans
#' @importFrom matrixStats   rowAlls   rowMaxs   rowMins   rowSds   rowVars
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
utils::globalVariables('N')         # data.table constructs