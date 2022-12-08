
#===============================================================================

#' Get/set analysis
#' @param object SummarizedExperiment
#' @param value list
#' @return analysis details (get) or updated object (set)
#' @rdname analysis
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file, plot = FALSE)
#' analysis(object)
#' @export
setGeneric("analysis", function(object) standardGeneric("analysis"))


#' @rdname analysis
#' @export
setMethod("analysis", signature("SummarizedExperiment"),
function(object) metadata(object)$analysis)


#' @rdname analysis
#' @export
setGeneric("analysis<-", function(object, value)  standardGeneric("analysis<-"))


#' @rdname analysis
setReplaceMethod("analysis", signature("SummarizedExperiment", "list"),
function(object, value){
    metadata(object)$analysis <- value
    object})


#' @title Get/Set expr values
#' @description Get/Set value matrix
#' @param object SummarizedExperiment
#' @param value ratio matrix (features x samples)
#' @return value matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' values(object)[1:3, 1:3]
#' values(object) <- 0
#' values(object)[1:3, 1:3]
#' @export
setGeneric('values',  function(object)   standardGeneric("values"))


#' @rdname values
setMethod("values", signature("SummarizedExperiment"),
function(object)   assays(object)[[1]])


#' @rdname values
#' @export
setGeneric('values<-',   function(object, value) standardGeneric("values<-"))


#' @rdname values
setReplaceMethod("values", signature("SummarizedExperiment", "matrix"),
function(object, value){
    assays(object)[[1]] <- value
    object })


#' @rdname values
setReplaceMethod("values", signature("SummarizedExperiment", "numeric"),
function(object, value){
    assays(object)[[1]][] <- value
    object })


#==============================================================================

#' Get fvar levels
#' @param  object  SummarizedExperiment
#' @param  fvar    feature variable
#' @return fvar values
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' head(flevels(object, 'feature_id'))
#' @export
flevels <- function(object, fvar){
    object %>%
    fvalues(fvar) %>%
    (function(x)if (is.factor(x)) levels(x) else unique(x))
}



#==============================================================================

#' @title Get/Set fnames
#' @description Get/Set feature names
#' @param object SummarizedExperiment, eSet, or EList
#' @param value character vector with feature names
#' @return feature name vector (get) or updated object (set)
#' @examples
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' fnames(object) %<>% paste0('protein_', .)
#' object
#' @rdname fnames
#' @export
setGeneric("fnames", function(object)   standardGeneric("fnames"))


#' @rdname fnames
setMethod("fnames", signature("SummarizedExperiment"),
function(object)   rownames(object))


#' @rdname fnames
#' @export
setGeneric("fnames<-",
function(object, value)   standardGeneric("fnames<-"))


#' @rdname fnames
setReplaceMethod("fnames", signature("SummarizedExperiment", "character"),
function(object, value){  
    fdata(object)$feature_id <- NULL
    rownames(object) <- value
    fdata(object)$feature_id <- value
    fdata(object)$feature_id %<>% factor()
    object
})



#==============================================================================

#' @title Get fvalues
#' @description Get fvar values
#' @param  object  SummarizedExperiment
#' @param  fvar    feature variable
#' @return fvar values
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' head(fvalues(object, 'feature_id'))
#' fvalues(object, NULL)
#' @export
fvalues <- function(object, fvar){

    # Return NULL output for NULL input
    if (is.null(fvar)) return(NULL)

    # Assert that fvar is present
    assert_is_subset(fvar, fvars(object))

    # Extract and return
    fdata(object)[[fvar]]
}

#' Get feature id variable
#' @param object SummarizedExperiment
#' @noRd
fid_var <- function(object) 'feature_id'


#' @rdname fid_var
#' @noRd
fid_values <- function(object) fvalues(object, 'feature_id')



#==============================================================================

#' @title Get/Set fvars
#' @description Get/Set feature variables
#' @param object SummarizedExperiment
#' @param value character vector with feature variables
#' @return feature variables vector (get) or updated object (set)
#' @examples
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' fvars(object)[1] %<>% paste0('1')
#' fvars(object)[1]
#' @rdname fvars
#' @export
setGeneric("fvars", function(object)   standardGeneric("fvars"))


#' @rdname fvars
setMethod("fvars", signature("SummarizedExperiment"),
function(object) names(rowData(object)))


#' @rdname fvars
#' @export
setGeneric("fvars<-", function(object, value)  standardGeneric("fvars<-") )


#' @rdname fvars
setReplaceMethod("fvars", signature("SummarizedExperiment", "character"),
function(object, value){ names(rowData(object)) <- value
                        object })

#==============================================================================

#' Get/Set sample/feature data
#' @param object SummarizedExperiment/MultiAssayExperiment
#' @param value data.frame/data.table
#' @return data.frame/data.table (get) or updated object (set)
#' @examples
#' # Read data
#'     require(magrittr)
#'     file <- download_data('billing16.proteingroups.txt')
#'     object <- read_maxquant_proteingroups(file)
#' # sdata/fdata
#'     sdata(object)[1:3, ]
#'     fdata(object)[1:3, ]
#'     sdata(object) %<>% cbind(a=1)
#'     fdata(object) %<>% cbind(a=1)
#'     sdata(object)[1:3, ]
#'     fdata(object)[1:3, ]
#' # sdt/fdt
#'     sdt(object)[1:3, ]
#'     fdt(object)[1:3, ]
#'     sdt(object) %<>% cbind(b=1)
#'     fdt(object) %<>% cbind(b=1)
#'     sdt(object)
#'     fdt(object)
#' @rdname fdata
#' @export
setGeneric('fdata',  function(object)  standardGeneric('fdata'))               # fdata

#' @rdname fdata
#' @export
setGeneric('sdata',  function(object)  standardGeneric('sdata'))               # sdata

#' @rdname fdata
#' @export
setGeneric('fdt',    function(object)  standardGeneric('fdt'))                 # fdt

#' @rdname fdata
#' @export
setGeneric('sdt',    function(object)  standardGeneric('sdt'))                 # sdt

#' @rdname fdata
setMethod('fdata',  signature('SummarizedExperiment'),                         # fdata se
function(object)  as(rowData(object), "data.frame"))
    # !! as.data.frame somehow somewhere performs a check.names

#' @rdname fdata
setMethod('sdata', signature('SummarizedExperiment'),                          # sdata se
function(object)  as(colData(object), "data.frame")) 
    # !! as.data.frame somehow somewhere performs a check.names

#' @rdname fdata
setMethod('sdata', signature('MultiAssayExperiment'),                          # sdata mae
function(object)  as(colData(object), "data.frame"))

#' @rdname fdata
setMethod('fdt',  signature('SummarizedExperiment'),                           # fdt se
function(object)  data.table(data.frame(rowData(object), check.names = FALSE)))

#' @rdname fdata
setMethod('sdt',  signature('SummarizedExperiment'),                           # sdt se
function(object)  data.table(data.frame(colData(object), check.names = FALSE)))

#' @rdname fdata
setMethod('sdt', signature('MultiAssayExperiment'),                            # sdt mae
function(object)  data.table(data.frame(colData(object))))

#' @rdname fdata
#' @export
setGeneric('fdata<-', function(object, value)  standardGeneric('fdata<-'))     # fdata<-

#' @rdname fdata
#' @export
setGeneric('sdata<-', function(object, value)  standardGeneric('sdata<-'))     # sdata<-

#' @rdname fdata
#' @export
setGeneric('fdt<-',   function(object, value)  standardGeneric('fdt<-'))       # fdt<-

#' @rdname fdata
#' @export
setGeneric('sdt<-',   function(object, value)  standardGeneric('sdt<-'))       # sdt<-

#' @rdname fdata
setReplaceMethod('fdata', signature('SummarizedExperiment', 'data.frame'),     # fdata<- se df
function(object, value){
    rowData(object) <- DataFrame(value, check.names = FALSE)
    object })

#' @rdname fdata
setReplaceMethod('sdata',  signature('SummarizedExperiment', 'data.frame'),    # sdata<- se df
function(object, value){
    colData(object) <- DataFrame(value, check.names = FALSE)
    object })

#' @rdname fdata
setReplaceMethod('sdata', signature('SummarizedExperiment', 'DataFrame'),      # sdata<- se DF
function(object, value){
    colData(object) <- value
    object })

#' @rdname fdata
setReplaceMethod('sdata', signature('MultiAssayExperiment', 'data.frame'),     # sdata<- mae df
function(object, value){
    colData(object) <- DataFrame(value, check.names = FALSE)
    object })

#' @rdname fdata
setReplaceMethod('sdata', signature('MultiAssayExperiment', 'DataFrame'),      # sdata<- mae DF
function(object, value){
    colData(object) <- value
    object })

#' @rdname fdata
setReplaceMethod('fdt', signature('SummarizedExperiment', 'data.table'),       # fdt<- se dt
function(object, value){
    rowData(object) <- DataFrame(value, check.names = FALSE, row.names = value$feature_id)
    object })

#' @rdname fdata
setReplaceMethod('sdt', signature('SummarizedExperiment', 'data.table'),       # sdt<- se dt
function(object, value){
    colData(object) <- DataFrame(value, check.names = FALSE, row.names = value$sample_id)
    object })

#' @rdname fdata
setReplaceMethod('sdt', signature('MultiAssayExperiment', 'data.table'),       # sdt<- mae dt
function(object, value){
    colData(object) <- DataFrame(value, check.names = FALSE, row.names = value$sample_id)
    object })



#=====================================================================

#' @title Get/Set snames
#' @description Get/Set sample names
#' @param object SummarizedExperiment
#' @param value string vector with sample names
#' @return sample names vector (get) or updated eSet (set)
#' @examples
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' head(snames(object))
#' head(snames(object) %<>% paste0('SAMPLE_', .))
#' @rdname snames
#' @export
setGeneric("snames",  function(object)   standardGeneric("snames"))


#' @rdname snames
setMethod('snames',
    signature("SummarizedExperiment"),
    function(object)   colnames(object))

#' @rdname snames
#' @export
setGeneric("snames<-", function(object, value)  standardGeneric("snames<-"))


#' @rdname snames
setReplaceMethod("snames", signature("SummarizedExperiment", "character"),
function(object, value){
    colnames(object)  <- value
    object$sample_id <- value
    object })



#=========================================================

#' @title Get slevels
#' @description Get svar levels
#' @param object SummarizedExperiment, eSet, or eList
#' @param svar sample var (character)
#' @return svar values (character)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' slevels(object, 'subgroup')
#' subgroup_levels(object)
#' @rdname slevels
#' @export
slevels <- function(object, svar){
    object %>%
    svalues(svar) %>%
    (function(x) if (is.factor(x)) levels(x) else sort(unique(x)))
}

#' @rdname slevels
#' @export
subgroup_levels <- function(object){
    slevels(object, 'subgroup')
}


#=========================================================

#' @title Get/Set svalues
#' @description Get/Set svar values
#' @param object SummarizedExperiment
#' @param svar   sample var (character)
#' @param value  value vector
#' @return character vector (get) or SummarizedExperiment (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' svalues(object, 'subgroup')
#' subgroup_values(object)
#' @rdname svalues
#' @export
svalues <- function(object, svar){
    if (is.null(svar)) return(NULL)
    sdata(object)[[svar]]
}

#' @rdname svalues
#' @export
subgroup_values <- function(object){
    svalues(object, 'subgroup')
}

#' @rdname svalues
#' @export
sampleid_values <- function(object){
    svalues(object, 'sample_id')
}


# Set
#====
#' @rdname svalues
#' @export
setGeneric('svalues<-',
function(object, svar, value)  standardGeneric('svalues<-'))


#' @rdname svalues
setReplaceMethod('svalues', signature('SummarizedExperiment','character',"ANY"),
function(object, svar, value){
    colData(object)[svar] <- value
    object })



#=========================================================================

#' @title Get/Set svars
#' @description Get/Set sample variables
#' @param object SummarizedExperiment
#' @param value string fector with variable names
#' @return sample variable names (get) or updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' svars(object)[1]
#'(svars(object)[1] %<>% paste0('1'))
#' @rdname svars
#' @export
setGeneric("svars", function(object) standardGeneric("svars") )

#' @rdname svars
setMethod("svars", signature("SummarizedExperiment"),
function(object)   names(colData((object))))

#' @rdname svars
setMethod("svars", signature("MultiAssayExperiment"),
function(object)   names(colData((object))))

#' @rdname svars
#' @export
setGeneric("svars<-", function(object, value)  standardGeneric("svars<-") )

#' @rdname svars
setReplaceMethod("svars", signature("SummarizedExperiment", "character"),
function(object, value){
    names(colData(object)) <- value
    object
})

#' @rdname svars
setReplaceMethod("svars", signature("MultiAssayExperiment", "character"),
function(object, value){
    names(colData(object)) <- value
    object
})


