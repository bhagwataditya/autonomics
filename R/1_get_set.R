
#===============================================================================
#' Get/set analysis
#' @param object SummarizedExperiment
#' @param value list
#' @return analysis details (get) or updated object (set)
#' @rdname analysis
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' analysis(object)
#' @export
setGeneric("analysis", function(object) standardGeneric("analysis"))

#' @rdname analysis
#' @export
setMethod(
    "analysis",
    signature("SummarizedExperiment"),
    function(object) metadata(object)$analysis)

#' @rdname analysis
#' @export
setGeneric("analysis<-", function(object, value)  standardGeneric("analysis<-"))

#' @rdname analysis
setReplaceMethod(
    "analysis",
    signature("SummarizedExperiment", "list"),
    function(object, value){
        metadata(object)$analysis <- value
        object})

#' @title Get/Set exprs
#' @description Get / Set exprs matrix
#' @param object SummarizedExperiment, ExpressionSet, EList
#' @param value ratio matrix (features x samples)
#' @return exprs matrix (get) or updated object (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' exprs(object) <- 0
#' object
#' @export
setGeneric('exprs',  function(object)   standardGeneric("exprs"))

#' @rdname exprs
setMethod(
    "exprs",
    signature("SummarizedExperiment"),
    function(object)   assays(object)$exprs)

#' @rdname exprs
#' @export
setGeneric('exprs<-',   function(object, value) standardGeneric("exprs<-"))

#' @rdname exprs
setReplaceMethod(
    "exprs",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$exprs <- value
        object })

#' @rdname exprs
setReplaceMethod(
    "exprs",
    signature("SummarizedExperiment", "numeric"),
    function(object, value){
        assays(object)$exprs[] <- value
        object })


#==============================================================================
#' @title Get/Set fdata
#' @description Get/Set feature data
#' @param object SummarizedExperiment, eSet, or EList
#' @param value data.frame
#' @return feature dataframe (get) or updated object (set)
#' @examples
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' head(fdata(object)) # Getter
#' fdata(object) %<>% cbind(z=1)
#' head(fdata(object)) # Setter
#' @export
setGeneric('fdata',   function(object)   standardGeneric('fdata'))

#' @rdname fdata
setMethod('fdata',  signature('SummarizedExperiment'),
    function(object) as(rowData(object), "data.frame"))
    # do not use as.data.frame !
    # that somehow somewhere performs a check.names operation!

#(1) as.data.frame(object@elementMetadata) doesn't handle check.names correctly!
#(2) rowData returns a DataFrame (which is not generic to EList objects)

#' @rdname fdata
#' @export
setGeneric('fdata<-',   function(object, value)  standardGeneric('fdata<-'))

#' @rdname fdata
setReplaceMethod(
    'fdata',
    signature('SummarizedExperiment', 'data.frame'),
    function(object, value){
        rowData(object) <- DataFrame(value, check.names = FALSE)
        object })


#==============================================================================
#' Get fvar levels
#' @param  object  SummarizedExperiment
#' @param  fvar    feature variable
#' @return fvar values
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' head(flevels(object, 'feature_name'))
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
#' object <- read_proteingroups(file)
#' fnames(object) %<>% paste0('PG', .)
#' object
#' @rdname fnames
#' @export
setGeneric("fnames", function(object)   standardGeneric("fnames"))

#' @rdname fnames
setMethod(
    "fnames",
    signature("SummarizedExperiment"),
    function(object)   rownames(object))

#' @rdname fnames
#' @export
setGeneric(
    "fnames<-",
    function(object, value)   standardGeneric("fnames<-"))

#' @rdname fnames
setReplaceMethod(
    "fnames",
    signature("SummarizedExperiment", "character"),
    function(object, value){  rownames(object) <- value
        object})


#==============================================================================
#' @title Get fvalues
#' @description Get fvar values
#' @param  object  SummarizedExperiment
#' @param  fvar    feature variable
#' @return fvar values
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' head(fvalues(object, 'feature_name'))
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
#' object <- read_proteingroups(file)
#' fvars(object)[1] %<>% paste0('1')
#' fvars(object)[1]
#' @rdname fvars
#' @export
setGeneric("fvars", function(object)   standardGeneric("fvars"))

#' @rdname fvars
setMethod(
    "fvars",
    signature("SummarizedExperiment"),
    function(object) names(rowData(object)))

#' @rdname fvars
#' @export
setGeneric("fvars<-", function(object, value)  standardGeneric("fvars<-") )

#' @rdname fvars
setReplaceMethod(
    "fvars",
    signature("SummarizedExperiment", "character"),
    function(object, value){ names(rowData(object)) <- value
                            object })



#==============================================================================
#' @title Get/Set sdata
#' @description Get/Set sample data
#' @param object SummarizedExperiment, eSet, or EList
#' @param value dataframe
#' @return sample dataframe (get) or updated object (set)
#' @examples
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' head(sdata(object))
#' head(sdata(object) %<>% cbind(z=1))
#' @rdname sdata
#' @export
setGeneric('sdata',
    function(object){ standardGeneric('sdata')})

#' @rdname sdata
setMethod('sdata', signature('SummarizedExperiment'),
    function(object){
        as(colData(object), "data.frame") })
        # !! as.data.frame somehow somewhere performs a check.names

#' @rdname sdata
#' @export
setGeneric('sdata<-', function(object, value)  standardGeneric('sdata<-'))

#' @rdname sdata
setReplaceMethod(
    'sdata',
    signature('SummarizedExperiment', 'data.frame'),
    function(object, value){
        colData(object) <- DataFrame(value, check.names = FALSE)
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
#' object <- read_proteingroups(file)
#' head(snames(object))
#' head(snames(object) %<>% paste0('SAMPLE_', .))
#' @rdname snames
#' @export
setGeneric("snames",  function(object)   standardGeneric("snames"))

#' @rdname snames
setMethod(
    'snames',
    signature("SummarizedExperiment"),
    function(object)   colnames(object))

#' @rdname snames
#' @export
setGeneric("snames<-", function(object, value)  standardGeneric("snames<-"))

#' @rdname snames
setReplaceMethod(
    "snames",
    signature("SummarizedExperiment", "character"),
    function(object, value){
        colnames(object)  <- value
        object })


#=========================================================
#' @title Get slevels
#' @description Get svar levels
#' @param object SummarizedExperiment, eSet, or eList
#' @param svar sample var (character)
#' @return svar values (character)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
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
#' object <- read_proteingroups(file)
#' svalues(object, 'subgroup')
#' subgroup_values(object)
#' @rdname svalues
#' @export
svalues <- function(object, svar){
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
setGeneric(
    'svalues<-',
    function(object, svar, value)  standardGeneric('svalues<-'))

#' @rdname svalues
setReplaceMethod(
    'svalues',
    signature('SummarizedExperiment', 'character', "ANY"),
    function(object, svar, value){
        colData(object)[svar] <- value
        object })


#=========================================================================
#' @title Get/Set svars
#' @description Get/Set sample variables
#' @param object SummarizedExperiment
#' @param value string fector with variable names
#' @return sample variable names (get) or updated eSet (set)
#' @examples
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' svars(object)[1]
#' (svars(object)[1] %<>% paste0('1'))
#' @rdname svars
#' @export
setGeneric("svars", function(object) standardGeneric("svars") )

#' @rdname svars
setMethod(
    "svars",
    signature("SummarizedExperiment"),
    function(object)   names(colData((object))))

#' @rdname svars
#' @export
setGeneric("svars<-", function(object, value)  standardGeneric("svars<-") )

#' @rdname svars
setReplaceMethod(
    "svars",
    signature("SummarizedExperiment", "character"),
    function(object, value){
        names(colData(object)) <- value
        object
    })



