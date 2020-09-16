
#===============================================================================
#' Get/set analysis
#' @param object SummarizedExperiment
#' @param value list
#' @return analysis details (get) or updated object (set)
#' @rdname analysis
#' @examples
#' file <- download_data('stemcells.proteinGroups.txt')
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


#==============================================================================
#' Get/Set contrast definitions
#' @param object SummarizedExperiment
#' @param value named string vector (see examples)
#' @return updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('stemcells.proteinGroups.txt')
#' object <- read_proteingroups(file)
#' contrastdefs1 <- c(EM_E = 'EM_E', BM_E = 'BM_E', BM_EM = 'BM_EM')
#' contrastdefs(object) <-     contrastdefs1   # conventional setter
#' object %>% set_contrastdefs(contrastdefs1)  # piping       setter
#' contrastdefs(object)                        # getter
#' @rdname contrastdefs
#' @export
setGeneric("contrastdefs",   function(object)   standardGeneric("contrastdefs"))

#' @rdname contrastdefs
#' @export
setMethod(
    "contrastdefs",
    signature("SummarizedExperiment"),
    function(object) metadata(object)$contrastdefs )

#' @rdname contrastdefs
#' @export
setGeneric(
    "contrastdefs<-",
    function(object, value)  standardGeneric("contrastdefs<-") )

#' @rdname contrastdefs
setReplaceMethod(
    "contrastdefs",
    signature("SummarizedExperiment", "character"),
    function(object, value){ metadata(object)$contrastdefs <- value; object})

#' @rdname contrastdefs
setReplaceMethod(
    "contrastdefs",
    signature("SummarizedExperiment", "NULL"),
    function(object, value){object})

#' @rdname contrastdefs
#' @export
set_contrastdefs <- function(object, value){
    contrastdefs(object) <- value
    object
}

#=========================================================================

#' @title Get/Set counts
#' @description Get / Set counts matrix
#' @param object SummarizedExperiment
#' @param value count matrix (features x samples)
#' @return count matrix (get) or updated object (set)
#' @examples
#' file <- download_data('stemcells.rnacounts.txt')
#' object <- read_counts(file)
#' counts(object) <- exprs(object)
#' counts(object)[1:3, 1:3]
#' @rdname counts
#' @export
setGeneric('counts',   function(object)   standardGeneric("counts"))

#' @rdname counts
setMethod(
    "counts",
    signature("SummarizedExperiment"),
    function(object)   assays(object)$counts)

#' @rdname counts
#' @export
setGeneric('counts<-',   function(object, value)   standardGeneric("counts<-"))

#' @rdname counts
setReplaceMethod(
    "counts",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$counts <- value
        object })

#' @rdname counts
setReplaceMethod(
    "counts",
    signature("SummarizedExperiment", "numeric"),
    function(object, value){
        assays(object)$counts[] <- value
        object })


#' @title Get/Set exprs
#' @description Get / Set exprs matrix
#' @param object SummarizedExperiment, ExpressionSet, EList
#' @param value ratio matrix (features x samples)
#' @return exprs matrix (get) or updated object (set)
#' @examples
#' file <- download_data('stemcells.proteinGroups.txt')
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
#' file <- download_data('stemcells.proteinGroups.txt')
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
#' file <- download_data('stemcells.proteinGroups.txt')
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
#' file <- download_data('stemcells.proteinGroups.txt')
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
#' file <- download_data('stemcells.proteinGroups.txt')
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
#' file <- download_data('stemcells.proteinGroups.txt')
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
#'@title Get/set is_imputed
#'@description Get/Set is_imputed
#'@param object SummarizedExperiment
#'@param value matrix
#'@return matrix (get) or updated object (set)
#'@examples
#' file <- download_data('stemcells.proteinGroups.txt')
#' object <- read_proteingroups(file)
#' sum(is_imputed(object))
#' @rdname is_imputed
#' @export
setGeneric("is_imputed",  function(object) standardGeneric("is_imputed") )

#' @rdname is_imputed
setMethod("is_imputed", signature("SummarizedExperiment"),  function(object){
    if ('is_imputed' %in% names(assays(object))){
        assays(object)$is_imputed
    } else {
        matrix(FALSE, nrow = nrow(object), ncol = ncol(object),
                dimnames = dimnames(object))
    }
})

#' @rdname is_imputed
#' @export
setGeneric(
    "is_imputed<-",
    function(object, value)  standardGeneric("is_imputed<-") )

#' @rdname is_imputed
setReplaceMethod(
    "is_imputed",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$is_imputed <- value; object})

#' @rdname is_imputed
setReplaceMethod(
    "is_imputed",
    signature("SummarizedExperiment", "NULL"),
    function(object, value){object})

#==============================================================================
#' @title Get/set limma results
#' @description Get/Set limma results
#' @param object SummarizedExperiment
#' @param value list
#' @return limma results (get) or updated object (set)
#' @export
setGeneric("limma", function(object)   standardGeneric("limma") )

#' @rdname limma
setMethod(
    "limma",
    signature("SummarizedExperiment"),
    function(object) S4Vectors::metadata(object)$limma )

#' @rdname limma
#' @export
setGeneric("limma<-", function(object, value)  standardGeneric("limma<-") )

#' @rdname limma
setReplaceMethod(
    "limma",
    signature("SummarizedExperiment", "array"),
    function(object, value){ metadata(object)$limma <- value; object})

#' @rdname limma
setReplaceMethod(
    "limma",
    signature("SummarizedExperiment", "NULL"),
    function(object, value){object})



#=========================================================================

#' @title Get/Set occupancies
#' @description Get / Set occupancies matrix
#' @param object SummarizedExperiment
#' @param value occupancy matrix (features x samples)
#' @return occpuancy matrix (get) or updated object (set)
#' @examples
#' proteinfile <- download_data('differentiation.proteinGroups.txt')
#' phosphofile <- download_data('differentiation.phosphoSites.txt')
#' object <- read_phosphosites(phosphofile, proteinfile)
#' exprs(object)[1:3, 1:3]
#' occupancies(object)[1:3, 1:3]
#' @rdname occupancies
#' @export
setGeneric('occupancies', function(object)   standardGeneric("occupancies"))

#' @rdname occupancies
setMethod(
    "occupancies",
    signature("SummarizedExperiment"),
    function(object)   assays(object)$occupancies)


#' @rdname occupancies
#' @export
setGeneric(
    'occupancies<-',
    function(object, value) standardGeneric("occupancies<-"))

#' @rdname occupancies
setReplaceMethod(
    "occupancies",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$occupancies <- value
        object })

#' @rdname occupancies
setReplaceMethod(
    "occupancies",
    signature("SummarizedExperiment", "numeric"),
    function(object, value){
        assays(object)$occupancies[] <- value
        object })


#==============================================================================
#' @title Get/Set sdata
#' @description Get/Set sample data
#' @param object SummarizedExperiment, eSet, or EList
#' @param value dataframe
#' @return sample dataframe (get) or updated object (set)
#' @examples
#' require(magrittr)
#' file <- download_data('stemcells.proteinGroups.txt')
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
#' file <- download_data('stemcells.proteinGroups.txt')
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
#' file <- download_data('stemcells.proteinGroups.txt')
#' object <- read_proteingroups(file)
#' slevels(object, 'subgroup')
#' subgroup_levels(object)
#' @rdname slevels
#' @export
slevels <- function(object, svar){
    object %>%
    svalues(svar) %>%
    (function(x) if (is.factor(x)) levels(x) else unique(x))
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
#' file <- download_data('stemcells.proteinGroups.txt')
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
#' file <- download_data('stemcells.proteinGroups.txt')
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



#=========================================================================
#' @title Get/Set weights
#' @description Get/Set weight matrix
#' @param object SummarizedExperiment
#' @param value ratio matrix (features x samples)
#' @param ... addtional params
#' @return weight matrix (get) or updated object (set)
#' @examples
#' file <- download_data('stemcells.proteinGroups.txt')
#' object <- read_proteingroups(file)
#' weights(object)[1:3, 1:2]
#' weights(object) <- 1; weights(object)[1:3, 1:2]
#' @rdname weights
#' @export
setGeneric('weights', function(object)   standardGeneric("weights"))

#' @rdname weights
setMethod(
    "weights",
    signature("SummarizedExperiment"),
    function(object)   assays(object)$weights)


#' @rdname weights
#' @export
setGeneric('weights<-', function(object, value) standardGeneric("weights<-"))

#' @rdname weights
setReplaceMethod(
    "weights",
    signature("SummarizedExperiment", "matrix"),
    function(object, value){
        assays(object)$weights <- value
        object })

#' @rdname weights
setReplaceMethod(
    "weights",
    signature("SummarizedExperiment", "numeric"),
    function(object, value){
        if (!'weights' %in% names(assays(object))){
            assays(object)$weights <- matrix(
            1, nrow=nrow(object), ncol=ncol(object), dimnames=dimnames(object))
        }
        assays(object)$weights[] <- value
        object })

