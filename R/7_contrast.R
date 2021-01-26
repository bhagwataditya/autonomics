#=============================================================================
#
#                   guess_sep
#                   guess_sep.SummarizedExperiment
#                   guess_sep.factor
#                   guess_sep.character
#                       has_identical_values
#                       is_max
#                           cequals
#
#=============================================================================

#' Convenient equals operator
#'
#' Performs x == y, but returns FALSE rather than NA for NA elements of x.
#' @param x numeric vector or scalar
#' @param y numeric scalar
#' @return logical vector
#' @examples
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' y <- 3
#' cequals(x, y)
#' @noRd
cequals <- function(x,y){
    result <- rep(FALSE, length(x)) %>% set_names(names(x))
    if (is.na(y)){
        result[ is.na(x)] <- TRUE
        result[!is.na(x)] <- FALSE
    } else {
        result[ is.na(x)] <- FALSE
        result[!is.na(x)] <- x[!is.na(x)] == y
    }
    result
}


#' Is maximal
#' @param x numeric vector
#' @return logical vector
#' @examples
#' x <- c(A=1,B=3,C=2,D=3, E=NA)
#' is_max(x)
#' @noRd
is_max <- function(x) cequals(x, max(x, na.rm = TRUE))


#' All elements of vector are identical
#' @param x vector
#' @return TRUE or FALSE
#' @examples
#' x <- c(2,2,1,2)
#' has_identical_values(x)
#' @noRd
has_identical_values <- function(x) length(unique(x))==1

#' Guess separator
#' @param x          character vector or SummarizedExperiment
#' @param var        svar or fvar
#' @param separators character vector: possible separators to look for
#' @param verbose    TRUE or FALSE
#' @param ...        used for proper S3 method dispatch
#' @return separator (string) or NULL (if no separator could be identified)
#' @examples
#' # charactervector
#'    x <- c('PERM_NON.R1[H/L]', 'PERM_NON.R2[H/L]', 'PERM_NON.R3[H/L]')
#'    guess_sep(x)
#'
#'    x <- c('WT untreated 1', 'WT untreated 2', 'WT treated 1')
#'    guess_sep(x)
#'
#'    x <- c('group1', 'group2', 'group3.R1')
#'    guess_sep(x)
#'
#' # SummarizedExperiment
#'    # file <- download_data('halama18.metabolon.xlsx')
#'    # object <- read_metabolon(file)
#'    # guess_sep(object)
#'
#'    # file <- download_data('billing16.proteingroups.txt')
#'    # object <- read_proteingroups(object)
#'    # guess_sep(object)
#' @export
guess_sep <- function (x, ...)  UseMethod("guess_sep", x)


#' @rdname guess_sep
#' @export
guess_sep.character <- function(
    x, separators = c('.', '_'), verbose = FALSE, ...
){
# Initialize
    . <- NULL
    sep_freqs <-Map(function(y) stri_split_fixed(x, y), separators) %>%
                lapply(function(y) vapply(y, length, integer(1)))            %>%
                extract( vapply(., has_identical_values, logical(1)))        %>%
                vapply(unique, integer(1))
# No separator detected - return NULL
    if (all(sep_freqs==1)){
        if (verbose) message(x[1],': no (consistent) separator. Returning NULL')
        return('NOSEP')   # no separator detected
    }
# Find best separator
    best_sep <- sep_freqs %>%
                extract(.!=1)  %>%
                extract(is_max(vapply(., extract, integer(1), 1)))   %>%
                names()
# Ambiguous separator - take first from tail
    if (length(best_sep)>1){
        pattern <- best_sep %>% paste0(collapse='') %>% paste0('[', ., ']')
        best_sep <- x[1] %>% stri_extract_last_regex(pattern)
    }
# Separator identified - return
    if (verbose) message("\t\tGuess sep: '", best_sep, "'")
    return(best_sep)
}


#' @rdname guess_sep
#' @export
guess_sep.factor <- function(x, ...)  guess_sep.character(levels(x))


#' @rdname guess_sep
#' @export
guess_sep.SummarizedExperiment <- function(
    x, var = 'sample_id', separators =  c('.', '_'), verbose = FALSE, ...
){
    assert_is_subset(var, c(svars(x), fvars(x)))
    (if (var %in% svars(x)) slevels(x, var) else flevels(x, var)) %>%
    guess_sep(separators = separators, verbose = verbose)
}


#=============================================================================
#
#                       nfactors
#                       split_extract
#
#=============================================================================

#' @export
#' @rdname split_extract
nfactors <- function(x, sep = guess_sep(x)){
    length(unlist(stri_split_fixed(x[1], sep)))
}

#' stri_split and extract
#' @param x string
#' @param i integer
#' @param sep string
#' @return character
#' @examples
#' require(magrittr)
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' x <- object$sample_id[1:5]
#' nfactors(x)
#' split_extract(x, 1:2)
#' split_extract(x, seq_len(nfactors(x)-1))
#' split_extract(x, nfactors(x))
#' @export
split_extract <- function(x, i, sep=guess_sep(x)){
    factors <- stri_split_fixed(x, sep)
    vapply(factors, function(y) paste0(y[i], collapse=sep), character(1))
}


#=============================================================================
#
#               merge_samplefile
#                   file_exists
#                       default_samplefile
#
#=============================================================================

# Deals properly with NULL values
# file.exists does not!
file_exists <- function(file){
    if (is.null(file))      return(FALSE)
    if (file.exists(file))  return(TRUE)
                            return(FALSE)
}

#' Default samplefile
#' @param file data file
#' @return sample file
#' @export
default_samplefile <- function(file){
    samplefile <- tools::file_path_sans_ext(file)
    samplefile %<>% paste0('.samples.txt')
    samplefile
}


#==============================================================================
#
#                       create_design
#                           single_subgroup
#                           are_factor
#
#==============================================================================

single_subgroup <- function(object){
    assert_is_subset('subgroup', svars(object))
    length(unique(object$subgroup))==1
}


are_factor <- function(df) vapply(df, is.factor, logical(1))


#' Create design
#'
#'  Create design matrix  for statistical analysis
#' @param object  SummarizedExperiment
#' @param formula formula with svars
#' @param verbose whether to message
#' @return design matrix
#' @examples
#' file <- download_data('billing16.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' create_design(object)
#'
#' object$subgroup <- 'billing16'
#' create_design(object)
#'
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' create_design(object)
#' create_design(object, ~ 0 + subgroup + Sex + T2D + age + bmi)
#' object$subgroup <- 'atkin18'
#' create_design(object)
#' @export
create_design <- function(object, formula = NULL, verbose = TRUE){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    if (is.null(formula)){
        formula <- if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup }
    assert_is_subset(all.vars(formula), svars(object))
    . <- NULL
# Ensure that subgroup vector is a factor to preserve order of levels
    object$subgroup %<>% factor() # required for '(Intercept)' -> factor1.level1
    for (var in setdiff(all.vars(formula), 'subgroup')){
        if (is.character(sdata(object)[[var]])){
            sdata(object)[[var]] %<>% factor()
        }
    }
# Create design matrix
    if (verbose)  cmessage('\t\tDesign: %s', deparse(formula))
    myDesign <- model.matrix(formula,  data = sdata(object))
# Rename "(Intercept)" column
    if (ncol(myDesign)==1){ colnames(myDesign) <- 'subgroup1'
                            return(myDesign) }
    factors <- svars(object)[are_factor(sdata(object))] # Rename intercept
    factor1 <- factors[1]                               # to factor1level1
    level1  <- levels(sdata(object)[[factor1]])[1]
    colnames(myDesign) %<>% gsub('(Intercept)', level1, ., fixed = TRUE)
# Rename regressors
    for (var in factors) colnames(myDesign) %<>% gsub(var, '', ., fixed = TRUE)
        # Fails for e.g. T2D = YES/NO: a meaningless column "YES" is created
        # For other cases it works wonderfully, so I keep it for now.
        # If it gives too many issues, roll back to doing the dropping only
        # for "subgroup" levels:
        #colnames(myDesign) %<>% gsub('subgroup', '', ., fixed=TRUE)
# Validify names
    colnames(myDesign) %<>% gsub(':', '..', ., fixed = TRUE)
    colnames(myDesign) %<>% make.names()
# Return
    return(myDesign)
}


#=============================================================================
#
#               subgroup_matrix
#                   split_subgroup_levels
#                       split_subgroup_values
#                           split_values
#
#=============================================================================

split_values <- function(x){
    sep <- guess_sep(x)
    dt <- data.table::data.table(x = x)
    dt[, data.table::tstrsplit(x, sep) ]
}

split_subgroup_values <- function(object){
    subgroupvalues <- subgroup_values(object)
    cbind(subgroup = subgroupvalues, split_values(subgroupvalues))
}

split_subgroup_levels <- function(object){
    subgrouplevels <- subgroup_levels(object)
    cbind(subgroup = subgrouplevels, split_values(subgrouplevels))
}


#' @rdname subgroup_matrix
#' @export
subgroup_array <- function(object){
    x <- subgroup_levels(object)
    sep <- guess_sep(object)
    #x %<>% sort()
    dt <- data.table(subgroup = factor(x, x))
    components <- dt[, tstrsplit(subgroup, sep, fixed=TRUE)]
    for (i in seq_len(ncol(components)))   components[[i]] %<>%
                                    factor(., levels=unique(.))
    dt %<>% cbind(components)
    data.table::setorderv(dt, rev(names(components)))
    levels  <- dt[, -1] %>% lapply(unique)
    #levels[1:2] %<>% rev()
    nlevels <- levels %>% vapply(length, integer(1))
    array(dt$subgroup, dim = nlevels, dimnames = levels)
}


#' Get subgroup matrix
#'
#' Arrange (subgroup)levels in matrix
#'
#' @param object SummarizedExperiment
#' @return matrix
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' subgroup_matrix(object)
#' @export
subgroup_matrix <- function(object){
    subgroup_array <- subgroup_array(object)
    if (length(dim(subgroup_array))==1)  return(matrix(subgroup_array,
        byrow=TRUE, nrow=1, dimnames=list(NULL, subgroup_array)))
    otherdims <- names(dim(subgroup_array)) %>% setdiff('V1')
    ncol1   <- Reduce('*', dim(subgroup_array)[otherdims])
    colnames1 <- dimnames(subgroup_array)[otherdims] %>%
                expand.grid()                        %>%
                apply(1, paste0, collapse='.')
    subgroupmat <- matrix(subgroup_array,
                        nrow = nrow(subgroup_array), ncol = ncol1,
                        dimnames=list(rownames(subgroup_array), colnames1))
    subgroupmat %>% extract(nrow(.):1, )
    #dt <- split_subgroup_levels(object)
    #subgroupmat <- as.matrix(data.table::dcast(
    #    dt, V1 ~ V2, value.var = 'subgroup'), rownames = 'V1')
    #subgroupmat %>% extract(rev(order(rownames(.))), order(colnames(.)))
}





#=============================================================================
#
#               contrast_subgroup_cols
#               contrast_subgroup_rows
#               col_contrast_defs
#               row_contrast_defs
#
#==============================================================================


#' Row/Col contrasts
#' @param object SummarizedExperiment
#' @return  matrix
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' subgroup_matrix(object)
#' contrast_subgroup_cols(object)
#' contrast_subgroup_rows(object)
#' @export
contrast_subgroup_cols <- function(object){
    subgroupmat <- subgroup_matrix(object)
    if (is_scalar(subgroupmat))  return(subgroupmat)
    if (ncol(subgroupmat)==1) return(matrix(, ncol=0, nrow=nrow(subgroupmat)))
    colcontrasts <- matrix(  sprintf('%s - %s',
                                    subgroupmat[, -1],
                                    subgroupmat[, -ncol(subgroupmat)]),
                            nrow = nrow(subgroupmat),
                            ncol = ncol(subgroupmat)-1)
    rownames(colcontrasts) <- rownames(subgroupmat)
    colnames(colcontrasts) <- sprintf('%s - %s',
                            colnames(subgroupmat)[-1],
                            colnames(subgroupmat)[-ncol(subgroupmat)])
    colcontrasts
}


#' @rdname contrast_subgroup_cols
#' @export
contrast_subgroup_rows <- function(object){
    subgroupmat <- subgroup_matrix(object)
    if (nrow(subgroupmat)==1) return(matrix(, nrow=0, ncol=ncol(subgroupmat)))
    rowcontrasts <- matrix(  sprintf('%s - %s',
                                    subgroupmat[-nrow(subgroupmat), ],
                                    subgroupmat[-1, ]),
                            nrow = nrow(subgroupmat)-1,
                            ncol = ncol(subgroupmat))
    colnames(rowcontrasts) <- colnames(subgroupmat)
    rownames(rowcontrasts) <- sprintf('%s - %s',
                            rownames(subgroupmat)[-nrow(subgroupmat)],
                            rownames(subgroupmat)[-1])
    rowcontrasts
}

contrast_subgroups <- function(object)  list(contrast_subgroup_cols(object),
                                            contrast_subgroup_rows(object))

#==============================================================================
#
#                            contrast
#
#==============================================================================

contrvec2mat  <- function(contrastdefs)  matrix(
                   contrastdefs, nrow=1, dimnames=list("", contrastdefs))

contrmat2list <- function(contrastdefs)  list(colcontrasts = contrastdefs)


#' Contrast
#'
#' Fit linear model and analyze contrasts
#'
#' Limma results can be easily accessed with limma(object).
#' This returns a list with components:
#' \itemize{
#'    \item {effect} matrix (ngene x ncontrast): effect sizes
#'    \item {rank}   matrix (ngene x ncontrast): effect ranks
#'    \item {t}      matrix (ngene x ncontrast): t    values (moderated t test)
#'    \item {se}     matrix (ngene x ncontrast): se   values (moderated t test)
#'    \item {p}      matrix (ngene x ncontrast): p    values (moderated t test)
#'    \item {fdr}    matrix (ngene x ncontrast): fdr  values (moderated t test)
#'    \item {bonf}   matrix (ngene x ncontrast): bonf values (moderated t test)
#'    \item {F}      vector (ngene)            : F    values (moderated F test)
#'    \item {F.p}    vector (ngene)            : p    values (moderated F test)
#' }
#' @param object       SummarizedExperiment
#' @param contrastdefs contrastdef vector / matrix / list:
#' \itemize{
#' \item{c("t1-t0", "t2-t1", "t3-t2")}
#' \item{matrix(c("WT.t1-WT.t0", "WT.t2-WT.t1", "WT.t3-WT.t2"), \cr
#'      c("KD.t1-KD.t0", "KD.t2-KD.t1", "KD.t3-KD.t2"), nrow=2, byrow=TRUE)}
#' \item{list(matrix(c("WT.t1-WT.t0", "WT.t2-WT.t1", "WT.t3-WT.t2"), \cr
#'      c("KD.t1-KD.t0", "KD.t2-KD.t1", "KD.t3-KD.t2"), nrow=2, byrow=TRUE), \cr
#'      matrix(c("KD.t0-WT.t0", "KD.t1-WT.t1", "KD.t2-WT.t2", "KD.t3-WT.t3"),\cr
#'      nrow=1, byrow=TRUE))}}
#' @param formula      formula to create design matrix (using svars)
#' @param block     block svar (or NULL)
#' @param plot         TRUE/FALSE
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file, plot=FALSE)
#' lmfit(object)
#'
#' file <- download_data('billing19.proteingroups.txt')
#' select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
#' select %<>% paste0('_STD')
#' object <- read_proteingroups(file, select_subgroups = select, plot = FALSE)
#' object %<>% lmfit()
#'
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_rnaseq_counts(file, plot=FALSE)
#' object$subgroup %<>% factor(sort(unique(.))[c(2:length(.), 1)])
#' object %<>% lmfit()
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' object %<>% lmfit()
#' @export
lmfit <- function(object, contrastdefs = NULL,
    formula   = NULL, block = NULL, verbose = TRUE, plot =  TRUE
){
# Initialize
    assert_is_all_of(object, 'SummarizedExperiment')
    if (is.null(contrastdefs)) contrastdefs <- contrast_subgroups(object)
    if (is.character(contrastdefs)) contrastdefs %<>% contrvec2mat()
    if (is.matrix(contrastdefs))    contrastdefs %<>% contrmat2list()
    design <- create_design(object, formula=formula)
    design(object)    <- design
    contrastdefs(object) <- contrastdefs
# Prepare block
    if (verbose)  cmessage('\t\tLmfit/contrast: %s%s',
           if(is.null(block))            '' else paste0(' | ', block),
           if(is.null(weights(object)))  '' else paste0(' (weights)'))
    if (!is.null(block)){
        assert_is_subset(block, svars(object))
        blockvar <- block
        block <- sdata(object)[[block]]
        if (is.null(metadata(object)$blockcor)){
            if (verbose)  cmessage('\t\t\t\tdupcor `%s`', blockvar)
            metadata(object)$blockcor <- duplicateCorrelation(
                exprs(object), design=design, block=block
            )$consensus.correlation }}
# Fit lm and compute contrasts
    fit <- suppressWarnings(lmFit(object = exprs(object[, rownames(design)]),
        design = design, block = block, correlation = metadata(object)$blockcor,
        weights = weights(object)))
    object %<>% add_contrast_results(fit)
# Plot/Return
    if (plot)  plot_contrastogram(object)
    cmessage_df('\t\t\t%s', extract_limma_summary(object))
    return(object)
}

#' @rdname lmfit
#' @export
add_limma <- function(...){
    .Deprecated('lmfit')
    lmfit(...)
}


vectorize_contrastdefs <- function(contrastdefs){
    unname(unlist(lapply(contrastdefs, function(x) na.exclude(c(t(x))))))
}


add_contrast_results <- function(object, fit){
# compute contrasts
    contrastmat <- makeContrasts(
                    contrasts = vectorize_contrastdefs(contrastdefs(object)),
                    levels    = design(object))
    fit %<>% contrasts.fit(contrasts = contrastmat)
    limma_quantities <- if (all(fit$df.residual==0)){ c('effect', 'rank')
                        } else { c('effect','rank','t','se','p','fdr','bonf') }
    limma(object) <- array( dim=c(nrow(fit),ncol(fit),length(limma_quantities)),
                            dimnames = list(feature  = rownames(fit),
                                            contrast = colnames(fit),
                                            quantity = limma_quantities))
    limma(object)[,,'effect'] <- fit$coefficients
    limma(object)[,,'rank'  ] <- apply(-abs(fit$coefficients), 2, rank)
# perform moderated t test
    if (!all(fit$df.residual==0)){
        fit %<>% eBayes()
        pp <- fit$p.value
        limma(object)[,,'t' ] <- fit$t
        limma(object)[,,'se'] <- sqrt(fit$s2.post) * fit$stdev.unscaled
        limma(object)[,,'p' ] <- pp
        limma(object)[,,'rank'] <- apply(pp, 2, rank)
        limma(object)[,,'fdr' ] <- apply(pp, 2, p.adjust, 'fdr')
        limma(object)[,,'bonf'] <- apply(pp, 2, p.adjust, 'bonferroni')
        fdata(object)$F   <- fit$F
        fdata(object)$F.p <- fit$F.p
    }
    return(object)
}


#==============================================================================
#
#                    limma & limma<-
#                    extract_limma_dt
#
#==============================================================================

#' @title Get/set design matrix
#' @param object SummarizedExperiment
#' @param value list
#' @return design (get) or SummarizedExperiment (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(file, invert_subgroups=inv, plot=FALSE)
#' design(object)
#' @export
setGeneric("design", function(object)   standardGeneric("design") )


#' @rdname design
#' @export
setGeneric("design<-",
function(object, value)  standardGeneric("design<-") )


#' @rdname design
setMethod("design", signature("SummarizedExperiment"),
function(object) metadata(object)$design)


#' @rdname design
setReplaceMethod("design", signature("SummarizedExperiment", "matrix"),
function(object, value){
    metadata(object)$design <- value
    object  })


#' @title Get/set contrastdefs
#' @param object SummarizedExperiment
#' @param value list
#' @return contrastdefs (get) or SummarizedExperiment (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(file, invert_subgroups=inv, plot=FALSE)
#' contrastdefs(object)
#' @export
setGeneric("contrastdefs", function(object)   standardGeneric("contrastdefs") )


#' @rdname contrastdefs
setMethod("contrastdefs", signature("SummarizedExperiment"),
function(object) metadata(object)$contrastdefs)


#' @rdname contrastdefs
#' @export
setGeneric("contrastdefs<-",
function(object, value)  standardGeneric("contrastdefs<-") )


#' @rdname contrastdefs
setReplaceMethod("contrastdefs", signature("SummarizedExperiment", "list"),
function(object, value){
    metadata(object)$contrastdefs <- value
    object  })


#' @title Get/set limma results
#' @description Get/Set limma results
#' @param object SummarizedExperiment
#' @param value list
#' @return limma results (get) or updated object (set)
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(file, invert_subgroups=inv, plot=FALSE)
#' dim(limma(object))
#' dim(limma(object[1:5, ]))
#' @export
setGeneric("limma", function(object)   standardGeneric("limma") )


#' @rdname limma
setMethod("limma", signature("SummarizedExperiment"),
function(object){
    limma_array <- metadata(object)$limma
    if (is.null(limma_array)) NULL else limma_array[
                                            fnames(object), , , drop=FALSE] })

#' @rdname limma
#' @export
setGeneric("limma<-", function(object, value)  standardGeneric("limma<-") )


#' @rdname limma
setReplaceMethod("limma", signature("SummarizedExperiment", "array"),
function(object, value){
    metadata(object)$limma <- value
    object  })


#' @rdname limma
setReplaceMethod("limma", signature("SummarizedExperiment", "NULL"),
function(object, value) object)



#' Extract limma quantity
#' @param object SummarizedExperiment
#' @param quantity 'effect', 'p', 'fdr', 'bonf'
#' @return melted data.table
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(file, invert_subgroups=inv, plot=FALSE)
#' extract_limma_quantity(object)
#' @noRd
extract_limma_quantity <- function(object, quantity='p'){
    fvars0 <- c('feature_id','feature_name','imputed')
    fvars0 %<>% intersect(fvars(object))
    dt <- data.table(fdata(object)[, fvars0, drop=FALSE])
    dt %<>% cbind(adrop(limma(object)[, , quantity, drop=FALSE], drop=3))
    data.table::melt.data.table(
        dt, id.vars = fvars0, variable.name = 'contrast', value.name = quantity)
}


merge_limma_quantities <- function(x, y){
    names0 <- c('feature_id','feature_name','imputed', 'contrast')
    names0 %<>% intersect(names(x)) %>% intersect(names(y))
    merge(x, y, by = names0)
}


#' Extract limma datatable
#' @param object SummarizedExperiment
#' @return melted data.table
#' @examples
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_proteingroups(file)
#' extract_limma_dt(object)
#'
#' file <- download_data('billing16.proteingroups.txt')
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(file, invert_subgroups=inv, plot=FALSE)
#' extract_limma_dt(object)
#' @noRd
extract_limma_dt <- function(object){
    Reduce( merge_limma_quantities,
            mapply( extract_limma_quantity,
                    quantity = c('effect', 'p', 'fdr', 'bonf'),
                    MoreArgs = list(object=object),
                    SIMPLIFY = FALSE ))
}

#' Extract contrast analysis summary
#' @param object SummarizedExperiment
#' @examples
#' # RNASEQCOUNTS
#'     # ~ 0 + subgroup
#'         file <- download_data('billing19.rnacounts.txt')
#'         object <- read_rnaseq_counts(file, plot=FALSE)
#'         extract_limma_summary(object)
#'
#'     # ~ 0 + subgroups | weights
#'         weights(object) <- NULL
#'         object %<>% lmfit(plot=FALSE)
#'         extract_limma_summary(object)
#'
#' # METABOLON
#'     # ~ 0 + subgroup
#'           file <- download_data('atkin18.metabolon.xlsx')
#'           object <- read_metabolon(file, plot=FALSE)
#'           extract_limma_summary(object)
#'
#'    # ~ 0 + subgroup | block
#'           svars(object) %<>% stri_replace_first_fixed('SUB', 'block')
#'           object %<>% lmfit(plot=FALSE)
#'           extract_limma_summary(object)
#'
#'    # ~ 0 + subgroup + t2d | block
#'           object %<>% lmfit(formula=~0+subgroup+T2D, plot=FALSE)
#'           extract_limma_summary(object)
#' @export
extract_limma_summary <- function(object){
    effect <- fdr <- NULL
    extract_limma_dt(object)[,
        .(ndown = sum(effect<0 & fdr<0.05, na.rm=TRUE),
          nup   = sum(effect>0 & fdr<0.05, na.rm=TRUE)),
        by='contrast']
}

#=============================================================================
#
#             plot_contrastogram
#                 compute_connections
#
#=============================================================================


true_names <- function(x) names(x[x])

compute_connections <- function(
    object, colors = make_colors(subgroup_levels(object), guess_sep(object))
){
# subgroup matrix, difference contrasts, limma
    pvalues <- limma(object)[, , 'p',      drop=FALSE]
    effects <- limma(object)[, , 'effect', drop=FALSE]
    nsignif <- apply(pvalues < 0.05, 2, sum, na.rm=TRUE)
                #colSums( pvalues < 0.05, na.rm=TRUE)  # BREAKS ON SINGLE CONTR!
    nup     <- apply(pvalues < 0.05 & effects>0, 2, sum, na.rm=TRUE)
                #colSums((pvalues < 0.05) & (effects > 0), na.rm=TRUE)
    ndown   <- apply(pvalues < 0.05 & effects<0, 2, sum, na.rm=TRUE)
                #colSums((pvalues < 0.05) & (effects < 0), na.rm=TRUE)
# Create diagram
    sep <- guess_sep(object)
    subgroupmatrix <- subgroup_matrix(object)
    subgrouplevels <- c(t(subgroupmatrix))
    arrowsizes <- arrowcolors <- matrix(0,
        nrow = length(subgrouplevels), ncol = length(subgrouplevels),
        dimnames = list(subgrouplevels, subgrouplevels))
    arrowlabels <- matrix("0", nrow = nrow(arrowsizes), ncol = ncol(arrowsizes),
                        dimnames = dimnames(arrowsizes))
# Add contrast numbers
    design    <- design(object)
    colcontrasts <- contrastdefs(object)[[1]]
    rowcontrasts <- contrastdefs(object)[[2]]
    contrastdefs <-  c( c(t(colcontrasts)), c(t(rowcontrasts)))
    contrastmat  <- makeContrasts(contrasts = contrastdefs, levels = design)
    for (contrastname in colnames(contrastmat)){
        contrastvector <- contrastmat[, contrastname]
        to   <- true_names(contrastvector>0)
        from <- if (any(contrastvector<0)) true_names(contrastvector<0) else to
        ns <- nsignif[[contrastname]]
        nu <- nup[[contrastname]]
        nd <- ndown[[contrastname]]
        arrowsizes[ to, from] <- nu#ns
        arrowsizes[ from, to] <- nd#ns
        arrowcolors[to, from] <- colors[[to]]
        arrowcolors[from, to] <- colors[[from]]
        arrowlabels[to, from] <- if (nu>0) nu else 0
                            #paste0(nu,  " %up% phantom(.)") else "phantom(.)"
        arrowlabels[from, to] <- if (nd>0) nd else 0
                            #paste0(nd," %down% phantom(.)") else "phantom(.)"
    }
# Return
    #arrowlabels[arrowcolors==0] <- "0"
    list(arrowsizes = arrowsizes,
        arrowcolors = arrowcolors,
        arrowlabels = arrowlabels)
}


#' Plot contrastogram
#' @param object SummarizedExperiment
#' @param colors named color vector (names = subgroups)
#' @param curve  arrow curvature
#' @return list returned by \code{\link[diagram]{plotmat}}
#' @examples
#' # subgroup matrix
#'    file <- download_data('halama18.metabolon.xlsx')
#'    object <- read_metabolon(file, plot=FALSE)
#'    plot_contrastogram(object)
#' # subgroup vector
#'     require(magrittr)
#'     file <-  download_data('billing19.proteingroups.txt')
#'     select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
#'     select %<>% paste0('_STD')
#'     object <- read_proteingroups(file, select_subgroups = select, plot=FALSE)
#'     object %<>% lmfit(plot=FALSE)
#'     plot_contrastogram(object, curve=0.8)
#' # subgroup vector
#'     file <-  download_data('fukuda20.proteingroups.txt')
#'     object <- read_proteingroups(file, plot=FALSE)
#'     plot_contrastogram(object)
#' # Ratios: self-contrasts
#'    file <- download_data('billing16.proteingroups.txt')
#'    invert <- c('EM_E', 'BM_E', 'BM_EM')
#'    object <- read_proteingroups(file, invert_subgroups=invert, plot=FALSE)
#'    plot_contrastogram(object)
#' @export
plot_contrastogram <- function(
    object,
    colors = make_colors(subgroup_levels(object), guess_sep(object)),
    curve  = 0.1
){
# Initialize
    V2 <- N <- NULL
# Prepare
    contrastogram_matrices <- compute_connections(object, colors = colors)
    arrowsizes  <- contrastogram_matrices$arrowsizes
    arrowcolors <- contrastogram_matrices$arrowcolors
    arrowlabels <- contrastogram_matrices$arrowlabels
    widths <- scales::rescale(arrowsizes, c(0.01,30))
# Plot
    dt <- split_subgroup_levels(object)
    nrow <- dt[, data.table::uniqueN(V2)]
    nperrow <- dt[, .N, by = 'V1'][, N]
    if (all(nperrow==1)) nperrow %<>% length()
    #dir.create('~/autonomicscache/contrastogram')
    #pdf('~/autonomicscache/contrastogram/directed_contrastogram.pdf',
    #width = 9, height = 9)
    arrowlabels %<>% as.data.frame()
    diagram::plotmat(A          = arrowlabels,
                    pos         = nperrow,
                    curve       = curve,
                    name        = rownames(arrowsizes),
                    relsize     = 1,
                    box.size    = 0.05,
                    box.col     = colors[rownames(arrowsizes)],
                    box.type    = 'square',
                    box.prop    = 0.8,
                    arr.lwd     = widths,
                    shadow.size = 0, # sqrt(arrowsizes)
                    arr.lcol    = arrowcolors,
                    arr.col     = arrowcolors,
                    arr.type    = 'triangle')
    #, arr.lcol = log2(1+diagram_matrix))
    #dev.off()
}



#==============================================================================
#
#          plot_volcano
#              make_volcano_dt
#                  top_down/top_up
#                      nmax/nmin
#
#==============================================================================

invwhich <- function(indices, totlength) is.element(seq_len(totlength), indices)

#' Return nth max (min) value in vector
#'
#' Orders a vector and returns n'th ordered value.
#' When vector length is smaller than n, returns last value.
#'
#' @param x numeric vector
#' @param n integer
#' @return value
#' @examples
#' nmax(c(1,2,3,4,5), 2)
#' nmin(c(1,2,3,4,5), 2)
#' @noRd
nmax <- function(x, n) sort(x, decreasing = TRUE) %>% extract(min(length(.), n))
nmin <- function(x, n) sort(x) %>% extract(min(length(.), n))

top_down <- function(effect, fdr, mlp, ntop){
    fdr_ok   <- fdr  < 0.05
    coef_ok  <- effect < -1
    coef_top <- if (any(fdr_ok)){  effect < nmin(effect[fdr_ok], ntop+1)
                } else {           rep(FALSE, length(effect))            }
    mlp_top  <- if (any(coef_ok)){ mlp  > nmax(mlp[coef_ok], ntop+1)
                } else {           rep(FALSE, length(effect))            }
    fdr_ok & coef_ok & (coef_top | mlp_top)
}

#' @examples
#' file <- download_data("billing16.proteingroups.txt")
#' invert_subgroups <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(file, invert_subgroups=invert_subgroups)
#' effect <-      limma(object)[,1,'effect']
#' fdr    <-      limma(object)[,1,'fdr']
#' mlp    <- -log(limma(object)[,1,'p'])
#' ntop   <- 3
#' table(top_up(effect, fdr, mlp, ntop))
#' table(top_down(effect, fdr, mlp, ntop))
#' @noRd
top_up <- function(effect, fdr, mlp, ntop){
    fdr_ok   <- fdr  < 0.05
    coef_ok  <- effect >  1
    coef_top <- if(any(fdr_ok)){  effect > nmax(effect[fdr_ok], ntop+1)
                } else {          rep(FALSE, length(effect)) }
    mlp_top  <- if (any(coef_ok)){ mlp > nmax(mlp[coef_ok], ntop+1)
                } else {           rep(FALSE, length(effect)) }
    fdr_ok & coef_ok & (coef_top | mlp_top)
}


melt_contrastdefs <- function(contrastdefmat){
    facetrow <- NULL
    contrastdefdt <- data.table(contrastdefmat, facetrow = "")
    if (!is.null(rownames(contrastdefmat))) contrastdefdt[,
                                        facetrow := rownames(contrastdefmat)]
    data.table::melt.data.table(
        contrastdefdt,
        id.vars       = 'facetrow',
        variable.name = 'facetcol',
        value.name    = 'contrast', by = 'contrast')
}


#' Create volcano datatable
#' @param object          SummarizedExperiment
#' @param contrastdefmat  contrastdef matrix
#' @param ntop            no of top features to be annotated
#' @return data.table
#' @examples
#' file <- download_data("billing16.proteingroups.txt")
#' inv <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(file, invert_subgroups=inv,
#'             contrastdefs = c('E_EM', 'E_BM', 'EM_BM'))
#' make_volcano_dt(object)
#' @export
make_volcano_dt <- function(
    object, contrastdefmat = contrastdefs(object)[[1]], ntop = 3
){
    effect <- p <- mlp <- topdown <- topup <- significance <- fdr <- NULL
    dt <- extract_limma_dt(object)
    dt %<>% merge(melt_contrastdefs(contrastdefmat), by = 'contrast')
    dt %<>% extract(!is.na(effect) & !is.na(p))
    dt[, mlp  := -log10(p)]

    # Prepare volcano datatable
    # Note: Using effect <= 0 (rather than effect <0) is required.
    # Otherwise the (very few) features with effect=0 will have no effect for
    # 'significance'
    by <- intersect(c('contrast', 'imputed'), names(dt))
    dt[,topdown := top_down(effect, fdr, mlp, ntop), by=by]
    dt[,topup   := top_up(  effect, fdr, mlp, ntop), by=by]
    dt[effect<=0,            significance := 'down']
    dt[effect> 0,            significance :=   'up']
    dt[effect<=0 & fdr<0.05, significance := 'down: fdr<0.05']
    dt[effect> 0 & fdr<0.05, significance :=   'up: fdr<0.05']
    dt[topdown==TRUE,        significance := 'down: top']
    dt[topup  ==TRUE,        significance :=   'up: top']
    dt[,significance := factor(significance, c(
        'down: top','down: fdr<0.05','down','up','up: fdr<0.05','up: top'))]
    dt[]
}

#' Plot volcano
#'
#' @param object           SummarizedExperiment
#' @param contrastdefmat   contrast layout matrix
#' @param label            fvar for labeling top features
#' @param ntop             number: n top features to be annotated
#' @return ggplot object
#' @examples
#' # proteingroup group ratios
#'     require(magrittr)
#'     file <- download_data("billing16.proteingroups.txt")
#'     inv <- c('EM_E', 'BM_E', 'BM_EM')
#'     object <- read_proteingroups(file, invert_subgroups=inv, plot=FALSE)
#'     plot_volcano(object)
#'     contrasts <- subgroup_matrix(object)
#'     object %<>% lmfit(contrasts = contrasts, plot = FALSE)
#'     plot_volcano(object)
#'
#' # proteingroup LFQ intensities
#'     file <- download_data('fukuda20.proteingroups.txt')
#'     object <- read_proteingroups(file, plot=FALSE)
#'     plot_volcano(object)
#'
#' # metabolon intensities: complex design
#'     file <- download_data('halama18.metabolon.xlsx')
#'     object <- read_metabolon(file, plot=FALSE)
#'     plot_volcano(object)
#'
#' # proteingroup internalstandard ratios
#'     file <-  download_data('billing19.proteingroups.txt')
#'     select <-  c('E00','E01', 'E02','E05','E15','E30', 'M00')
#'     select %<>% paste0('_STD')
#'     object <- read_proteingroups(file, select_subgroups = select, plot=FALSE)
#'     plot_volcano(object)
#' @export
plot_volcano <- function(
    object, contrastdefmat = contrastdefs(object)[[1]],
    label = feature_name, ntop = 3
){
# Assert
    assert_is_all_of(object, "SummarizedExperiment")
    assert_is_matrix(contrastdefmat)
    topup <- topdown <- effect <- mlp <- facetrow <- facetcol <- NULL
    label <- enquo(label)
# Prepare
    plotdt <- make_volcano_dt(object, contrastdefmat, ntop = ntop)
    txtdt  <- copy(plotdt)[topup==TRUE | topdown==TRUE]
    colorvalues <-c(hcl(h=  0, l=c(20, 70, 100), c=100),
                    hcl(h=120, l=c(100, 70, 20), c=100))
    names(colorvalues) <- levels(plotdt$significance)
# Plot
    imputed <- NULL # fallback when plotdt misses "imputed"
    significance <- NULL
    p <- ggplot(plotdt) +
        facet_grid(
            rows = if (all(stri_isempty(plotdt$facetrow))) {NULL} else {
                vars(facetrow)},
            cols = if (all(stri_isempty(plotdt$facetcol))) {NULL} else {
                vars(facetcol)},
            scales = 'fixed') +
        geom_point( aes(x=effect, y=mlp, color = significance, shape=imputed),
                    na.rm = TRUE)
    if (!quo_is_null(label)){
        p <- p + geom_text_repel(
                        data = txtdt,
                        aes(x=effect, y=mlp, label=!!label, color=significance),
                        #hjust = 'outward',
                        na.rm = TRUE,
                        show.legend = FALSE)}#,
                        #direction = 'x'
    p + theme_bw() +
        scale_color_manual(values = colorvalues, name = NULL) +
        xlab(expression(log[2](FC))) +
        ylab(expression(-log[10](p))) +
        ggtitle('volcano')#+
        #guides(color=FALSE)
}
