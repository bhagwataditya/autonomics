#==============================================================================
#
#                is_valid_contrast
#                select_valid_contrast
#                validify_contrasts
#
#==============================================================================


#' Is a valid contrast?
#' @param contrast contrast
#' @param design   design
#' @return logical
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' design <- create_design(object)
#' contrast <- default_contrasts(object)[1]
#' is_valid_contrast(contrast, design)
#' @noRd
is_valid_contrast <- function(contrast, design){
  contrast %<>% as.character()
  contrast %<>% stri_replace_all_fixed(' ', '')
  contrast %<>% stri_replace_all_regex('(/[0-9]+)', '') %>%
                stri_replace_all_fixed('(', '')         %>%
                stri_replace_all_fixed(')', '')
  terms <- strsplit(contrast, '[-+ ]+') %>% unlist() %>% unname()
  all(terms %in% colnames(design))
}


#' Select valid contrasts
#' @param  contrasts vector with contrast definitions
#' @param  design    design matrix
#' @param  verbose   whether or not to report
#' @return subset of contrasts
#' @noRd
select_valid_contrasts <- function(contrasts, design, verbose = TRUE){
  selector <- contrasts %>% vapply(is_valid_contrast, logical(1), design)
  if (verbose){
    cmessage('\t\tKeep %d valid contrasts (out of %d)',
            sum(selector), length(selector))
  }
  contrasts[selector]
}


#' Validy subgroups and contrasts
#' @param design    design matrix
#' @param contrasts contrasts vector
#' @return validified contrasts
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object  <- read_proteingroups(file)
#' design  <- create_design(object)
#' contrasts <- default_contrasts(object)
#' validify_contrasts(contrasts, design)
#' @noRd
validify_contrasts <- function(contrasts, design){

  # Assert valid inputs
   assert_is_matrix(design)
   assert_is_numeric(design)

  # Validify contrast names
  if (!has_names(contrasts)){
    names(contrasts) <- make.names(contrasts)
  }
  assert_has_no_duplicates(names(contrasts))

  # Select valid contrasts
  contrasts %<>% select_valid_contrasts(design)
  contrasts
}


#==============================================================================
#
#                       create_contrastmat
#
#==============================================================================


#' Create contrast matrix
#' @param object        SummarizedExperiment
#' @param contrasts  vector of contrast definitions
#' @param design     design matrix
#' @return contrast matrix
#' @examples
#' # PROTEINGROUPS
#' file <- download_data('billing16.proteingroups.txt')
#' invert_subgroups <- c('BM_EM', 'EM_E', 'BM_E')
#' object <- read_proteingroups(file, invert_subgroups = invert_subgroups)
#' create_contrastmat(object)
#'
#' # RNACOUNTS
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_counts(file)
#' create_contrastmat(object, c(EM0.8_0 = 'EM.8 - EM00'))
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' create_contrastmat(object)
#' @export
create_contrastmat <- function(
    object, contrasts = default_contrasts(object),
    design = importomics::design(object) # be explicit to disambiguate!
){
    if (!has_names(contrasts)) names(contrasts) <-make.names(contrasts)
    if (length(contrasts) == 0)  return(NULL)
    assert_has_no_duplicates(names(contrasts))
    makeContrasts(contrasts = contrasts, levels = design) %>%
    set_colnames(names(contrasts))
}



#=============================================================================
#
#               matrixify_subgroups
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


#' Arrayify (subgroup) levels
#'
#' Arrange (subgroup) levels in array
#'
#' @param x  subgroup levels (character vector)
#' @return array
#' @examples
#' array_levels(x = c('wt',             'kd'))
#' array_levels(x = c('wt.t0',     'wt.t1',     'kd.t0',     'kd.t1'))
#' array_levels(x = c('wt.t0.uM0', 'wt.t1.uM0', 'kd.t0.uM0', 'kd.t1.uM0',
#'                    'wt.t0.uM5', 'wt.t1.uM5', 'kd.t0.uM5', 'kd.t1.uM5'))
#' @noRd
arrayify_subgroups <- function(x){
    x %<>% sort()
    dt <- data.table(subgroup = x)
    components <- dt[, tstrsplit(subgroup, guess_sep(x), fixed=TRUE)]
    dt %<>% cbind(components)
    data.table::setorderv(dt, rev(names(components)))
    levels  <- dt[, -1] %>% lapply(unique)
    nlevels <- levels %>% vapply(length, integer(1))
    array(dt$subgroup, dim = nlevels, dimnames = levels)
}

#' Matrixify subgroups
#'
#' Arrange (subgroup) levels in matrix
#'
#' @param x  subgroup levels (character vector)
#' @return matrix
#' @examples
#' matrixify_subgroups(x = c('wt', 'kd'))
#' matrixify_subgroups(x = c('wt.t0', 'wt.t1', 'kd.t0', 'kd.t1'))
#' matrixify_subgroups(x = c('wt.t0.uM0', 'wt.t1.uM0', 'kd.t0.uM0', 'kd.t1.uM0',
#'                    'wt.t0.uM5', 'wt.t1.uM5', 'kd.t0.uM5', 'kd.t1.uM5'))
#' file <- download_data('halama18.metabolon.xlsx')
#' matrixify_subgroups(x = subgroup_levels(read_metabolon(file, plot=FALSE)))
#' @noRd
matrixify_subgroups <- function(x){
    subgroup_array <- arrayify_subgroups(x)
    if (length(dim(subgroup_array))==1)  return(matrix(subgroup_array,
                      byrow=TRUE, nrow=1, dimnames=list(NULL, subgroup_array)))
    otherdims <- names(dim(subgroup_array)) %>% setdiff('V1')
    ncol1   <- Reduce('*', dim(subgroup_array)[otherdims])
    colnames1 <- dimnames(subgroup_array)[otherdims] %>%
                expand.grid()                        %>%
                apply(1, paste0, collapse='.')
    subgroup_matrix <- matrix(subgroup_array,
           nrow = nrow(subgroup_array), ncol = ncol1,
           dimnames=list(rownames(subgroup_array), colnames1))
    subgroup_matrix %>% extract(nrow(.):1, )
    #dt <- split_subgroup_levels(object)
    #subgroup_matrix <- as.matrix(data.table::dcast(
    #    dt, V1 ~ V2, value.var = 'subgroup'), rownames = 'V1')
    #subgroup_matrix %>% extract(rev(order(rownames(.))), order(colnames(.)))
}





#=============================================================================
#
#               contrast_columns
#               contrast_rows
#
#==============================================================================


#' Contrast columns
#'
#' Contrast columns within a row
#'
#' Contrast timepoints within conc
#'
#' @param subgroup_matrix subgroup matrix (nconc x ntime)
#' @return      contrast matrix: nconc x (ntime-1)
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' subgroup_matrix <- matrixify_subgroups(subgroup_levels(object))
#' contrast_columns(subgroup_matrix)
#' @noRd
contrast_columns <- function(subgroup_matrix, symbol = ' - '){

    contrastmat <- matrix(  sprintf('%s%s%s',
                                    subgroup_matrix[, -1],
                                    symbol,
                                    subgroup_matrix[, -ncol(subgroup_matrix)]),
                            nrow = nrow(subgroup_matrix),
                            ncol = ncol(subgroup_matrix)-1)

    rownames(contrastmat) <- rownames(subgroup_matrix)
    colnames(contrastmat) <- sprintf('%s - %s',
                            colnames(subgroup_matrix)[-1],
                            colnames(subgroup_matrix)[-ncol(subgroup_matrix)])
    contrastmat
}


#' Contrast rows
#'
#' Contrast rows within a column
#'
#' Contrast concentrations within timepoint
#'
#' @param subgroup_matrix subgroup matrix: nconc x ntime
#' @return      contrast matrix: (nconc-1) x ntime
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' subgroup_matrix <- matrixify_subgroups(subgroup_levels(object))
#' contrast_rows(subgroup_matrix)
#' @noRd
contrast_rows <- function(subgroup_matrix, symbol = ' - '){

    contrastmat <- matrix(  sprintf('%s%s%s',
                                  subgroup_matrix[-nrow(subgroup_matrix), ],
                                  symbol,
                                  subgroup_matrix[-1, ]),
                            nrow = nrow(subgroup_matrix)-1,
                            ncol = ncol(subgroup_matrix))

    colnames(contrastmat) <- colnames(subgroup_matrix)
    rownames(contrastmat) <- sprintf('%s - %s',
                            rownames(subgroup_matrix)[-nrow(subgroup_matrix)],
                            rownames(subgroup_matrix)[-1])
    contrastmat
}


#=============================================================================
#
#               aggregate_col_contrasts
#               aggregate_row_contrasts
#                   aggregate_contrasts
#
#==============================================================================


#' Aggregate contrasts
#'
#' Aggregate row contrasts across columns (or column contrasts across rows)
#'
#' @param contrastmat contrast matrix
#' @param dim 1 (aggregate across rows) or 2 (aggregate across columns)
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' (x <- matrixify_subgroups(subgroup_levels(object)))
#' aggregate_contrasts(contrast_rows(x), 1)    # some conc across t
#' aggregate_contrasts(contrast_rows(x), 2)    # concentrations at t
#' aggregate_contrasts(contrast_columns(x), 2) # some time across conc
#' aggregate_contrasts(contrast_columns(x), 1) # times at conc
#' @noRd
aggregate_contrasts <- function(contrastmat, dim){
    apply(contrastmat,
          dim,
          function(x){
              paste0(sprintf('(%s)/%d', x, length(x)), collapse = ' + ')})
}

aggregate_col_contrasts <- function(contrastmat){
    aggregate_contrasts(contrastmat, 2)
}

aggregate_row_contrasts <- function(conc_contrastmat){
    aggregate_contrasts(conc_contrastmat, 1)
}


#=============================================================================
#
#               create_diff_contrasts
#
#==============================================================================


#' Create diff contrasts
#'
#' Create vector with difference contrasts
#'
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' create_diff_contrasts(object)
#' @noRd
create_diff_contrasts <- function(object){
# Single subgroup
    subgroups <- subgroup_levels(object)
    if (assertive::is_scalar(subgroups)) return(
        structure(subgroups, names=subgroups))
# Row contrasts
    subgroup_matrix <- matrixify_subgroups(subgroup_levels(object))
    if (nrow(subgroup_matrix)>1){
        row_contrmat   <- contrast_rows(subgroup_matrix)
        row_contrnames <- contrast_rows(subgroup_matrix, '__')
        contrasts  <- structure(c(row_contrmat), names = c(row_contrnames))
    }
  # Column contrasts
    if (ncol(subgroup_matrix)>1){
        col_contrmat   <- contrast_columns(subgroup_matrix)
        col_contrnames <- contrast_columns(subgroup_matrix, '__')
        contrasts %<>% c(structure(c(col_contrmat), names = c(col_contrnames)))
    }
# Aggregated contrasts
    # aggregated_col_contrasts <- aggregate_col_contrasts(col_contrmat)
    # aggregated_row_contrasts <- aggregate_row_contrasts(row_contrmat)
    # contrasts %<>% c(aggregated_col_contrasts)
    # contrasts %<>% c(aggregated_row_contrasts))
# Return
    contrasts
}


#=============================================================================
#
#               default_contrasts
#
#==============================================================================


validify_contrast_names <- function(x){
   x %>% gsub(' ', '',  ., fixed = TRUE) %>%
         gsub('-', '_', ., fixed = TRUE) %>%
         make.names()
}


#' Default contrasts
#' @param object SummarizedExperiment
#' @return named character vector: contrast definitions
#' @examples
#' # STEM CELL COMPARISON
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' invert_subgroups <- c('EM_BM', 'EM_E', 'E_BM')
#' object <- read_proteingroups(file, invert_subgroups=invert_subgroups)
#' default_contrasts(object)
#'
#' # GLUTAMINASE
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' default_contrasts(object)
#' @export
default_contrasts <- function(object){

    # Ratios themselves for ratio data
    if (contains_ratios(object)){
        message('\tGenerate contrasts from ratios')
        return(
            subgroup_levels(object) %>% set_names(validify_contrast_names(.)))

    # Difference contrasts for abundance data
    } else {
        message('\tGenerate difference contrasts')
        contrasts <- create_diff_contrasts(object)
    }

    # Return
    contrasts
}


#==============================================================================
#
#                            add_limma
#
#==============================================================================


#' Add limma
#'
#' Run limma and add results
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
#' @param object        SummarizedExperiment
#' @param contrasts  contrast vector, preferably named
#'                   (automatically generated names are not always intuitive)
#' @param formula   formula to create design matrix (using svars)
#' @return Updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('atkin18.somascan.adat')
#' object <- read_somascan(file)
#' add_limma(object)
#'
#' file <- download_data('billing19.proteingroups.txt')
#' rm_subgroups<-c('BLANK_BM00','BM00_BM00','EM01_EM00','EM05_EM02','EM30_EM15')
#' object <- read_proteingroups(file, rm_subgroups = rm_subgroups)
#' contrasts <- c(EM01_EM00 = 'EM01_STD - EM00_STD')
#' object %<>% add_limma(contrasts)
#' sum(limma(object)[,'EM01_EM00','bonf'] < 0.05, na.rm = TRUE)
#'
#' file <- download_data('billing19.rnacounts.txt')
#' object <- read_counts(file)
#' contrasts <- c(EM01_EM00 = 'EM01 - EM00')
#' object %<>% add_limma(contrasts)
#' sum(limma(object)[,'EM01_EM00','bonf'] < 0.05, na.rm = TRUE)
#' @export
add_limma <- function(object, contrasts = default_contrasts(object),
    formula = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup){
# Assert
    design <- create_design(object, formula=!!enquo(formula))
    assert_is_matrix(design)
    assert_is_numeric(design)
    assert_is_identical_to_true(unname(ncol(object)) == nrow(design))
    if (is.null(contrasts))    return(object)
    contrasts %<>% validify_contrasts(design)
    if (length(contrasts)==0) return(object)
# Set block and correlation if required
    cmessage('\t\tRun limma')
    block <- NULL; correlation <- NULL
    my_sdata <- sdata(object)
    if (has_complete_block_values(object)){
        cmessage("\t\tBlock on svar 'block'")
        block <- my_sdata$block
        correlation  <- duplicateCorrelation(exprs(object), design,
                                     block = block)[['consensus.correlation']]}
# Fit lm and compute contrasts
    contrastmat <- create_contrastmat(object, contrasts, design)
    fit <- suppressWarnings(lmFit(object = exprs(object), design = design,
                                  block = block, correlation = correlation,
                                  weights = weights(object)))
    fit %<>% contrasts.fit(contrasts = contrastmat)
    limma_quantities <- if (all(fit$df.residual==0)){ c('effect', 'rank')
                        } else { c('effect','rank','t','se','p','fdr','bonf')}
    limma(object) <- array( dim=c(nrow(fit),ncol(fit),length(limma_quantities)),
                            dimnames = list(feature  = rownames(fit),
                                            contrast = colnames(fit),
                                            quantity = limma_quantities))
    limma(object)[,,'effect'] <- fit$coefficients
    limma(object)[,,'rank'  ] <- apply(-abs(fit$coefficients), 2, rank)
# Perform moderated t test
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
#                    limma<-
#
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




#=============================================================================
#
#             plot_contrastogram
#                 compute_connections
#                     default_color_values2
#                         default_color_values
#                             default_color_var
#
#=============================================================================


#' default color_var
#' @param object SummarizedExperiment
#' @return default value of color_var
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' rm_subgroups <- c('EM_E', 'BM_E', 'EM_BM')
#' object <- read_proteingroups(file, rm_subgroups = rm_subgroups)
#' default_color_var(object)
#'
#' file <- download_data('billing16.somascan.adat')
#' object <- read_somascan(file)
#' default_color_var(object)
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' default_color_var(object)
#' @export
default_color_var <- function(object){
   if (     'block'    %in% svars(object))  'block'
   else if ('subgroup' %in% svars(object))  'subgroup'
   else                                      NULL
}


#' Default color values
#' @param object     SummarizedExperiment
#' @param color_var  string: svar mapped to color
#' @param show       logical
#' @return default color values vector
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' rm_subgroups <- c('EM_E', 'BM_E', 'EM_BM')
#' object <- read_proteingroups(file, rm_subgroups = rm_subgroups)
#' default_color_values(object, show = TRUE)
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' default_color_values(object, show = TRUE)
#' @export
default_color_values <- function(
   object,
   color_var = default_color_var(object),
   show = FALSE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(color_var, svars(object))
# sep
    sep              <- guess_sep(object, color_var)
    color_var_levels <- slevels(object, color_var)
# Make colors
    color_values <- make_colors(color_var_levels, sep, show)
    color_values[color_var_levels]
}


default_color_values2 <- function(object){
    default_color_values(object)[
        c(t(matrixify_subgroups(subgroup_levels(object))))]
}


compute_connections <- function(
    object, contrasts, subgroup_colors = default_color_values2(object)
){
# subgroup matrix, difference contrasts, limma
    pvalues <- limma(object)[, , 'p']
    effects <- limma(object)[, , 'effect']
    nsignif <- colSums( pvalues < 0.05, na.rm=TRUE)
    nup     <- colSums((pvalues < 0.05) & (effects > 0), na.rm=TRUE)
    ndown   <- colSums((pvalues < 0.05) & (effects < 0), na.rm=TRUE)
# Create diagram
    subgrouplevels <- c(t(matrixify_subgroups(subgroup_levels(object))))
    sizes <- colors <- matrix(0,
        nrow = length(subgrouplevels), ncol = length(subgrouplevels),
        dimnames = list(subgrouplevels, subgrouplevels))
    labels <- matrix("0", nrow = nrow(sizes), ncol = ncol(sizes),
                     dimnames = dimnames(sizes))
# Add column contrasts
    subgroup_matrix <- matrixify_subgroups(subgroup_levels(object))
    for (i in nrow(subgroup_matrix):1){
    for (j in 2:ncol(subgroup_matrix)){
        from <- subgroup_matrix[i, j-1]
        to   <- subgroup_matrix[i, j  ]
        ctr <- paste0(to, '__', from)
        ns <- nsignif[[ctr]]; nu <- nup[[ctr]]; nd <- ndown[[ctr]]
        sizes[ to, from] <- ns
        colors[to, from] <- subgroup_colors[[to]]
        labels[to, from] <- if (nu>0) paste0(nu,  " %up% phantom(.)") else "phantom(.)"
        labels[from, to] <- if (nd>0) paste0(nd," %down% phantom(.)") else "phantom(.)"}}
# Add row contrasts
    for (j in 1:ncol(subgroup_matrix)){
    for (i in nrow(subgroup_matrix):2){
        from <- subgroup_matrix[i,   j]
        to   <- subgroup_matrix[i-1, j]
        ctr  <- paste0(to, '__', from)
        ns <- nsignif[[ctr]]; nu <- nup[[ctr]]; nd <- ndown[[ctr]]
        sizes[ to, from] <- ns
        colors[to, from] <- subgroup_colors[[to]]
        labels[to, from] <- if (nd>0) paste0(nd," %down% phantom(.)") else "phantom(.)"
        labels[from, to] <- if (nu>0) paste0(nu,  " %up% phantom(.)") else "phantom(.)"}}
# Return
    #labels[colors==0] <- "0"
    list(sizes = sizes, colors = colors, labels = labels)
}


#' Plot contrastogram
#' @param object SummarizedExperiment
#' @param subgroup_colors named color vector (names = subgroups)
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' plot_contrastogram(object)
#'
#' file <- download_data('billing16.proteingroups.txt')
#' invert_subgroups <- c('EM_E', 'BM_E', 'BM_EM')
#' object <- read_proteingroups(file, invert_subgroups = invert_subgroups)
#' plot_contrastogram(object)
#' @export
plot_contrastogram <- function(
  object, contrasts = default_contrasts(object),
  subgroup_colors = default_color_values2(object)
){
    # Initialize
    V2 <- .N <- N <- NULL

    # Perform limma
    object %<>% add_limma()
    contrastogram_matrices <- compute_connections(object, contrasts, subgroup_colors = subgroup_colors)
    sizes  <- contrastogram_matrices$sizes
    colors <- contrastogram_matrices$colors
    labels <- contrastogram_matrices$labels
    widths <- scales::rescale(sizes, c(0.01,30))

    # Plot diagram
    dt <- split_subgroup_levels(object)
    nrow <- dt[, data.table::uniqueN(V2)]
    nperrow <- dt[, .N, by = 'V1'][, N]
    #dir.create('~/importomicscache/contrastogram')
    #pdf('~/importomicscache/contrastogram/directed_contrastogram.pdf',
    #width = 9, height = 9)
    labels %<>% as.data.frame()
    diagram::plotmat(labels, nperrow, relsize = 1, box.size = 0.05,
    #diagram::plotmat(sizes, nperrow, relsize = 1, box.size = 0.05,
        name = rownames(sizes), box.col = subgroup_colors,
        box.type = 'square', arr.lwd = widths, shadow.size =0, # sqrt(sizes)
        arr.lcol = colors, arr.col = colors)
    #, arr.lcol = log2(1+diagram_matrix))
    #dev.off()
}



