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


add_contrastmat <- function(
    object,
    contrasts = default_contrasts(object),
    design    = design(object)
){
    contrasts(object) <- create_contrastmat(object, contrasts, design)
    object
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
add_limma <- function(
    object,
    contrasts = default_contrasts(object),
    formula = if (single_subgroup(object)) ~ 1 else ~ 0 + subgroup
){
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
        correlation  <- exprs(object) %>%
                        duplicateCorrelation(design, block = block) %>%
                        extract2('consensus.correlation')}
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


#=============================================================================

#' Matrixify subgroups
#'
#' Organize subgroups in 2d matrix
#'
#' @param object SummarizedExperiment
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' matrixify_subgroups(object)
#' @noRd
matrixify_subgroups <- function(object){
    dt <- split_subgroup_levels(object)
    subgroup_matrix <- as.matrix(data.table::dcast(
        dt, V1 ~ V2, value.var = 'subgroup'), rownames = 'V1')
    subgroup_matrix
}

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


#=============================================================================


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
#' subgroup_matrix <- matrixify_subgroups(object)
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
#' subgroup_matrix <- matrixify_subgroups(object)
#' contrast_rows(subgroup_matrix)
#' @noRd
contrast_rows <- function(subgroup_matrix, symbol = ' - '){

    contrastmat <- matrix(  sprintf('%s%s%s',
                                  subgroup_matrix[-1, ],
                                  symbol,
                                  subgroup_matrix[-nrow(subgroup_matrix), ]),
                            nrow = nrow(subgroup_matrix)-1,
                            ncol = ncol(subgroup_matrix))

    colnames(contrastmat) <- colnames(subgroup_matrix)
    rownames(contrastmat) <- sprintf('%s - %s',
                            rownames(subgroup_matrix)[-1],
                            rownames(subgroup_matrix)[-nrow(subgroup_matrix)])
    contrastmat
}


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
#' (x <- matrixify_subgroups(object))
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

aggregate_column_contrasts <- function(contrastmat){
    aggregate_contrasts(contrastmat, 2)
}

aggregate_row_contrasts <- function(conc_contrastmat){
    aggregate_contrasts(conc_contrastmat, 1)
}


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

    # Subgroup matrix
    subgroup_matrix <- matrixify_subgroups(object)

    # Contrast matrix
    column_contrastmat   <- contrast_columns(subgroup_matrix)
    column_contrastnames <- contrast_columns(subgroup_matrix, '__')
    row_contrastmat      <- contrast_rows(   subgroup_matrix)
    row_contrastnames    <- contrast_rows(   subgroup_matrix, '__')

    # Contrast vector
    column_contrasts <- structure(c(column_contrastmat),
                                    names = c(column_contrastnames))
    row_contrasts    <- structure( c(row_contrastmat),
                                    names = c(row_contrastnames))
    #aggregated_column_contrasts <- aggregate_column_contrasts(
    #                                    column_contrastmat)
    #aggregated_row_contrasts    <- aggregate_row_contrasts(
    #                                    row_contrastmat)

    # Return
    c(  column_contrasts,
        row_contrasts)#,
        #3aggregated_column_contrasts,
        #aggregated_row_contrasts)
}

#==============================================================================

#' Default contrasts
#' @param object SummarizedExperiment
#' @return named character vector: contrast definitions
#' @examples
#' # STEM CELL COMPARISON
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' default_contrasts(object)
#'
#' # GLUTAMINASE
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' default_contrasts(object)
#' @export
default_contrasts <- function(object){

    # Extract subgroup levels
    subgroup_values <- sdata(object)$subgroup
    if (is.character(subgroup_values))  subgroup_values %<>% factorify()

    # Ratios themselves for ratio data
    if (contains_ratios(object)){
        message('\tGenerate contrasts from ratios')
        contrasts <- levels(subgroup_values) %>%
                    set_names(validify_contrast_names(.))

    # Difference contrasts for abundance data
    } else {
        message('\tGenerate difference contrasts')
        contrasts <- create_diff_contrasts(object)
    }

    # Return
    contrasts
}


validify_contrast_names <- function(x){
   x %>% gsub(' ', '',  ., fixed = TRUE) %>%
         gsub('-', '_', ., fixed = TRUE) %>%
         make.names()
}

# add_limma2 <- function(object, contrasts = create_diff_contrasts(object)){
#     S4Vectors::metadata(object)$contrasts <- contrasts
#     object %<>% autonomics.find::add_limma(contrasts = contrasts)
#     object
# }


#==============================================================================


#' Plot contrastogram
#' @param object SummarizedExperiment
#' @param directed TRUE or FALSE: whether to distinguish up and downregulations
#' @param subgroup_colors named color vector (names = subgroups)
#' @param sparse TRUE or FALSE
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' plot_contrastogram(object)
#' plot_contrastogram(object, directed = FALSE)
#' @export
plot_contrastogram <- function(object, directed = TRUE,
    subgroup_colors = default_color_values2(object), sparse = FALSE
){
    # Initialize
    V2 <- .N <- N <- NULL

    # Perform limma
    object %<>% add_limma()
    contrastogram_matrices <- compute_connections(
        object, directed = directed, subgroup_colors = subgroup_colors)
    connection_sizes  <- contrastogram_matrices$connection_sizes
    connection_colors <- contrastogram_matrices$connection_colors

        widths <- scales::rescale(connection_sizes, c(0.01,30))
        if (sparse) connection_colors[connection_sizes/nrow(object)<0.50] <- "0"

    # Plot diagram
    dt <- split_subgroup_levels(object)
    nrow <- dt[, data.table::uniqueN(V2)]
    nperrow <- dt[, .N, by = 'V1'][, N]
    #dir.create('~/importomicscache/contrastogram')
    #pdf('~/importomicscache/contrastogram/directed_contrastogram.pdf',
    #width = 9, height = 9)
    diagram::plotmat(connection_sizes, nperrow, relsize = 1, box.size = 0.05,
        name = rownames(connection_sizes), box.col = subgroup_colors,
        box.type = 'square', arr.lwd = widths, # sqrt(connection_sizes)
        arr.lcol = connection_colors, arr.col = connection_colors)
    #, arr.lcol = log2(1+diagram_matrix))
    #dev.off()
}

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
    default_color_values(object)[subgroup_levels(object)]
}


compute_connections <- function(
    object, directed = FALSE, subgroup_colors = default_color_values2(object)
){
# subgroup matrix, difference contrasts, limma
    pvalues <- limma(object)[, , 'p']
    effects <- limma(object)[, , 'effect']
    is_significant <-  pvalues < 0.05
    is_up          <- (pvalues < 0.05) & (effects > 0)
    is_down        <- (pvalues < 0.05) & (effects < 0)

    n <- apply(if (directed) is_up else is_significant, 2,sum, na.rm = TRUE)
    ndown <- apply(is_down, 2, sum, na.rm = TRUE)
# Create diagram
    subgrouplevels <- subgroup_levels(object)
    connection_sizes <- connection_colors <- matrix(
        0,
        nrow = length(subgrouplevels),
        ncol = length(subgrouplevels),
        dimnames = list(subgrouplevels, subgrouplevels))
# Add column contrasts
    subgroup_matrix <- matrixify_subgroups(object)
    for (i in 1:nrow(subgroup_matrix)){
        for (j in 2:ncol(subgroup_matrix)){
            from <- subgroup_matrix[i, j-1]
            to   <- subgroup_matrix[i, j]
            connection_sizes[to, from]         <- n[[paste0(to, '__', from)]]
            connection_colors[to, from] <- subgroup_colors[[to]]
            if (directed){
                connection_sizes[from, to] <- ndown[[paste0(to, '__', from)]]
                connection_colors[from, to] <- subgroup_colors[[to]]
            }
        }
    }
# Add row contrasts
    for (j in 1:ncol(subgroup_matrix)){
        for (i in 2:ncol(subgroup_matrix)){
            from <- subgroup_matrix[i-1, j]
            to   <- subgroup_matrix[i, j]
            connection_sizes[to, from]         <- n[[paste0(to, '__', from)]]
            connection_colors[to, from] <- subgroup_colors[[to]]
            if (directed){
                connection_sizes[from, to] <- ndown[[paste0(to, '__', from)]]
                connection_colors[from, to] <- subgroup_colors[[to]]
            }
        }
    }
# Return
    list(connection_sizes = connection_sizes,
        connection_colors = connection_colors)
}





#==============================================================================
#' Get/Set contrasts
#' @param object SummarizedExperiment
#' @param value named string vector (see examples)
#' @return updated SummarizedExperiment
#' @examples
#' require(magrittr)
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_proteingroups(file)
#' contrasts1 <- c(EM_E = 'EM_E', BM_E = 'BM_E', BM_EM = 'BM_EM')
#' contrasts(object) <-     contrasts1   # conventional setter
#' object %>% set_contrasts(contrasts1)  # piping       setter
#' contrasts(object)                        # getter
#' @rdname contrasts
#' @export
setGeneric("contrasts",   function(object)   standardGeneric("contrasts"))

#' @rdname contrasts
#' @export
setMethod(
    "contrasts",
    signature("SummarizedExperiment"),
    function(object) metadata(object)$contrasts )

#' @rdname contrasts
#' @export
setGeneric(
    "contrasts<-",
    function(object, value)  standardGeneric("contrasts<-") )

#' @rdname contrasts
setReplaceMethod(
    "contrasts",
    signature("SummarizedExperiment", "character"),
    function(object, value){ metadata(object)$contrasts <- value; object})

#' @rdname contrasts
setReplaceMethod(
    "contrasts",
    signature("SummarizedExperiment", "NULL"),
    function(object, value){object})

#' @rdname contrasts
#' @export
set_contrasts <- function(object, value){
    contrasts(object) <- value
    object
}


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



