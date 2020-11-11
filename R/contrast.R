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
aggregate_contrasts <- function(contrastmat, dim){
    apply(contrastmat,
          dim,
          function(x){
              paste0(sprintf('(%s)/%d', x, length(x)), collapse = ' + ')})
}

aggregate_column_contrasts <- function(contrastmat){
    aggregate_contrasts(contrastmat, 2)
}

aggregate_row_contrasts <- function(conc_contrast_matrix){
    aggregate_contrasts(conc_contrast_matrix, 1)
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

   # If contrasts present in object, use these
   if (!is.null(autonomics.import::contrastdefs(object))){
      message('\tUse contrasts in autonomics.import::contrastdefs(object)')
      return(autonomics.import::contrastdefs(object))
   }

   # Extract subgroup levels
   subgroup_values <- autonomics.import::sdata(object)$subgroup
   if (is.character(subgroup_values))  subgroup_values %<>%
                                            autonomics.support::factorify()

   # Ratios themselves for ratio data
   if (autonomics.import::contains_ratios(object)){
      message('\tGenerate contrasts from ratios')
      contrasts <- levels(subgroup_values) %>%
                    magrittr::set_names(validify_contrast_names(.))

   # Difference contrasts for abundance data
   } else {
      message('\tGenerate difference contrasts')
      contrasts <- create_diff_contrasts(object)
   }

   # Return
   contrasts
}

# add_limma2 <- function(object, contrastdefs = create_diff_contrasts(object)){
#     S4Vectors::metadata(object)$contrastdefs <- contrastdefs
#     object %<>% autonomics.find::add_limma(contrastdefs = contrastdefs)
#     object
# }


#===============================================================================


#' Plot contrastogram
#' @param object SummarizedExperiment
#' @param directed TRUE or FALSE: whether to distinguish up and downregulations
#' @param subgroup_colors named color vector (names = subgroups)
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file)
#' plot_contrastogram(object)
#' plot_contrastogram(object, directed = FALSE)
#' @export
plot_contrastogram <- function(
    object, directed = TRUE, subgroup_colors = default_color_values2(object),
    sparse = FALSE
){
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
    #pdf('~/importomicscache/contrastogram/directed_contrastogram.pdf', width = 9, height = 9)
    diagram::plotmat(connection_sizes, nperrow, relsize = 1, box.size = 0.05,
        name = rownames(connection_sizes), box.col = subgroup_colors,
        box.type = 'square', arr.lwd = widths, # sqrt(connection_sizes)
        arr.lcol = connection_colors, arr.col = connection_colors) #, arr.lcol = log2(1+diagram_matrix))
    #dev.off()
}


default_color_values2 <- function(object){
    autonomics.plot::default_color_values(object)[autonomics.import::subgroup_levels(object)]
}


compute_connections <- function(
    object, directed = FALSE, subgroup_colors = default_color_values2(object)
){
# subgroup matrix, difference contrasts, limma
    pvalues <- autonomics.import::limma(object)[, , 'p']
    effects <- autonomics.import::limma(object)[, , 'effect']
    is_significant <-  pvalues < 0.05
    is_up          <- (pvalues < 0.05) & (effects > 0)
    is_down        <- (pvalues < 0.05) & (effects < 0)

    n <- apply(if (directed) is_up else is_significant, 2,sum, na.rm = TRUE)
    ndown <- apply(is_down, 2, sum, na.rm = TRUE)
# Create diagram
    subgrouplevels <- autonomics.import::subgroup_levels(object)
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

