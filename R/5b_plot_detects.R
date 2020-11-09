#==============================================================================
#
#                        plot_detects_per_subgroup
#
#==============================================================================

#' Plot detects per subgroup
#'
#' Plot number of detects, partial detects, and nondetects per subgroup
#'
#' @param object SummarizedExperiment
#' @param group  svar (symbol)
#' @return ggplot object
#' @examples
#' # STEMCELLS
#'     file <- download_data('billing16.proteingroups.txt')
#'     invert_subgroups <- c('E_EM', 'E_BM', 'EM_BM')
#'     object <- read_proteingroups(file, invert_subgroups = invert_subgroups,
#'                                   impute = FALSE, plot = FALSE)
#'     plot_detects_per_subgroup(object)
#' @export
plot_detects_per_subgroup <- function(object, group = subgroup){
# Assert
    assert_is_all_of(object, 'SummarizedExperiment')
    group_var <- as_string(ensym(group))
    assert_is_subset(group_var, svars(object))
    variable <- subgroup <- value <- . <- NULL
# Prepare datatable
    split_objects  <- split_by_svar(object, group_var)
    imps     <- vapply(split_objects, count_imputes,    numeric(1))
    nons     <- vapply(split_objects, count_nondetects, numeric(1))
    partials <- vapply(split_objects, count_pardetects, numeric(1)) - nons
    fulls    <- nrow(object) - partials - imps - nons
    plot_dt  <- data.table(subgroup = factor(names(split_objects),
                                            slevels(object, group_var)))
    if (any(fulls   !=0))            plot_dt %<>% extract(, fulls   := fulls)
    if (any(partials!=0))            plot_dt %<>% extract(, partials:= partials)
    if (any(nons!=0) | any(imps!=0)) plot_dt %<>% extract(, nons := nons + imps)
    plot_dt %<>% data.table::melt(id.vars = 'subgroup')
# Set order of variable levels
    variable_levels <- c('nons', 'partials', 'fulls') %>%
                        extract(.%in% plot_dt$variable)
    plot_dt$variable %<>% factor(variable_levels)
# Plot
    title <- 'fulldetects  |  partialdetects'
    if ( all(imps==0) & !all(nons==0))  title %<>% paste0('  |  nondetects')
    if (!all(imps==0) &  all(nons==0))  title %<>% paste0('  |  imputes')
    plot_dt$subgroup %<>% factor(rev(levels(.)))

    plot_data(plot_dt, geom_col,
                x = subgroup, y = value, fill = subgroup, group = variable,
                fixed = list(position = position_stack(), color="black")) +
    #ggplot(plot_dt, aes(x = subgroup, y = value, fill = subgroup,
    #                    group = variable)) +
    ggtitle(title) + theme_bw() +
    #geom_col(color = 'black', position = position_stack()) +
    #scale_fill_manual(values = colorscale) +
    geom_text(aes(label=value, x = subgroup),
                    position = position_stack(vjust=0.5),
                    size = rel(3)) +
    theme(#axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust=1),
                panel.grid.major.x = element_blank(),
                panel.border       = element_blank(),
                plot.title         = element_text(hjust = 0.5)) +
    xlab(NULL) +
    ylab(NULL) +
    coord_flip()
}

count_imputes <- function(object){
    sum(rowAlls(is_imputed(object)))
}

count_nondetects <- function(object){
    sum(rowAlls(is.na(exprs(object))))
}

count_pardetects <- function(object){
    sum(rowAnys(is.na(exprs(object))))
}

