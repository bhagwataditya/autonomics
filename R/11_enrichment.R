#' msig collections
#' @param organism 'human' or 'mouse'
#' @examples
#' msigcollections('human')
#' msigcollections('mouse')
#' @export
msigcollections <- function(organism){
    switch(
        organism, 
        human = MSIGHUMAN,
        mouse = MSIGMOUSE
    )
}


MSIGHUMAN <- c(
    positional          = 'C1',                   # positional    : - genes in same cytogenetic band 
    perturbation        = 'C2:CGP',               # curated       : + genes with altered expression upon (chemical/genetic) perturbation
    otherpathway        = 'C2:CP',                #
    biocarta            = 'C2:CP:BIOCARTA',       #                 + genes in same biocarta pathway
    kegglegacy          = 'C2:CP:KEGG_LEGACY',    #                 -                   kegg pathway (legacy)
    keggmedicus         = 'C2:CP:KEGG_MEDICUS',   #                 +           kegg medicus pathway
    pid                 = 'C2:CP:PID',            #                 -        prot interac db pathway
    reactome            = 'C2:CP:REACTOME',       #                 +               reactome pathway
    wiki                = 'C2:CP:WIKIPATHWAYS',   #                 +                   wiki pathway
    mirlegacy           = 'C3:MIR:MIR_LEGACY',    # regulatory    : + genes with shared mir site (mirdb 6.0, 2020)
    mir                 = 'C3:MIR:MIRDB',         #                 -                            (legacy)
    tfbs                = 'C3:TFT:GTRD',          #                 +                   tfb site (Gene Transcr Reg Db, 2021)
    tfbslegacy          = 'C3:TFT:TFT_LEGACY',    #                 -                            (legacy, 2005)
    cancersubtype       = 'C4:3CA',               # computational : ? genes upregulated in tumor subpopulation  (2023)
    cancerneighbourhood = 'C4:CGN',               #                 ? genes in same cancer driver neighbourhood (2005)
    cancercondition     = 'C4:CM',                #                 ? genes with altered expression in (same) cancer condition (2004)
    gobp                = 'C5:GO:BP',             # ontological   : + genes in same biological process
    gocc                = 'C5:GO:CC',             #                 +             cellular compartment
    gomf                = 'C5:GO:MF',             #                 +     with same molecular function
    humanphenotype      = 'C5:HPO',               #                 + with altered expression in (same) human phenotype
    cancerpathway       = 'C6',                   # oncogenic     : ? genes in cancer-dysregulated pathway
    immuneperturbation  = 'C7:IMMUNESIGDB',       # immunologic   : ? genes with altered expression upon (same) immune system perturbation
    vaccination         = 'C7:VAX',               #                 ? genes with altered expression upon (same) vaccination
    celltype            = 'C8',                   # celltype      : + genes with celltype-specific expression
    hallmark            = 'H'                     # hallmark      : + shared hallmark
 )

MSIGMOUSE <- c(
    positional          = 'M1',                   # positional    : - genes in same cytogenetic band 
    perturbation        = 'M2:CGP',               # curated       : + genes with altered expression upon (chemical/genetic) perturbation
    biocarta            = 'M2:CP:BIOCARTA',       #                 + genes in same biocarta pathway
    reactome            = 'M2:CP:REACTOME',       #                 +               reactome pathway
    wiki                = 'M2:CP:WIKIPATHWAYS',   #                 +                   wiki pathway
    tfbs                = 'M3:GTRD',              #                 +                   tfb site (Gene Transcr Reg Db, 2021)
    mir                 = 'M3:MIRDB',             #                 -                            (legacy)
    gobp                = 'M5:GO:BP',             # ontological   : + genes in same biological process
    gocc                = 'M5:GO:CC',             #                 +             cellular compartment
    gomf                = 'M5:GO:MF',             #                 +     with same molecular function
    tumorphenotype      = 'M5:MPT',               #                 + tumor phenotype ontology 
    celltype            = 'C8',                   # celltype      : + genes with celltype-specific expression
    hallmark            = 'MH'                    # hallmark      : + shared hallmark
)

#' Get path to msigfile
#' @param organism 'human' or 'mouse'
#' @param year      number
#' @param release   number
#' @examples
#' msigfile('human')
#' msigfile('mouse')
#' @export
msigfile <- function(organism, year = 2023, release = 2){
    assert_scalar_subset(organism, c('human', 'mouse'))
    dir <- file.path(R_user_dir('autonomics', 'cache'), 'msigdb')
    file <- sprintf('msigdb_v%d.%d.%s.db', year, release, switch(organism, human = 'Hs', mouse = 'Mm'))
    file <- file.path(dir, file)
    file
}

#' Read msig datatable
#' @param organism    'human' or 'mouse'
#' @param file         sqlite file
#' @param collections  subset of names(MSIGCOLLECTIONS)
#' @examples
#' dir <- file.path(tools::R_user_dir('autonomics', 'cache'), 'msigdb')
#' read_msigdt(file = file.path(dir, 'msigdb_v2023.2.Hs.db'))
#' read_msigdt(file = file.path(dir, 'msigdb_v2023.2.Mm.db'))
#' @export
read_msigdt <- function(
    organism    = 'human',
    file        = msigfile(organism), 
    collections = c('gobp', 'gomf', 'gocc', 'reactome', 'wiki')
){
# Assert
    if (!requireNamespace('DBI',     quietly = TRUE))  message("BiocManager::install('DBI'). Then re-run.")
    if (!requireNamespace('RSQLite', quietly = TRUE))  message("BiocManager::install('RSQLite'). Then re-run.")
    if (!file.exists(file))  return(NULL)
    gene_set_id <- gene_symbol_id <- id <- symbol <- NULL
    standard_name <- collection_name <- collection <- NULL
# Read
    con <- DBI::dbConnect(RSQLite::SQLite(), file)
    dt1 <- data.table(DBI::dbReadTable(con, 'gene_set_gene_symbol'))
    dt2 <- data.table(DBI::dbReadTable(con, 'gene_symbol'         ))
    dt3 <- data.table(DBI::dbReadTable(con, 'gene_set'            ))
    dt1 %<>% extract(, .(setid  = gene_set_id, geneid = gene_symbol_id))
    dt2 %<>% extract(, .(geneid = id, gene = symbol))
    dt3 %<>% extract(, .(setid  = id, set = standard_name, collection = collection_name))
    DBI::dbDisconnect(con)
    
    msigdt <- dt1
    msigdt %<>% merge(dt2, by = 'geneid')
    msigdt %<>% merge(dt3, by = 'setid')
    msigdt %<>% extract(, c('collection', 'set', 'gene'), with = FALSE)
# Return
    msigdt %<>% extract(collection %in% msigcollections(organism)[collections])
    msigdt
}

# # Lower-level function
#     # selected 
#         file <- download_data('atkin.somascan.adat')
#         object <- read_somascan(file)
#         detected <- fdt(object)$EntrezGeneSymbol %>% split_extract_fixed(' ', 1) %>% unique()
#     # universe and set
#         coldt <- read_msigdt()
#         universe <- unique(coldt[, gene])  # 19 669
#         set <- coldt[set == 'GOBP_GLYCOLYTIC_PROCESS', gene]
#     # enrichment
#         .enrichdt <- .enrichment(detected, set)
.enrichment <- function(selected, set, universe){
    dt <- data.table(`in`         = nintersect(universe, set),
                      in.selected = nintersect(selected, set),
                      out         =   nsetdiff(universe, set),
                      selected    =     length(selected))
    dt[, p.selected := 1 - phyper(in.selected, `in`, out, selected)]
    dt[]
}

.enrichmentVERBOSE <- function(selected, set, universe){
    dt <- data.table(`in`         = nintersect(universe, set),
                      in.selected = nintersect(selected, set),
                      out         =   nsetdiff(universe, set),
                      selected    =     length(selected),
                in.selected.genes = collapse(intersect(selected, set), ' '))
    dt[, p.selected := 1 - phyper(in.selected, `in`, out, selected)]
    dt[]
}



#' Analyze enrichment
#' 
#' Are up / downregulated genes enriched in geneset?
#' 
#' @param object   \code{SummarizedExperiment}
#' @param coldt    \code{data.table}, e.g. \code{\link{read_msigdt}}
#' @param by       \code{fvar}
#' @param sep      \code{string} or \code{NULL}
#' @param contrast \code{string} in \code{coefs(object)}
#' @param fit      \code{'limma'}, \code{'lm'}, \code{'lme'}, \code{'lmer'}, \code{'wilcoxon'}
#' @param p       pvalue cutoff
#' @param n       no of detected genes required (for geneset to be examined)
#' @param verbose TRUE or FALSE
#' @param fast    TRUE or FALSE : use fast implementation (or modular)?

#' @examples
#' file <- download_data('atkin.somascan.adat')
#' object <- read_somascan(file, fit = 'limma')
#' fvars(object) %<>% gsub('EntrezGeneSymbol', 'gene', .)
#' coldt <- read_msigdt(collections = 'gobp')
#' coldt
#' enrichdt  <- enrichment(object, coldt, by = 'gene', sep = ' ')
#' enrichdt2 <- enrichment(object, coldt, by = 'gene', sep = ' ', fast = FALSE)
#' @details
#' Four enrichment analyses per geneset using the Fisher Exact Test (see four pvalues).
#' Results are returned in a data.table
#' \tabular{rl}{
#'    in                \tab :                  genes      in pathway \cr
#'    in.detected       \tab :         detected genes      in pathway \cr
#'    in.updown         \tab : up/downregulated genes      in pathway \cr
#'    in.up(.genes)     \tab :      upregulated genes      in pathway \cr
#'    in.down(.genes)   \tab :    downregulated genes      in pathway \cr
#'    out               \tab :                  genes outside pathway \cr
#'    detected          \tab :         detected genes (in + out)      \cr
#'    updown            \tab : up/downregulated genes (in + out)      \cr
#'    up                \tab :      upregulated genes (in + out)      \cr
#'    down              \tab :    downregulated genes (in + out)      \cr
#'    p.contrast.up     \tab : prob to randomly select this many (or more)   upregulated genes (among detected genes)       \cr
#'    p.contrast.down   \tab : prob to randomly select this many (or more) downregulated genes (among detected genes)       \cr
#'    p.contrast.updown \tab : prob to randomly select this many (or more) up OR downregulated genes (among detected genes) \cr
#'    p.contrast.naive  \tab : prob to randomly select this many (or more) up OR downregulated genes (among all genes)      \cr
#'    p.detected        \tab : prob to randomly select this many (or more) detected genes (among all genes)
#' }
#' @importFrom stats phyper
#' @export
enrichment <- function(
    object, 
    coldt, 
    by       = 'gene', 
    sep, 
    contrast = default_coefs(object)[1], 
    fit      = fits(object)[1], 
    p        = 0.05, 
    n        = 3, 
    verbose  = TRUE, 
    fast     = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    if (is.null(coldt)) return(object)
    assert_is_data.table(   coldt)
    assert_scalar_subset(by, fvars(object))
    assert_scalar_subset(by, names(coldt ))
    assert_all_are_non_missing_nor_empty_character(fdt(object)[[by]])
    assert_all_are_non_missing_nor_empty_character(      coldt[[by]])
    if (!is.null(sep))  assert_is_a_string(sep)
    assert_scalar_subset(contrast, coefs(object))
    assert_scalar_subset(fit, fits(object))
    genes0 <- fdt(object)[[by]]
    fdt(object)[[by]] %<>% split_extract_fixed(sep, 1)
    gene <- in.up <- in.detected <- in.down <- in.updown <- `in` <- out <- NULL
# Constants
          all <- unique(fdt(object)[[by]]) %>% union(coldt[, unique(gene)])                  # 17 987  all
     detected <- unique(fdt(object)[[by]])                                                   #  1 084  detected
       updown <-    pfeatures(object, fit = fit, coef = contrast, fvar = by, p = p)          #     59  updown
           up <-   upfeatures(object, fit = fit, coef = contrast, fvar = by, p = p)          #     14  up
         down <- downfeatures(object, fit = fit, coef = contrast, fvar = by, p = p)          #     45  down
    NOTdetected <- setdiff(all, detected)                                                    # 16 903  not detected
    NOTupdown   <- setdiff(detected, updown)                                                 #   1025  not updown
    NOTup       <- setdiff(detected, up)                                                     #   1070  not up
    NOTdown     <- setdiff(detected, down)                                                   #   1039  not down
    upvar     <- sprintf('p.%s.upDETECTED',     contrast)
    downvar   <- sprintf('p.%s.downDETECTED',   contrast)
    updownvar <- sprintf('p.%s.updownDETECTED', contrast)
    naivevar  <- sprintf('p.%s.updownGENOME',   contrast)
    detectedvar <- 'p.detectedGENOME'
# Message
    if (verbose){ 
        if (fast)  cmessage('\tAnalyze enrichment (fast)') else cmessage('\tAnalyze enrichment (modular)')
        cmessage( '\t\tIs pathway enriched in       UP genes (among DETECTED genes) ?')
        cmessage( '\t\tIs pathway enriched in     DOWN genes (among DETECTED genes) ?')
        cmessage( '\t\tIs pathway enriched in  UP/DOWN genes (among DETECTED genes) ?')
        cmessage( '\t\tIs pathway enriched in  UP/DOWN genes (among GENOME   genes) ?')
        cmessage( '\t\tIs pathway enriched in DETECTED genes (among GENOME   genes) ?')  }
# Enrichment
    if (fast){
        enrichdt <- coldt[,  .( `in`             =          nintersect(all,      gene),       # in
                                 in.detected     =          nintersect(detected, gene),       # in: detected
                                 in.updown       =          nintersect(updown,   gene),       # in: updown
                                 in.up           =          nintersect(up,       gene),       # in: up
                                 in.down         =          nintersect(down,     gene),       # in: down
                                 in.up.genes     =  collapse(intersect(up,       gene), ' '), #
                                 in.down.genes   =  collapse(intersect(down,     gene), ' '), #
                                 out             =            nsetdiff(all,      gene),       # out
                                 detected        = length(detected),                          # detected
                                 updown          = length(updown),                            # updown
                                 up              = length(up),                                # up
                                 down            = length(down)),                             # down
                        by = 'set'  ]
        enrichdt[, (      upvar) := 1- phyper(in.up,        in.detected, detected - in.detected,     up),  by = 'set']
        enrichdt[, (    downvar) := 1- phyper(in.down,      in.detected, detected - in.detected,   down),  by = 'set']
        enrichdt[, (  updownvar) := 1- phyper(in.updown,    in.detected, detected - in.detected, updown),  by = 'set']
        enrichdt[, (   naivevar) := 1- phyper(in.updown,   `in`,         out,                    updown),  by = 'set']
        enrichdt[, (detectedvar) := 1- phyper(in.detected, `in`,         out,                  detected),  by = 'set']
    } else {
              upDT <- coldt[, .enrichmentVERBOSE(  up, gene, detected), by = 'set']
            downDT <- coldt[, .enrichmentVERBOSE(down, gene, detected), by = 'set']
          updownDT <- coldt[, .enrichment(     updown, gene, detected), by = 'set']
          updownGN <- coldt[, .enrichment(     updown, gene,      all), by = 'set']
        detectedGN <- coldt[, .enrichment(   detected, gene,      all), by = 'set']
        assert_are_identical(upDT$set,     downDT$set)
        assert_are_identical(upDT$set,   updownDT$set)
        assert_are_identical(upDT$set,   updownGN$set)
        assert_are_identical(upDT$set, detectedGN$set)
        enrichdt <- data.table( set           =       upDT$set,
                               `in`           = detectedGN$`in`, 
                                in.detected   = detectedGN$in.selected,
                                in.updown     =   updownDT$in.selected,
                                in.up         =       upDT$in.selected,
                                in.down       =     downDT$selected,
                                in.up.genes   =       upDT$in.selected.genes,
                                in.down.genes =     downDT$in.selected.genes,
                                out           = detectedGN$out,
                                detected      = detectedGN$selected,
                                updown        =   updownDT$selected,
                                up            =       upDT$selected,
                                down          =     downDT$selected )
        enrichdt[, (      upvar) :=       upDT$p.selected]
        enrichdt[, (    downvar) :=     downDT$p.selected]
        enrichdt[, (  updownvar) :=   updownDT$p.selected]
        enrichdt[, (   naivevar) :=   updownGN$p.selected]
        enrichdt[, (detectedvar) := detectedGN$p.selected]
    }
# Return
    usefuldt  <- enrichdt[in.detected >= n]
    uselessdt <- enrichdt[in.detected <  n]
    usefuldt %<>% extract(order(get(updownvar)))
    usefuldt[]
}

    nsetdiff <- function(x, y)    length(base::setdiff(  x, y))
  nintersect <- function(x, y)    length(base::intersect(x, y))
    collapse <- function(x, sep)  paste0(x, collapse = sep)
