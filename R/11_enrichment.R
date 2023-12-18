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
#' dir <- file.path(R_user_dir('autonomics', 'cache'), 'msigdb')
#' read_msigdt(file = file.path(dir, 'msigdb_v2023.2.Hs.db'))
#' read_msigdt(file = file.path(dir, 'msigdb_v2023.2.Mm.db'))
#' @importFrom DBI      dbConnect  dbReadTable
#' @importFrom RSQLite  SQLite
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
# Read
    con <- dbConnect(SQLite(), file)
    dt1 <- data.table(dbReadTable(con, 'gene_set_gene_symbol'))[, .(setid  = gene_set_id, geneid     = gene_symbol_id)]
    dt2 <- data.table(dbReadTable(con, 'gene_symbol'         ))[, .(geneid = id,          gene = symbol        )]
    dt3 <- data.table(dbReadTable(con, 'gene_set'            ))[, .(setid  = id, set = standard_name, collection = collection_name)]
    msigdt <- dt1
    msigdt %<>% merge(dt2, by = 'geneid')
    msigdt %<>% merge(dt3, by = 'setid')
    msigdt %<>% extract(, c('collection', 'set', 'gene'), with = FALSE)
# Return
    msigdt %<>% extract(collection %in% msigcollections(organism)[collections])
    msigdt
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
#' @param p        p value cutoff
#' @param n        no of detected genes required (for geneset to be examined)
#' @examples
#' file <- download_data('atkin.somascan.adat')
#' object <- read_somascan(file, fit = 'limma')
#' fvars(object) %<>% stri_replace_first_fixed('EntrezGeneSymbol', 'gene')
#' coldt <- read_msigdt(collections = 'gobp')
#' coldt
#' enrichment(object, coldt, by = 'gene', sep = ' ')
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
#' @export
enrichment <- function(
    object, 
    coldt, 
    by       = 'gene', 
    sep, 
    contrast = default_coefs(object)[1], 
    fit      = fits(object)[1], 
    p        = 0.05, 
    n        = 3
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
# Genome
          all <- unique(fdt(object)[[by]]) %>% union(coldt[, unique(gene)])                  # 17 987  all
     detected <- unique(fdt(object)[[by]])                                                   #  1 084  detected
       updown <-    pfeatures(object, fit = fit, coef = contrast, fvar = by, p = p)          #     59  updown
           up <-   upfeatures(object, fit = fit, coef = contrast, fvar = by, p = p)          #     14  up
         down <- downfeatures(object, fit = fit, coef = contrast, fvar = by, p = p)          #     45  down
    NOTdetected <- setdiff(all, detected)                                                    # 16 903  not detected
    NOTupdown   <- setdiff(detected, updown)                                                 #   1025  not updown
    NOTup       <- setdiff(detected, up)                                                     #   1070  not up
    NOTdown     <- setdiff(detected, down)                                                   #   1039  not down
    
      enrichdt <- coldt[,  .(                                                                 #
                                `in`             =          nintersect(all,      gene),       # in
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
    upvar     <- sprintf('p.%s.up',     contrast)
    downvar   <- sprintf('p.%s.down',   contrast)
    updownvar <- sprintf('p.%s.updown', contrast)
    naivevar  <- sprintf('p.%s.naive',  contrast)
    detectedvar <- 'p.detected'
    enrichdt[, (      upvar) := 1- phyper(in.up,        in.detected, detected - in.detected, updown),  by = 'set']
    enrichdt[, (    downvar) := 1- phyper(in.down,      in.detected, detected - in.down,     updown),  by = 'set']
    enrichdt[, (  updownvar) := 1- phyper(in.updown,    in.detected, detected - in.up,       updown),  by = 'set']
    enrichdt[, (   naivevar) := 1- phyper(in.updown,   `in`,         out,                    updown),  by = 'set']
    enrichdt[, (detectedvar) := 1- phyper(in.detected, `in`,         out,                  detected),  by = 'set']
    usefuldt  <- enrichdt[in.detected >= n]
    uselessdt <- enrichdt[in.detected <  n]
    usefuldt %<>% extract(order(get(updownvar)))
    usefuldt[]
}

    nsetdiff <- function(x, y)    length(base::setdiff(  x, y))
  nintersect <- function(x, y)    length(base::intersect(x, y))
    collapse <- function(x, sep)  paste0(x, collapse = sep)
