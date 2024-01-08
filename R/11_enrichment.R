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


utils::globalVariables(c('in', 'in.selected', 'out', 'selected', 'p.selected'))
# # Lower-level function
#     # selected 
#         file <- download_data('atkin.somascan.adat')
#         object <- read_somascan(file)
#         detected <- fdt(object)$EntrezGeneSymbol %>% split_extract_fixed(' ', 1) %>% unique()
#     # universe and set
#         pathwaydt <- read_msigdt()
#         universe <- unique(pathwaydt[, gene])  # 19 669
#         set <- pathwaydt[set == 'GOBP_GLYCOLYTIC_PROCESS', gene]
#     # enrichment
#         .enrichdt <- .enrichment(detected, set)
.enrichment <- function(selected, pathway, universe){
    dt <- data.table(`in`         =  count_in(universe, pathway),
                      in.selected =  count_in(selected, pathway),
                      out         = count_out(universe, pathway),
                      selected    = length(selected))
    dt[, p.selected := 1 - phyper(in.selected, `in`, out, selected)]
    dt[]
}

.enrichmentVERBOSE <- function(selected, pathway, universe){
    dt <- data.table(`in`         =  count_in(universe, pathway),
                      in.selected =  count_in(selected, pathway),
                      out         = count_out(universe, pathway),
                      selected    =     length(selected),
                in.selected.genes = collapse(intersect(selected, pathway), ' '))
    dt[, p.selected := 1 - phyper(in.selected, `in`, out, selected)]
    dt[]
}


#' Abstract model fit
#' @param object            SummarizedExperiment
#' @param fit               character vector
#' @param coef              character vector
#' @param significancevar  'p' or 'fdr'
#' @param significance      fraction : pvalue cutoff
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file, fit = 'limma', coef = 't3')
#' object %<>% abstract_fit()
#' fdt(object)
#' @export
abstract_fit <- function(
    object, fit = fits(object), coef = coefs(object), 
    significancevar = 'p', significance = 0.05
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(fit,   fits(object))
    assert_is_subset(coef, coefs(object))
# Abstract
    for ( curfit in fit){
    for (curcoef in coef){
        abstractvar <- paste(curcoef, curfit, sep = FITSEP)
            pvalues <- modelvec(object, 'p',      fit = curfit, coef = curcoef)
       effectvalues <- modelvec(object, 'effect', fit = curfit, coef = curcoef)
        fdt(object)[[ abstractvar ]] <- 'flat'
        fdt(object)[[ abstractvar ]][ pvalues<significance  &  effectvalues<0 ] <- 'down' 
        fdt(object)[[ abstractvar ]][ pvalues<significance  &  effectvalues>0 ] <- 'up' 
        fdt(object)[[ abstractvar ]] %<>% factor(c('flat', 'up', 'down'))
    }}
    object
}


#' group by level
#' @param x named logical/character/factor
#' @param ... S3 dispatch
#' @return unnamed character
#' @examples
#' t1 <- c( KLF5 = 'up',  F11 = 'up', RIG = 'flat',   ABT1 = 'down')
#' dt <- data.table( gene = c( 'KL5', 'F11', 'RIG',  'ABT1' ), 
#'                     t1 = c( 'up',  'up',  'flat', 'down' ) )
#' group_by_level(t1)                #  character
#' group_by_level(factor(t1))        #     factor
#' group_by_level(dt, 't1', 'gene')  # data.table
#' @export
group_by_level <- function(x, ...)  UseMethod('group_by_level', x)


#' @rdname group_by_level
#' @export
group_by_level.character <- function(x){
    fun <- function(y) names(x)[x == y]
    Map(fun, unique(x))                          
}

#' @rdname group_by_level
#' @export
group_by_level.factor <- function(x){
    y <- x
    y %<>% as.character() 
    y %<>% set_names(names(x))
    group_by_level.character(y)
}

#' @rdname group_by_level
#' @export
group_by_level.data.table <- function(x, var, idvar){
    y <- x[[var]]
    names(y) <- x[[idvar]]
    group_by_level(y)
}

#' logical to factor
#' 
#' switch between factor representations
#' @param x           vector
#' @param var         string
#' @param falseprefix string
#' @examples
#' t1up <- c( TRUE,   FALSE,  TRUE)
#' t1   <- c('flat', 'down', 'up' )  %>%  factor(., .)
#' t1up
#' logical2factor(t1up)
#' factor2logical(t1)
#' @export
logical2factor <- function(
    x, 
    truelevel = get_name_in_parent(x),
    falselevel = paste0('not', truelevel)
){
# Assert
    assert_is_logical(x)
# Convert
    y <- rep(falselevel, length(x))
    y[x] <- truelevel
    y%<>% factor(c(truelevel, falselevel)) # level1 = (noninteresting) baseline
# Return
    y
}

#' @rdname logical2factor
#' @export
factor2logical <- function(x){
# Assert
    assert_is_factor(x)
# Convert
    y <- rep(FALSE, length(x))
    y[x==x[1]] <- TRUE
    y
}


#' Enrichment analysis
#' 
#' Are selected genes enriched in pathway?
#'
#' @param object    \code{SummarizedExperiment}
#' @param pathwaydt pathway \code{data.table}
#' @param var       selection fvar
#' @param levels    selection levels
#' @param genevar   gene fvar
#' @param genesep   gene separator (string)
#' @param verbose   whether to msg
#' @param genes     whether to report genes
#' @examples
#' # Read
#'     file <- download_data('atkin.somascan.adat')
#'     object <- read_somascan(file, fit = 'limma', coefs = 't1')
#'     fvars(object) %<>% gsub('EntrezGeneSymbol', 'gene', .)
#'     object %<>% abstract_fit()
#'     fdt(object)$flat <- fdt(object)$`t1~limma`
#'     fdt(object)$flat %<>% factor2logical() %>% logical2factor('flat', 'updown')
#' # Three usecases
#'     pathwaydt <- read_msigdt(collections = 'gobp')
#'     enrichdt1 <- enrichment(object, pathwaydt, var = abstractvar(object))                                      # 2:n factor 
#'     enrichdt2 <- enrichment(object, pathwaydt, var = 'flat')                                                   #     logical
#'     enrichdt3 <- enrichment(object, pathwaydt, var = abstractvar(object), levels = c('flat', 'down', 'up'))    # 1:n factor
#'     enrichdt4 <- enrichment(object, pathwaydt, var = 'flat', levels = c('flat', 'updown'))                     #     logical
#' # Alternative implementation
#'     enrichdt5 <-  altenrich(object, pathwaydt)   # alternative implementation
#'     cols <- intersect(names(enrichdt1), names(enrichdt5))
#'     all(enrichdt1[, cols, with = FALSE]  ==  
#'         enrichdt5[, cols, with = FALSE])       # identical
#' @details
#' Four enrichment analyses per geneset using the Fisher Exact Test (see four pvalues).
#' Results are returned in a data.table
#' \tabular{rl}{
#'    in                \tab :                  genes      in pathway \cr
#'    in.det            \tab :         detected genes      in pathway \cr
#'    in.sel            \tab : up/downregulated genes      in pathway \cr
#'    in.up(.genes)     \tab :      upregulated genes      in pathway \cr
#'    in.down(.genes)   \tab :    downregulated genes      in pathway \cr
#'    out               \tab :                  genes outside pathway \cr
#'    det               \tab :         detected genes (in + out)      \cr
#'    sel               \tab : up/downregulated genes (in + out)      \cr
#'    up                \tab :      upregulated genes (in + out)      \cr
#'    down              \tab :    downregulated genes (in + out)      \cr
#'    p.coef.upDET      \tab : prob to randomly select this many (or more)   upregulated genes (among detected genes)       \cr
#'    p.coef.downDET    \tab : prob to randomly select this many (or more) downregulated genes (among detected genes)       \cr
#'    p.coef.selDET     \tab : prob to randomly select this many (or more) up OR downregulated genes (among detected genes) \cr
#'    p.coef.selGEN     \tab : prob to randomly select this many (or more) up OR downregulated genes (among genome   genes) \cr
#'    p.detGEN          \tab : prob to randomly select this many (or more) detected genes (among genome genes)
#' }
#' @importFrom stats phyper
#' @export
enrichment <- function(
    object,
    pathwaydt,
    fit      = fits(object)[1],
    coef     = coefs(object, fit = fit)[1],
    var      = abstractvar(object, fit = fit, coef = coef),
    levels   = fdt(object)[[var]] %>% base::levels() %>% extract(-1),
    genevar  = 'gene', 
    genesep  = '[ ,;]',
    n        = 3,
    verbose  = TRUE,
    genes    = FALSE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_data.table(pathwaydt)
    if (!is.null(fit ))  assert_scalar_subset(fit,  fits(object))
    if (!is.null(coef))  assert_scalar_subset(coef, coefs(object))
    assert_scalar_subset(var, fvars(object))
    assert_is_factor(fdt(object)[[var]])
    # assert_is_vector(levels(fdt(object)[[var]]))
    assert_is_subset(levels, levels(fdt(object)[[var]]))
    assert_scalar_subset( genevar, names(pathwaydt))
    assert_all_are_non_missing_nor_empty_character(fdt(object)[[genevar]])
    assert_all_are_non_missing_nor_empty_character(  pathwaydt[[genevar]])
    if (!is.null(genesep))  assert_is_a_string(genesep)
    assert_is_a_bool(verbose)
    if (verbose){
        cmessage('\tAre pathways enriched in `%s` genes ?', var)
        cmessage("\t\tfdt(.)$%s  %%<>%%  split_extract_regex('%s', 1)", genevar, genesep) }
    fdt(object)[[genevar]] %<>% split_extract_regex(genesep, 1)
    `in` <- in.detected <- in.selected <- detected <- selected <- out <- NULL
# Sets
    SETall  <- unique(fdt(object)[[genevar]]) %>% union(pathwaydt[, unique(gene)])
    SETdet  <- unique(fdt(object)[[genevar]])
    SETSsel <- group_by_level(fdt(object), var, genevar) %>% extract(levels)
# Counts
    pad <- function(x) x %>% stringi::stri_pad_left(8)
    if (verbose){
           npathway <- length(unique(pathwaydt$set))
        ncollection <- length(unique(pathwaydt$collection))
       cmessage('\t\t%d pathways from %d collection(s)', npathway, ncollection)
        message( c(
            sprintf(        '\t\tAre %s genes enriched in pathway (among DETECTED) ?', pad(levels[ 1])),
            sprintf(      '\n\t\t    %s genes enriched in pathway (among DETECTED) ?', pad(levels[-1])),
            sprintf('\n\t\t    selected genes enriched in pathway (among DETECTED) ?'), 
            sprintf('\n\t\t    selected genes enriched in pathway (among   GENOME) ?'), 
            sprintf('\n\t\t    detected genes enriched in pathway (among   GENOME) ?')))
    }
    enrichdt <- pathwaydt[ , c(                                 # Tried to drop `(in.)sel` and sum (in.)up and (in.)down later
       `in`    = SETall                 %>%     count_in(gene), # That saves two computations, which always improves speed further.
        in.det = SETdet                 %>%     count_in(gene), # But it ignores the fact that multiple proteins can map to the same gene.
        in.sel = Reduce(union, SETSsel) %>%     count_in(gene), # Which leads to double counting. So rolled back that approach
                 SETSsel                %>%     count_in(gene)       %>% set_names(paste0('in.', names(.))),
                 SETSsel                %>%  collapse_in(gene, ' ')  %>% set_names(paste0('in.', names(.), '.genes')),
        out    = SETall                 %>%    count_out(gene),
        det    = SETdet                 %>%     length(),
        sel    = Reduce(union, SETSsel) %>%     length(),
                 SETSsel                %>%     lapply(length)
    ),  by = 'set' ]
# pvalues
    #one <- length(levels) == 1
    twoplus <- length(levels) >= 2
    nminus  <- length(levels) < length(flevels(object, var))
    invar <- function(levs)  sprintf('in.%s', levs)
                                enrichdt[ , ( sprintf( 'p.%s.DET', levels[1] )) := 1 - phyper( get(invar(levels[1])),  in.det,  det-in.det,  get(levels[1])), by = 'set' ]
                                enrichdt[ , ( sprintf( 'p.%s.GEN', levels[1] )) := 1 - phyper( get(invar(levels[1])),  in.det,         out,  get(levels[1])), by = 'set' ]
    for (level in levels[-1]){  enrichdt[ , ( sprintf( 'p.%s.DET', level     )) := 1 - phyper( get(invar(level    )),  in.det,  det-in.det,  get(level    )), by = 'set' ]
                                enrichdt[ , ( sprintf( 'p.%s.GEN', level     )) := 1 - phyper( get(invar(level    )),  in.det,         out,  get(level    )), by = 'set' ] }
    if (twoplus & nminus)       enrichdt[ , ( sprintf( 'p.sel.DET'           )) := 1 - phyper(                in.sel,  in.det,  det-in.det,  sel           ), by = 'set' ]
    if (twoplus & nminus)       enrichdt[ , ( sprintf( 'p.sel.GEN'           )) := 1 - phyper(                in.sel, `in`,            out,  sel           ), by = 'set' ]
                                enrichdt[ , ( sprintf( 'p.det.GEN'           )) := 1 - phyper(                in.det, `in`,            out,  det           ), by = 'set' ]
# Return                                
    if ( !(twoplus & nminus))   enrichdt[, c('in.sel', 'sel') := NULL]  # n=0 -> p=0 always (1 way  only to choose 0 from 0)
    enrichdt %<>% extract(in.det >= n)                             
    pvar1 <- c('p.sel.DET', sprintf('p.%s.DET', levels[1]))        # NOTE `p.sel.DET` available only when  1 < length(levels) < nlevel
    pvar1 %<>% intersect(names(enrichdt))
    pvar1 %<>% extract(1)                                          
    enrichdt %<>% extract(order(get(pvar1)))
    if (!genes)  enrichdt %<>% extract(, .SD, .SDcols = !patterns('genes'))
    enrichdt[]
        # Note : n=0 -> p=0 always (1 way  only to choose 0 from 0)
}       #        n=1 -> p=0 always (1 way  only to choose 1 from 1)
        #        n=2 -> p=0 often  (2 ways only to choose 2 from 2)
    
#' Alternative Enrichment Analysis
#' @details
#' This is an alternative enrichent analysis implementation.
#' It is more modular: uses four times \code{.enrichment(VERBOSE)?} as backend.
#' But also four times slower than \code{enrichment}, so not recommended.
#' It is retaind for testing purposes.
#' @details This alternative enrichment implementation
#' @param object     \code{SummarizedExperiment}
#' @param pathwaydt  \code{data.table}, e.g. \code{\link{read_msigdt}}
#' @param genevar    \code{gene fvar}
#' @param genesep    \code{string} or \code{NULL}
#' @param coef       \code{string} in \code{coefs(object)}
#' @param fit        \code{'limma'}, \code{'lm'}, \code{'lme'}, \code{'lmer'}, \code{'wilcoxon'}
#' @param significancevar 'p' or 'fdr'
#' @param significance     significance cutoff
#' @param effectsize       effectsize   cutoff
#' @param n        no of detected genes required (for geneset to be examined)
#' @param genes    whether to record genes
#' @param verbose  whether to msg
#' @seealso [enrichment()]
#' @importFrom stats phyper
#' @export
altenrich <- function(
    object, 
    pathwaydt, 
    genevar         = 'gene', 
    genesep         = '[ ,;]',
    coef            = default_coefs(object)[1], 
    fit             = fits(object)[1],
    significancevar = 'p',
    significance    = 0.05,
    effectsize      = 0,
    n               = 3, 
    genes           = FALSE,
    verbose         = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    if (is.null(pathwaydt)) return(object)
    assert_is_data.table(   pathwaydt)
    assert_scalar_subset(genevar, fvars(object))
    assert_scalar_subset(genevar, names(pathwaydt ))
    assert_all_are_non_missing_nor_empty_character(fdt(object)[[genevar]])
    assert_all_are_non_missing_nor_empty_character(  pathwaydt[[genevar]])
    if (!is.null(genesep))  assert_is_a_string(genesep)
    assert_scalar_subset(coef, coefs(object))
    assert_scalar_subset(fit, fits(object))
    genes0 <- fdt(object)[[genevar]]
    fdt(object)[[genevar]] %<>% split_extract_regex(genesep, 1)
    gene <- in.up <- in.detected <- in.down <- in.selected <- `in` <- out <- NULL
# Function
    object %<>% abstract_fit(fit = fit, coef = coef, significancevar = significancevar)
    abstractvar <- abstractvar(object, fit = fit, coef = coef)
# Constants
     all <- unique(fdt(object)[[genevar]])
     all %<>% union(pathwaydt[, unique(gene)])                                 # 17 987  all
     det <- unique(fdt(object)[[genevar]])                                     #  1 084  detected
     sel <- fdt(object)[ get(abstractvar) %in% c('down', 'up') ][[ genevar ]]  #     59  selected
    down <- fdt(object)[ get(abstractvar) %in% c('down'      ) ][[ genevar ]]  #     45  down
      up <- fdt(object)[ get(abstractvar) %in% c('up'        ) ][[ genevar ]]  #     14  up
          upvar <- sprintf('p.%s.up.DET',   coef)
        downvar <- sprintf('p.%s.down.DET', coef)
    selectedvar <- sprintf('p.%s.sel.DET',  coef)
       naivevar <- sprintf('p.%s.sel.GEN',  coef)
    detectedvar <- 'p.det.GEN'
# Message
    if (verbose){ 
        cmessage('\tAnalyze enrichment (modular)')
        cmessage( '\t\tIs pathway enriched in       UP genes (among DETECTED genes) ?')
        cmessage( '\t\tIs pathway enriched in     DOWN genes (among DETECTED genes) ?')
        cmessage( '\t\tIs pathway enriched in  UP/DOWN genes (among DETECTED genes) ?')
        cmessage( '\t\tIs pathway enriched in  UP/DOWN genes (among GENOME   genes) ?')
        cmessage( '\t\tIs pathway enriched in DETECTED genes (among GENOME   genes) ?')  }
# Enrichment
         upDT <- pathwaydt[ , .enrichmentVERBOSE(   up, gene, det ), by = 'set' ]
       downDT <- pathwaydt[ , .enrichmentVERBOSE( down, gene, det ), by = 'set' ]
        selDT <- pathwaydt[ , .enrichment(         sel, gene, det ), by = 'set' ]
        selGN <- pathwaydt[ , .enrichment(         sel, gene, all ), by = 'set' ]
        detGN <- pathwaydt[ , .enrichment(         det, gene, all ), by = 'set' ]
        assert_are_identical(upDT$set, downDT$set)
        assert_are_identical(upDT$set,  selDT$set)
        assert_are_identical(upDT$set,  selGN$set)
        assert_are_identical(upDT$set,  detGN$set)
        enrichdt <- data.table( set           =    upDT$set,
                               `in`           =   detGN$`in`, 
                                in.det        =   detGN$in.selected,
                                in.sel        =   selDT$in.selected,
                                in.up         =    upDT$in.selected,
                                in.down       =  downDT$in.selected,
                                in.up.genes   =    upDT$in.selected.genes,
                                in.down.genes =  downDT$in.selected.genes,
                                out           =   detGN$out,
                                det           =   detGN$selected,
                                sel           =   selDT$selected,
                                up            =    upDT$selected,
                                down          =  downDT$selected )
        enrichdt[, (      upvar) :=   upDT$p.selected]
        enrichdt[, (    downvar) := downDT$p.selected]
        enrichdt[, (selectedvar) :=  selDT$p.selected]
        enrichdt[, (   naivevar) :=  selGN$p.selected]
        enrichdt[, (detectedvar) :=  detGN$p.selected]
# Return
    enrichdt %<>% extract(in.det >= n)
    enrichdt %<>% extract(order(get(selectedvar)))
    if (!genes)  enrichdt %<>% extract(, .SD, .SDcols = !patterns('genes'))
    enrichdt[]
}

    collapse <- function(x, sep)  paste0(x, collapse = sep)


#' Count/Collapse in/outside intersection
#' @param x character OR list
#' @param y character 
#' @param ... used for S3 dispatch
#' @return number OR numeric
#' @examples
#' # Sets
#'    contrast1 <- c('a', 'b', 'c', 'd')
#'      pathway <- c('c', 'd', 'e', 'f')
#'    contrast2 <- c('e', 'f', 'g', 'h')
#'
#' # Count outside
#'    count_out(contrast1, pathway)
#'    count_out(list(contrast1 = contrast1, contrast2 = contrast2), pathway)
#'
#' # Count inside
#'    count_in(contrast1, pathway)
#'    count_in(list(contrast1 = contrast1, contrast2 = contrast2), pathway)
#'
#' # Collapse inside
#'    collapse_in(contrast1, pathway, sep = ' ')
#'    collapse_in(list(contrast1 = contrast1, contrast2 = contrast2), pathway, sep = ' ')
#' @export
count_in <- function(x, ...)  UseMethod('count_in', x)

#' @rdname count_in
#' @export
count_in.character <- function(x, y)  length(intersect(x, y))

#' @rdname count_in
#' @export
count_in.factor <- function(x, y)     length(intersect(x, y))

#' @rdname count_in
#' @export
count_in.list <- function(x, y)    lapply(x, count_in, y = y)
                            # dont vapply, list form required in data.table env

#' @rdname count_in
#' @export
collapse_in <- function(x, ...) UseMethod('collapse_in')

#' @rdname count_in
#' @export
collapse_in.character <- function(x, y, sep)  collapse(intersect(x, y), sep)

#' @rdname count_in
#' @export
collapse_in.factor    <- function(x, y, sep)  collapse(intersect(x, y), sep)

#' @rdname count_in
#' @export                         
collapse_in.list <- function(x, y, sep)   lapply(x, collapse_in, y = y, sep = sep)
                                   # dont vapply, list form required in data.table env

#' @rdname count_in
#' @export
count_out <- function(x, ...)  UseMethod('count_out', x)

#' @rdname count_in
#' @export
count_out.character <- function(x, y)  length(setdiff(x, y))

#' @rdname count_in
#' @export
count_out.factor    <- function(x, y)  length(setdiff(x, y))

#' @rdname count_in
#' @export
count_out.list <- function(x, y)    lapply(x, count_out, y = y)
                             # dont vapply, list form required in data.table env
