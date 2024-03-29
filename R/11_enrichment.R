

#' Human/Mouse Msigdb Collections
#' @export
MSIGCOLLECTIONSHUMAN <- c( 
    'C1',                   # positional    : - genes in same cytogenetic band 
    'C2:CGP',               # curated       : + genes with altered expression upon (chemical/genetic) perturbation
    'C2:CP',                #
    'C2:CP:BIOCARTA',       #                 + genes in same biocarta pathway
    'C2:CP:KEGG_LEGACY',    #                 -                   kegg pathway (legacy)
    'C2:CP:KEGG_MEDICUS',   #                 +           kegg medicus pathway
    'C2:CP:PID',            #                 -        prot interac db pathway
    'C2:CP:REACTOME',       #                 +               reactome pathway
    'C2:CP:WIKIPATHWAYS',   #                 +                   wiki pathway
    'C3:MIR:MIR_LEGACY',    # regulatory    : + genes with shared mir site (mirdb 6.0, 2020)
    'C3:MIR:MIRDB',         #                 -                            (legacy)
    'C3:TFT:GTRD',          #                 +                   tfb site (Gene Transcr Reg Db, 2021)
    'C3:TFT:TFT_LEGACY',    #                 -                            (legacy, 2005)
    'C4:3CA',               # computational : ? genes upregulated in tumor subpopulation  (2023)
    'C4:CGN',               #                 ? genes in same cancer driver neighbourhood (2005)
    'C4:CM',                #                 ? genes with altered expression in (same) cancer condition (2004)
    'C5:GO:BP',             # ontological   : + genes in same biological process
    'C5:GO:CC',             #                 +             cellular compartment
    'C5:GO:MF',             #                 +     with same molecular function
    'C5:HPO',               #                 + with altered expression in (same) human phenotype
    'C6',                   # oncogenic     : ? genes in cancer-dysregulated pathway
    'C7:IMMUNESIGDB',       # immunologic   : ? genes with altered expression upon (same) immune system perturbation
    'C7:VAX',               #                 ? genes with altered expression upon (same) vaccination
    'C8',                   # celltype      : + genes with celltype-specific expression
    'H'                     # hallmark      : + shared hallmark
 )


#' @rdname MSIGCOLLECTIONSHUMAN
#' @export
MSIGCOLLECTIONSMOUSE <- c( 
    'M1',                   # positional    : - genes in same cytogenetic band 
    'M2:CGP',               # curated       : + genes with altered expression upon (chemical/genetic) perturbation
    'M2:CP:BIOCARTA',       #                 + genes in same biocarta pathway
    'M2:CP:REACTOME',       #                 +               reactome pathway
    'M2:CP:WIKIPATHWAYS',   #                 +                   wiki pathway
    'M3:GTRD',              #                 +                   tfb site (Gene Transcr Reg Db, 2021)
    'M3:MIRDB',             #                 -                            (legacy)
    'M5:GO:BP',             # ontological   : + genes in same biological process
    'M5:GO:CC',             #                 +             cellular compartment
    'M5:GO:MF',             #                 +     with same molecular function
    'M5:MPT',               #                 + tumor phenotype ontology 
    'C8',                   # celltype      : + genes with celltype-specific expression
    'MH'                    # hallmark      : + shared hallmark
)


#' local msigdb dir
#' @export
MSIGDIR <- file.path(R_user_dir('autonomics', 'cache'), 'msigdb')


#' @title       list files
#' @description list.files for programming
#' @details
#' Adds a small layer on list.files.
#' Returning NULL rather than character(0) when no files.
#' Making it better suited for programming.
#' @param dir         directory
#' @param full.names  TRUE or FALSE
#' @export
list_files <- function(dir, full.names){
    y <- list.files(dir, full.names = full.names)
    if (length(y)==0)  return(NULL)  else  return(y)
}


#' Read msigdb datatable
#' @param file         msigdb file: one of the files in dir(MSIGDB).
#' @param collections  subset of names(MSIGCOLLECTIONS)
#' @examples
#' read_msigdt()
#' @export
read_msigdt <- function(
           file = list_files(MSIGDIR, full.names = TRUE)[1], 
    collections = if (is.null(file))  NULL else 
                  switch( basename(file) %>% substr(nchar(.)-4, nchar(.)-3) , 
                          Hs = c( 'C2:CP:REACTOME', 'C5:GO:BP', 'C5:GO:MF', 'C5:GO:CC' ), 
                          Mm = c( 'M2:CP:REACTOME', 'M5:GO:BP', 'M5:GO:MF', 'M5:GO:CC' ))
){
# Assert
    if (is.null(file)){
        cmessage("\t\tVisit https://www.gsea-msigdb.org/gsea/downloads.jsp")
        cmessage("\t\tScrolldown. Locate SQLite database (not jason or xml!)")
        cmessage("\t\tDownload Human or Mouse SQLite database")
        cmessage("\t\tCreate %s", MSIGDIR)
        cmessage("\t\tUnzip into this dir")  
        cmessage("\t\tNow rerun `read_msigdt()`")
        return(NULL)
    }
    if (!requireNamespace('DBI',     quietly = TRUE))  message("BiocManager::install('DBI'). Then re-run.")
    if (!requireNamespace('RSQLite', quietly = TRUE))  message("BiocManager::install('RSQLite'). Then re-run.")
    assert_all_are_existing_files(file)
    assert_is_subset(collections, c(MSIGCOLLECTIONSHUMAN, MSIGCOLLECTIONSMOUSE))
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
# Merge
    msigdt <- dt1
    msigdt %<>% merge(dt2, by = 'geneid')
    msigdt %<>% merge(dt3, by = 'setid')
    msigdt %<>% extract(, c('collection', 'set', 'gene'), with = FALSE)
# Return
    msigdt %<>% extract(collection %in% collections)
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

#' guess fitsep
#' @param featuredt data.table
#' @return string
#' @examples 
#' file <- download_data('fukuda20.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file)
#' guess_fitsep(fdt(object))
#' guess_fitsep(fdt(fit_limma(object)))
#' @export
guess_fitsep <- function(featuredt){
    idx <- names(featuredt) %>% stri_detect_regex('^effect')
    if (all(!idx))  return(NULL)
    sep <- names(featuredt)[idx][1] %>% substr(7,7)
    return(sep)
}


#' Abstract model fit
#' @param object            SummarizedExperiment
#' @param sep               string
#' @param fit               character vector
#' @param coef              character vector
#' @param significancevar  'p' or 'fdr'
#' @param significance      fraction : pvalue cutoff
#' @return SummarizedExperiment
#' @examples
#' file <- download_data('atkin.metabolon.xlsx')
#' object <- read_metabolon(file, fit = 'limma', coef = 't3')
#' fdt(object)
#' fdt(abstract_fit(object))
#' @export
abstract_fit <- function(
             object, 
                sep = guess_fitsep(fdt(object)),
                fit = fits(fdt(object)), 
               coef = coefs(fdt(object), fit = fit), 
    significancevar = 'p', 
       significance = 0.05
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_subset(fit,   fits(fdt(object)))
    assert_is_subset(coef, coefs(fdt(object)))
# Abstract
    for ( curfit in fit){
    for (curcoef in coef){
        abstractvar <- paste(curcoef, curfit, sep = sep)
            pvalues <- modelvec(fdt(object), 'p',      fit = curfit, coef = curcoef)
       effectvalues <- modelvec(fdt(object), 'effect', fit = curfit, coef = curcoef)
        fdt(object)[[ abstractvar ]] <- 'flat'
        fdt(object)[[ abstractvar ]][ pvalues<significance  &  effectvalues<0 ] <- 'down' 
        fdt(object)[[ abstractvar ]][ pvalues<significance  &  effectvalues>0 ] <- 'up' 
        fdt(object)[[ abstractvar ]] %<>% factor(c('flat', 'up', 'down'))
    }}
    object
}


#' group by level
#' @param x named logical/character/factor
#' @param var    string
#' @param idvar  string
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
group_by_level <- function(x, ...)  UseMethod('group_by_level')


#' @rdname group_by_level
#' @export
group_by_level.character <- function(x, ...){
    fun <- function(y) names(x)[x == y]
    Map(fun, unique(x))                          
}

#' @rdname group_by_level
#' @export
group_by_level.factor <- function(x, ...){
    y <- x
    y %<>% as.character() 
    y %<>% set_names(names(x))
    group_by_level.character(y)
}

#' @rdname group_by_level
#' @export
group_by_level.data.table <- function(x, var, idvar, ...){
    y <- x[[var]]
    names(y) <- x[[idvar]]
    group_by_level(y)
}

#' logical to factor
#' @param x     logical vector
#' @param true  string : truelevel
#' @param false string : falselevel
#' @return factor
#' @examples
#' t1up <- c( TRUE,   FALSE,  TRUE)
#' t1   <- c('flat', 'down', 'up' )  %>%  factor(., .)
#' t1up
#' logical2factor(t1up)
#' factor2logical(t1)
#' @export
logical2factor <- function(
    x, 
    true = get_name_in_parent(x),
    false = paste0('not', true)
){
# Assert
    assert_is_logical(x)
# Convert
    y <- rep(false, length(x))
    y[x] <- true
    y%<>% factor(c(true, false)) # level1 = (noninteresting) baseline
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
#' @param fit       string
#' @param coef      string
#' @param var       selection fvar
#' @param levels    selection levels
#' @param genevar   gene fvar
#' @param genesep   gene separator (string)
#' @param n         number
#' @param verbose   whether to msg
#' @param genes     whether to report genes
#' @examples
#' # Read
#'     pathwaydt <- read_msigdt(collections = 'C5:GO:BP')
#'     file <- download_data('atkin.somascan.adat')
#'     object <- read_somascan(file, fit = 'limma', coefs = 't1')
#'     fvars(object) %<>% gsub('EntrezGeneSymbol', 'gene', .)
#'     object %<>% abstract_fit()
#'     var <- abstractvar(fdt(object))
#'     enrichdt1 <- enrichment(object, pathwaydt, var = var)                                   # 2:n factor 
#'     enrichdt2 <- enrichment(object, pathwaydt, var = var, levels = c('flat', 'down', 'up')) # 1:n factor
#'     enrichdt3 <-  altenrich(object, pathwaydt)                                              # alternative implementation
#'     cols <- intersect(names(enrichdt1), names(enrichdt3))
#'     all(enrichdt1[, cols, with = FALSE]  ==  enrichdt3[, cols, with = FALSE])   # identical
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
          fit = fits(fdt(object))[1],
         coef = coefs(fdt(object), fit = fit)[1],
          var = abstractvar(fdt(object), fit = fit, coef = coef),
       levels = fdt(object)[[var]] %>% base::levels() %>% extract(-1),
      genevar = 'gene', 
      genesep = '[ ,;]',
            n = 3,
      verbose = TRUE,
        genes = FALSE
){
# Assert
    assert_is_valid_sumexp(object)
    assert_is_data.table(pathwaydt)
    if (!is.null(fit ))  assert_scalar_subset(fit,  fits(fdt(object)))
    if (!is.null(coef))  assert_scalar_subset(coef, coefs(fdt(object)))
    # assert_is_vector(levels(fdt(object)[[var]]))
    assert_scalar_subset(var, fvars(object), .xname = get_name_in_parent(var))
    assert_is_factor(fdt(object)[[var]])
    assert_is_subset(levels, levels(fdt(object)[[var]]))
    assert_scalar_subset( genevar, names(pathwaydt))
    if (any(is_missing_or_empty_character(fdt(object)[[genevar]]))){
        cmessage("\t\tFirst run: `object %%<>%% filter_features(%s!='')`", genevar)
        cmessage("\t\t Then run: `enrichdt <- enrichment(object)`")
        return(NULL)
    }
    assert_all_are_non_missing_nor_empty_character(  pathwaydt[[genevar]])
    if (!is.null(genesep))  assert_is_a_string(genesep)
    assert_is_a_bool(verbose)
    if (verbose)  cmessage('\tAre pathways enriched in `%s` genes ?', var)
    if (verbose)  cmessage("\t\tfdt(.)$%s  %%<>%%  split_extract_regex('%s', 1)", genevar, genesep)
    fdt(object)[[genevar]] %<>% split_extract_regex(genesep, 1)
    in.det <- in.sel <- det <- sel <- out <- gene <- NULL
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
#' @param coef       \code{string} in \code{coefs(fdt(object))}
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
            genevar = 'gene', 
            genesep = '[ ,;]',
               coef = default_coefs(fdt(object))[1], 
                fit = fits(fdt(object))[1],
    significancevar = 'p',
       significance = 0.05,
         effectsize = 0,
                  n = 3, 
              genes = FALSE,
            verbose = TRUE
){
# Assert
    assert_is_valid_sumexp(object)
    if (is.null(pathwaydt)) return(NULL)
    assert_is_data.table(   pathwaydt)
    assert_scalar_subset(genevar, fvars(object))
    assert_scalar_subset(genevar, names(pathwaydt ))
    assert_all_are_non_missing_nor_empty_character(fdt(object)[[genevar]])
    assert_all_are_non_missing_nor_empty_character(  pathwaydt[[genevar]])
    if (!is.null(genesep))  assert_is_a_string(genesep)
    assert_scalar_subset(coef, coefs(fdt(object)))
    assert_scalar_subset(fit, fits(fdt(object)))
    genes0 <- fdt(object)[[genevar]]
    fdt(object)[[genevar]] %<>% split_extract_regex(genesep, 1)
    gene <- in.up <- in.det <- in.down <- in.sel <- `in` <- out <- NULL
# Function
    object %<>% abstract_fit(fit = fit, coef = coef, significancevar = significancevar)
    abstractvar <- abstractvar(fdt(object), fit = fit, coef = coef)
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
#' @param x   character OR list
#' @param y   character 
#' @param sep string
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
count_in.character <- function(x, y, ...)  length(intersect(x, y))

#' @rdname count_in
#' @export
count_in.factor <- function(x, y, ...)     length(intersect(x, y))

#' @rdname count_in
#' @export
count_in.list <- function(x, y, ...)    lapply(x, count_in, y = y)
                            # dont vapply, list form required in data.table env

#' @rdname count_in
#' @export
collapse_in <- function(x, ...) UseMethod('collapse_in')

#' @rdname count_in
#' @export
collapse_in.character <- function(x, y, sep, ...)  collapse(intersect(x, y), sep)

#' @rdname count_in
#' @export
collapse_in.factor    <- function(x, y, sep, ...)  collapse(intersect(x, y), sep)

#' @rdname count_in
#' @export                         
collapse_in.list <- function(x, y, sep, ...)   lapply(x, collapse_in, y = y, sep = sep)
                                   # dont vapply, list form required in data.table env

#' @rdname count_in
#' @export
count_out <- function(x, ...)  UseMethod('count_out', x)

#' @rdname count_in
#' @export
count_out.character <- function(x, y, ...)  length(setdiff(x, y))

#' @rdname count_in
#' @export
count_out.factor    <- function(x, y, ...)  length(setdiff(x, y))

#' @rdname count_in
#' @export
count_out.list <- function(x, y, ...)    lapply(x, count_out, y = y)
                             # dont vapply, list form required in data.table env
