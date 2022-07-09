
#' Read diann proteingroups/precursors
#' 
#' @param file report.tsv file
#' @param quantity 'MaxLFQ' or 'Quantity'
#' @return data.table
#' @examples 
#' file <- download_data('szymanski22.report.tsv')
#' PR <- read_diann_precursors(file)     #    Precursor x Run (long)
#' PG <- read_diann_proteingroups(file)  # Proteingroup x Run (long)
#' PG[Quantity==top1] # matches      : 25 962 proteingroups
#' PG[Quantity!=top1] # doesnt match :  4 515 proteingroups
#' run <- 'IPT_HeLa_1_DIAstd_Slot1-40_1_9997'
#' PR[Protein.Group=='Q96JP5;Q96JP5-2' & Run == run][rev(order(Precursor.Quantity))] # 8884  ==   8884
#' PR[Protein.Group=='P36578'          & Run == run][rev(order(Precursor.Quantity))] # 8966  != 407978
#' @export
read_diann_precursors <- function(file){
    cols <- c('Genes', 'Protein.Group', 'Precursor.Id', 'Run', 
              'PG.MaxLFQ', 'Precursor.Quantity', 'PG.Quantity')
    fread(file, select = cols)
}

#' @rdname read_diann_precursors
#' @export
read_diann_proteingroups <- function(file){
    dt <- read_diann_precursors(file)
    dt[,        PG.Quantity := suppressWarnings(as.numeric(       PG.Quantity))]
    dt[, Precursor.Quantity := suppressWarnings(as.numeric(Precursor.Quantity))]
    dt[, .(MaxLFQ   = PG.MaxLFQ[1],
           Quantity = PG.Quantity[1], 
           top1     =     rev(sort(Precursor.Quantity))[1],
           top3     = sum(rev(sort(Precursor.Quantity))[1:3], na.rm = TRUE),
           Sum      = sum(         Precursor.Quantity,        na.rm = TRUE)), 
       by = c('Protein.Group', 'Run')]
}
