# MaxLFQ: diann-r versus diann-cpp
#---------------------------------
    require(magrittr)
    require(data.table)
    
    # Read
    url <- 'https://bitbucket.org/graumannlabtools/autonomics/downloads/dilution.report.tsv'
    file <- file.path(tempdir(), basename(url))
    # download.file(url, destfile = file, mode = 'wb')
    dt <- fread(file)
    dt %<>% extract(unique(Protein.Names)[1:2], on = 'Protein.Names')
    dt$File.Name %<>% factor()
    levels(dt$File.Name) %<>% substr(nchar(.)-2, nchar(.)-2)
    levels(dt$File.Name) %<>% paste0('_', .)
    dt$File.Name %<>% as.character()
    
    # cpp MaxLFQ
    cmat <- unique(dt[, .(Protein.Names, File.Name, PG.MaxLFQ)])
    cmat %<>% dcast.data.table(Protein.Names ~ File.Name, value.var = 'PG.MaxLFQ')
    
    # r MaxLFQ
    x <- diann::diann_maxlfq(dt)
    x %<>% extract(, names(cmat)[-1])
    x %<>% extract(cmat$Protein.Names, )
    
    # Seem to not match
    cmat
    x

# diann: Top-1
#-------------
    require(magrittr)
    require(data.table)
    
    # Read
    url <- 'https://bitbucket.org/graumannlabtools/autonomics/downloads/dilution.report.tsv'
    file <- file.path(tempdir(), basename(url))
    # download.file(url, destfile = file, mode = 'wb')
    dt <- fread(file)
    dt$File.Name %<>% factor()
    levels(dt$File.Name) %<>% substr(nchar(.)-2, nchar(.)-2)
    levels(dt$File.Name) %<>% paste0('_', .)
    dt$File.Name %<>% as.character()
    
    dt %<>% extract(, c('Protein.Names', 'File.Name', 'Precursor.Id', 'Precursor.Normalised', 'Precursor.Quantity', 'PG.Quantity'), with = FALSE)
    dt %<>% extract('EIF3J_HUMAN', on = 'Protein.Names')
    dt %<>% extract(c('_0', '_3', '_6', '_9'), on = 'File.Name')
    
    dt$Precursor.Quantity %<>% as.numeric()
    dt$Precursor.Normalised %<>% as.numeric()
    dt$PG.Quantity %<>% as.numeric()

    dt %<>% extract(order(Protein.Names, File.Name, -Precursor.Quantity))
    dt[, PG.Top1 := max(Precursor.Quantity), by = c('Protein.Names', 'File.Name')]
    dt    
    
