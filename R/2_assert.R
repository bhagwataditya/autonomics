
#==============================================================================
# has/contains


#' Does object have some svalues
#' @param object SummarizedExperiment
#' @param svar   sample var
#' @return logical
#' @examples
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
#' has_some_svalues(object, 'Group')
#' @noRd
has_some_svalues <- function(object, svar){
    if (is.null(svar))                          return(FALSE)
    if (!svar %in% autonomics::svars(object))   return(FALSE)
    if (all(is.na(svalues(object,svar)) | 
        svalues(object, svar)==''))             return(FALSE)
    return(TRUE)
}



#' Does object contain ratio values?
#' @param object SummarizedExperiment
#' @return logical
#' @examples
#' file <- download_data('billing16.proteingroups.txt')
#' object <- read_maxquant_proteingroups(file, plot = FALSE)
#' contains_ratios(object)
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot = FALSE)
#' contains_ratios(object)
#' @noRd
contains_ratios <- function(object){
    quantity <- metadata(object)$quantity
    if (is.null(quantity)) return(FALSE)
    return(stri_detect_fixed(quantity, 'Ratio'))
}


#==============================================================================
#
#                        assert_is_valid_sumexp
#
#==============================================================================


has_valid_fnames <- function(x, .xname = get_name_in_parent(x)){

    if (is.null(fnames(x))){
        return(false('fnames(%s) are NULL', .xname))}

    if (!all(fnames(x) == fdata(x)$feature_id)){
        return(false('fnames(%s) != fdata(%s)$feature_id', .xname, .xname))}

    #if (!all(fnames(x) == rownames(values(x)))){
    #    return(false('fnames(%s) != rownames(values(%s))', .xname, .xname))}

    #if (!all(fnames(x) == rownames(fdata(x)))){
    #    return(false('fnames(%s) != rownames(fdata(%s))', .xname, .xname))}

    TRUE
}


has_valid_snames <- function(x, .xname = get_name_in_parent(x)){

    if (is.null(snames(x))){
        return(false('snames(%s) are NULL', .xname))}

    if (!all(snames(x) == sdata(x)$sample_id)){
        return(false('snames(%s) != sdata(%s)$sample_id', .xname, .xname))}

    #if (!all(snames(x) == colnames(values(x)))){
    #    return(false('snames(%s) != colnames(values(%s))', .xname, .xname))}

    #if (!all(snames(x) == rownames(sdata(x)))){
    #    return(false('snames(%s) != colnames(sdata(%s))', .xname, .xname))}

    TRUE
}




#' Is valid SummarizedExperiment
#' @param x SummarizedExperiment
#' @param .xname see assertive.base::get_name_in_parent
#' @return TRUE or FALSE
#' @noRd
is_valid_sumexp <- function(x, .xname = get_name_in_parent(x)){
    if (!(ok <- assertive::is2(x, "SummarizedExperiment")))  return(ok)
    if (!(ok <- has_valid_fnames(x, .xname = .xname)))       return(ok)
    if (!(ok <- has_valid_snames(x, .xname = .xname)))       return(ok)
    TRUE
}


#' Assert that x is a valid SummarizedExperiment
#'
#' @param x SummarizedExperiment
#' @param .xname see assertive.base::get_name_in_parent
#' @return TRUE or FALSE
#' @examples
#' # VALID
#'     file <- download_data('halama18.metabolon.xlsx')
#'     x <- read_metabolon(file, plot = FALSE)
#'     assert_is_valid_sumexp(x)
#' # NOT VALID
#'     rownames(SummarizedExperiment::colData(x)) <- NULL
#'     # assert_is_valid_sumexp(x)
#' @export
assert_is_valid_sumexp <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_valid_sumexp, x, .xname = get_name_in_parent(x))
}

#' Is diann, fragpipe, proteingroups, phosphosites file?
#' @param x      file 
#' @param .xname name of x
#' @return NULL
#' @examples
#' file <- NULL
#' is_diann_file(file)
#' is_fragpipe_file(file)
#' is_proteingroups_file(file)
#' is_phosphosites(file)
#' 
#' file <- 3
#' is_diann_file(file)
#' is_fragpipe_file(file)
#' is_proteingroups_file(file)
#' is_phosphosites_file(file)
#' 
#' file <- 'blabla.tsv'
#' is_diann_file(file)
#' is_fragpipe_file(file)
#' is_proteingroups_file(file)
#' is_phosphosites_file(file)
#' 
#' file <- download_data('multiorganism.combined_protein.tsv')
#' is_diann_file(file)
#' is_fragpipe_file(file)
#' is_proteingroups_file(file)
#' is_phosphosites_file(file)
#' 
#' file <- download_data('dilution.report.tsv')
#' is_diann_file(file)
#' is_fragpipe_file(file)
#' is_proteingroups_file(file)
#' is_phosphosites_file(file)
#' 
#' file <- download_data('fukuda20.proteingroups.txt')
#' is_diann_file(file)
#' is_fragpipe_file(file)
#' is_proteingroups_file(file)
#' is_phosphosites_file(file)
#' 
#' file <- download_data('billing19.phosphosites.txt')
#' is_diann_file(file)
#' is_fragpipe_file(file)
#' is_proteingroups_file(file)
#' is_phosphosites_file(file)
#' @export
is_diann_file <- function(x, .xname = get_name_in_parent(x)){
    if (is.null(x)){                        false('%s is NULL',                  .xname)
    } else if (!is_a_string(x)){            false('%s is not a string',          .xname)
    } else if (!is_existing_file(x)){       false('%s does not exist',           .xname)
    } else if (col1(x) != 'File.Name'){     false('col1(%s) != "File.Name"',     .xname)
    } else if (col2(x) != 'Run'){           false('col2(%s) != "Run"',           .xname)
    } else if (col3(x) != 'Protein.Group'){ false('col3(%s) != "Protein.Group"', .xname)
    } else {                                TRUE 
    }
}

#' @rdname is_diann_file
#' @export
is_fragpipe_file <- function(x, .xname = get_name_in_parent(x)){
    if (is.null(x)){                      false('%s is NULL',                    .xname)
    } else if (!is_a_string(x)){          false('%s is not a string',            .xname)
    } else if (!is_existing_file(x)){     false('%s does not exist',             .xname)
    } else if (col1(x) != 'Protein'){     false('col1(%s) != "Protein"',         .xname)
    } else if (col2(x) != 'Protein ID'){  false('col2(%s) != "Protein ID"',      .xname)
    } else if (col3(x) != 'Entry Name'){  false('col3(%s) != "Entry Name"',      .xname)
    } else {                                TRUE 
    }
}

#' @rdname is_diann_file
#' @export
is_proteingroups_file <- function(x, .xname = get_name_in_parent(x)){
    if (is.null(x)){                                false('%s is NULL',                         .xname)
    } else if (!is_a_string(x)){                    false('%s is not a string',                 .xname)
    } else if (!is_existing_file(x)){               false('%s does not exist',                  .xname)
    } else if (col1(x) != 'Protein IDs'){           false('col1(%s) != "Protein IDs"',          .xname)
    } else if (col2(x) != 'Majority protein IDs'){  false('col2(%s) != "Majority protein ID"',  .xname)
    } else if (col3(x) != 'Peptide counts (all)'){  false('col3(%s) != "Peptide counts (all)"', .xname)
    } else {                                TRUE 
    }
}

#' @rdname is_diann_file
#' @export
is_phosphosites_file <- function(x, .xname = get_name_in_parent(x)){
    if (is.null(x)){                                    false('%s is NULL',                              .xname)
    } else if (!is_a_string(x)){                        false('%s is not a string',                      .xname)
    } else if (!is_existing_file(x)){                   false('%s does not exist',                       .xname)
    } else if (col1(x) != 'Proteins'){                  false('col1(%s) != "Proteins"',                  .xname)
    } else if (col2(x) != 'Positions within proteins'){ false('col2(%s) != "Positions within proteins"', .xname)
    } else if (col3(x) != 'Leading proteins'){          false('col3(%s) != "Leading proteins"',          .xname)
    } else {                                            TRUE 
    }
}

#' @rdname is_diann_file
#' @export
assert_diann_file <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_diann_file, x, .xname = .xname)
}

#' @rdname is_diann_file
#' @export
assert_fragpipe_file <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_fragpipe_file, x, .xname = .xname)
}

#' @rdname is_diann_file
#' @export
assert_proteingroups_file <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_proteingroups_file, x, .xname = .xname)
}

#' @rdname is_diann_file
#' @export
assert_phosphosites_file <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_phosphosites_file, x, .xname = .xname)
}


#--------

#' Is fastadt or NULL
#' @param x   fasta data.table
#' @param .xname string
#' @examples 
#' fastafile <- download_data('uniprot_hsa_20140515.fasta')    
#' x <- read_fastahdrs(fastafile)
#' is_fastadt_or_null(x)
#' @export
is_fastadt_or_null <- function(x, .xname = get_name_im_parent(x)){
    if (is.null(fastadt))                  return(TRUE)
    if (!is.data.table(fastadt))           return(false('%s is not a data.table', .xname))
    if (names(fastadt)[1] != 'reviewed')   return(false('col1(%s) != "reviewed"', .xname))
    if (names(fastadt)[2] != 'protein')    return(false('col2(%s) != "protein"',  .xname))
    return(TRUE)
}

#' @rdname is_fastadt_or_null
#' @export
assert_fastadt_or_null <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_fastadt_or_null, x, .xname = .xname)
}

