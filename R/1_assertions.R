
#==============================================================================
# has/contains


#' Does object have some svalues
#' @param object SummarizedExperiment
#' @param svar   sample var
#' @return logical
#' @examples
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' has_some_svalues(object, 'subgroup')
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
#' file <- system.file('extdata/billing19.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' contains_ratios(object)
#'
#' file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#' object <- read_metabolon(file)
#' contains_ratios(object)
#' @noRd
contains_ratios <- function(object)  any(grepl('[Rr]atio', assayNames(object)))


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
#' @param .xname see get_name_in_parent
#' @return TRUE or FALSE
#' @noRd
is_valid_sumexp <- function(x, .xname = get_name_in_parent(x)){
    if (!(ok <- is2(x, "SummarizedExperiment")))  return(ok)
    if (!(ok <- has_valid_fnames(x, .xname = .xname)))       return(ok)
    if (!(ok <- has_valid_snames(x, .xname = .xname)))       return(ok)
    TRUE
}


#' Assert that x is a valid SummarizedExperiment
#'
#' @param x SummarizedExperiment
#' @param .xname see get_name_in_parent
#' @return TRUE or FALSE
#' @examples
#' # VALID
#'     file <- system.file('extdata/atkin.metabolon.xlsx', package = 'autonomics')
#'     x <- read_metabolon(file)
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
#' is_diann_report(file)
#' is_fragpipe_tsv(file)
#' is_maxquant_proteingroups(file)
#' is_maxquant_phosphosites(file)
#'
#' file <- 3
#' is_diann_report(file)
#' is_fragpipe_tsv(file)
#' is_maxquant_proteingroups(file)
#' is_maxquant_phosphosites(file)
#'
#' file <- 'blabla.tsv'
#' is_diann_report(file)
#' is_fragpipe_tsv(file)
#' is_maxquant_proteingroups(file)
#' is_maxquant_phosphosites(file)
#'
#' file <- download_data('multiorganism.combined_protein.tsv')
#' is_diann_report(file)
#' is_fragpipe_tsv(file)
#' is_maxquant_proteingroups(file)
#' is_maxquant_phosphosites(file)
#'
#' file <- download_data('dilution.report.tsv')
#' is_diann_report(file)
#' is_fragpipe_tsv(file)
#' is_maxquant_proteingroups(file)
#' is_maxquant_phosphosites(file)
#'
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' is_diann_report(file)
#' is_fragpipe_tsv(file)
#' is_maxquant_proteingroups(file)
#' is_maxquant_phosphosites(file)
#'
#' file <- system.file('extdata/billing19.phosphosites.txt', package = 'autonomics')
#' is_diann_report(file)
#' is_fragpipe_tsv(file)
#' is_maxquant_proteingroups(file)
#' is_maxquant_phosphosites(file)
#' @export
is_diann_report <- function(x, .xname = get_name_in_parent(x)){
    if (is.null(x)){                        false('%s is NULL',                  .xname)
    } else if (!is_a_string(x)){            false('%s is not a string',          .xname)
    } else if (!is_existing_file(x)){       false('%s does not exist',           .xname)
    } else if (col1(x) != 'File.Name'){     false('col1(%s) != "File.Name"',     .xname)
    } else if (col2(x) != 'Run'){           false('col2(%s) != "Run"',           .xname)
    } else if (col3(x) != 'Protein.Group'){ false('col3(%s) != "Protein.Group"', .xname)
    } else {                                TRUE
    }
}

#' @rdname is_diann_report
#' @export
is_fragpipe_tsv <- function(x, .xname = get_name_in_parent(x)){
    if (is.null(x)){                      false('%s is NULL',                    .xname)
    } else if (!is_a_string(x)){          false('%s is not a string',            .xname)
    } else if (!is_existing_file(x)){     false('%s does not exist',             .xname)
    } else if (col1(x) != 'Protein'){     false('col1(%s) != "Protein"',         .xname)
    } else if (col2(x) != 'Protein ID'){  false('col2(%s) != "Protein ID"',      .xname)
    } else if (col3(x) != 'Entry Name'){  false('col3(%s) != "Entry Name"',      .xname)
    } else {                                TRUE
    }
}

#' @rdname is_diann_report
#' @export
is_maxquant_proteingroups <- function(x, .xname = get_name_in_parent(x)){
    if (is.null(x)){                                false('%s is NULL',                         .xname)
    } else if (!is_a_string(x)){                    false('%s is not a string',                 .xname)
    } else if (!is_existing_file(x)){               false('%s does not exist',                  .xname)
    } else if (col1(x) != 'Protein IDs'){           false('col1(%s) != "Protein IDs"',          .xname)
    } else if (col2(x) != 'Majority protein IDs'){  false('col2(%s) != "Majority protein ID"',  .xname)
    } else if (col3(x) != 'Peptide counts (all)'){  false('col3(%s) != "Peptide counts (all)"', .xname)
    } else {                                TRUE
    }
}

#' @rdname is_diann_report
#' @export
is_maxquant_phosphosites <- function(x, .xname = get_name_in_parent(x)){
    if (is.null(x))                               return(false('%s is NULL',                              .xname))
    if (!is_a_string(x))                          return(false('%s is not a string',                      .xname))
    if (!is_existing_file(x))                     return(false('%s does not exist',                       .xname))
    if (col1(x) != 'Proteins')                    return(false('col1(%s) != "Proteins"',                  .xname))
    if (col2(x) != 'Positions within proteins')   return(false('col2(%s) != "Positions within proteins"', .xname))
    if (col3(x) != 'Leading proteins')            return(false('col3(%s) != "Leading proteins"',          .xname))
    return(TRUE)
}

#' @rdname is_diann_report
#' @export
is_compounddiscoverer_output <- function(x, .xname = get_name_in_parent(x)){
  if (is.null(x))                  { false('%s is NULL',         .xname)
  } else if (!is_a_string(x))      { false('%s is not a string', .xname)
  } else if (!is_existing_file(x)) { false('%s does not exist',  .xname)
  } else if (col1(x) != 'Name')    { false('col1(%s) != "Name"', .xname)
  } else                           { TRUE
  }
}

#' @rdname is_diann_report
#' @export
assert_diann_report <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_diann_report, x, .xname = .xname)
}

#' @rdname is_diann_report
#' @export
assert_fragpipe_tsv <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_fragpipe_tsv, x, .xname = .xname)
}

#' @rdname is_diann_report
#' @export
assert_maxquant_proteingroups <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_maxquant_proteingroups, x, .xname = .xname)
}

#' @rdname is_diann_report
#' @export
assert_maxquant_phosphosites <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_maxquant_phosphosites, x, .xname = .xname)
}

#' @rdname is_diann_report
#' @export
assert_compounddiscoverer_output <- function(x, .xname = get_name_in_parent(x)){
  assert_engine(is_compounddiscoverer_output, x, .xname = .xname)
}

#--------

#' Is fastadt
#' @param x   fasta data.table
#' @param .xname string
#' @examples
#' fastafile <- system.file('extdata/uniprot_hsa_20140515.fasta', package = 'autonomics')
#' x <- read_uniprotdt(fastafile)
#' # is_fastadt(x)  # slow
#' @export
is_fastadt <- function(x, .xname = get_name_in_parent(x)){
    if (!is.data.table(x))       return(false('%s is not a data.table', .xname))
    if (names(x)[1] != 'dbid')   return(false('col1(%s) != "uniprot"', .xname))
    return(TRUE)
}

#' @rdname is_fastadt
#' @export
assert_fastadt <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_fastadt, x, .xname = .xname)
}


#---------------

#' Is scalar subset
#' @param x scalar
#' @param y SummarizedExperiment
#' @param .xname name of x
#' @param .yname name of y
#' @examples
#' file <- system.file('extdata/fukuda20.proteingroups.txt', package = 'autonomics')
#' object <- read_maxquant_proteingroups(file)
#' is_scalar_subset('subgroup',     svars(object))
#' is_scalar_subset('subject',      svars(object))
#' assert_scalar_subset('subgroup', svars(object))
#' @export
is_scalar_subset <- function(x, y, .xname = get_name_in_parent(x), .yname = get_name_in_parent(y)){
    if (!(ok <- is_scalar(x, .xname = .xname)))                       return(ok)
    if (!(ok <- is_subset(x, y, .xname = .xname, .yname = .yname))){
        return(false("%s is not in %s", .xname, .yname))
    }
    return(TRUE)
}

#' @rdname is_scalar_subset
#' @export
assert_scalar_subset <- function(x, y, .xname = get_name_in_parent(x), .yname = get_name_in_parent(y)){
    assert_engine(is_scalar_subset, x, y, .xname = .xname, .yname = .yname)
}

#-------------------

#' Is positive number
#' @param x number
#' @param .xname name of x
#' @return TRUE or false
#' @examples
#' is_positive_number( 3)
#' is_positive_number(-3)
#' is_positive_number( 0)
#' is_weakly_positive_number(0)
#' assert_positive_number(3)
#' @export
is_positive_number <- function(x, .xname = get_name_in_parent(x)){
    if (!is_a_number(x, .xname = .xname))                     return(false('%s is not a number',  .xname))
    if (!is_greater_than(x, 0,.xname = .xname))               return(false('%s <= 0', .xname))
    return(TRUE)
}

#' @rdname is_positive_number
#' @export
assert_positive_number <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_positive_number, x, .xname = .xname)
}

#' @rdname is_positive_number
#' @export
is_weakly_positive_number <- function(x, .xname = get_name_in_parent(x)){
    if (!is_a_number(x, .xname = .xname))                     return(false('%s is not a number',  .xname))
    if (!is_greater_than_or_equal_to(x, 0,.xname = .xname))   return(false('%s < 0', .xname))
    return(TRUE)
}


#' @rdname is_positive_number
#' @export
assert_weakly_positive_number <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_weakly_positive_number, x, .xname = .xname)
}

#---------------------

#' Is fraction
#' @param  x      number
#' @param .xname  string
#' @return TRUE or false
#' @examples
#' is_fraction(0.1)          # YES
#' is_fraction(1)            # YES
#' is_fraction(1.2)          # NO - more than 1
#' is_fraction(c(0.1, 0.2))  # NO - vector
#' @export
is_fraction <- function(x, .xname = get_name_in_parent(x)){
    if (!(ok <- is_a_number(x, .xname = .xname)))                       return(ok)
    if (!is_in_closed_range(x, lower = 0, upper = 1, .xname = .xname))  return(false('%s is not a fraction', .xname))
    TRUE
}

#' @rdname is_fraction
#' @export
assert_is_fraction <- function(x, .xname = get_name_in_parent(x)){
    assert_engine(is_fraction, x, .xname = .xname)
}


#---------------------


#' Variable has multiple levels?
#' @param  x      vector, data.table or SummarizedExperiment
#' @param .xname  string
#' @param  y      string
#' @param .yname  string
#' @param  ...    required for s3 dispatch
#' @return  TRUE or false
#' @examples
#' # numeric
#'     a <- numeric();                               has_multiple_levels(a)
#'     a <- c(1, 1);                                 has_multiple_levels(a)
#'     a <- c(1, 2);                                 has_multiple_levels(a)
#' # character
#'     a <- character();                             has_multiple_levels(a)
#'     a <- c('A', 'A');                             has_multiple_levels(a)
#'     a <- c('A', 'B');                             has_multiple_levels(a)
#' # factor
#'     a <- factor();                                has_multiple_levels(a)
#'     a <- factor(c('A', 'A'));                     has_multiple_levels(a)
#'     a <- factor(c('A', 'B'));                     has_multiple_levels(a)
#' # data.table
#'     dt <- data.table(a = factor());               has_multiple_levels(dt, 'b')
#'     dt <- data.table(a = factor());               has_multiple_levels(dt, 'a')
#'     dt <- data.table(a = factor());               has_multiple_levels(dt, 'a')
#'     dt <- data.table(a = factor(c('A', 'A')));    has_multiple_levels(dt, 'a')
#'     dt <- data.table(a = factor(c('A', 'B')));    has_multiple_levels(dt, 'a')
#' # sumexp
#'     object <- matrix(1:9, nrow = 3)
#'     rownames(object) <- sprintf('f%d', 1:3)
#'     colnames(object) <- sprintf('s%d', 1:3)
#'     object <- list(exprs = object)
#'     object %<>% SummarizedExperiment::SummarizedExperiment()
#'     object$subgroup <- c('A', 'A', 'A');          has_multiple_levels(object, 'group')
#'     object$subgroup <- c('A', 'A', 'A');          has_multiple_levels(object, 'subgroup')
#'     object$subgroup <- c('A', 'B', 'A');          has_multiple_levels(object, 'subgroup')
#' @export     
has_multiple_levels <- function(x, ...)  UseMethod('has_multiple_levels')


#' @rdname has_multiple_levels
#' @export
has_multiple_levels.character <- function(x, .xname = get_name_in_parent(x), ...){
    n <- length(unique(x))
    if (! n > 1)  return(false('%s has only %d level(s)', .xname, n))
    TRUE
}


#' @rdname has_multiple_levels
#' @export
has_multiple_levels.factor <- function(x, .xname = get_name_in_parent(x), ...){
    has_multiple_levels.character(x = x, .xname = .xname)
}


#' @rdname has_multiple_levels
#' @export
has_multiple_levels.numeric <- function(x, .xname = get_name_in_parent(x), ...){
    has_multiple_levels.character(x = x, .xname = .xname)
}


#' @rdname has_multiple_levels
#' @export
has_multiple_levels.data.table <- function(
    x,   # data.table
    y,   # var
    .xname = get_name_in_parent(x),
    .yname = get_name_in_parent(y), ...
){
    if (!(ok <- is_scalar_subset(y, names(x), .xname = .yname, .yname = .xname)))  return(ok)
    if (!(ok <- has_multiple_levels.factor(  x[[y]], .xname = .yname)))            return(ok)
    TRUE
}


#' @rdname has_multiple_levels
#' @export
has_multiple_levels.SummarizedExperiment <- function(
     x,  # sumexp
     y,  # svar
    .xname = get_name_in_parent(x),
    .yname = get_name_in_parent(y), ...
){
    if(!(ok <- has_multiple_levels.data.table(
                    sdt(x), y, .xname = .xname, .yname = .yname)))  return(ok)
    TRUE
}


#---------------------


#' Is valid formula
#' @param x formula
#' @param y SummarizedExperiment
#' @param .xname string
#' @param .yname string
#' @return TRUE or false
#' @examples 
#' object <- matrix(1:9, nrow = 3)
#' rownames(object) <- sprintf('f%d', 1:3)
#' colnames(object) <- sprintf('s%d', 1:3)
#' object <- list(exprs = object)
#' object %<>% SummarizedExperiment::SummarizedExperiment()
#' object$group    <- 'group0'
#' object$subgroup <- c('A', 'B', 'C')
#' svars(object)
#'     is_valid_formula( 'condition',   object)   # not formula
#'     is_valid_formula( ~condition,    object)   # not svar
#'     is_valid_formula( ~group,        object)   # not multilevel
#'     is_valid_formula( ~subgroup,     object)   # TRUE
#'     is_valid_formula( ~0+subgroup,   object)   # TRUE
#'     is_valid_formula( ~1,            object)   # TRUE
#' assert_valid_formula( ~subgroup,     object)
#' @export
is_valid_formula <- function(
    x,  # formula
    y,  # object
    .xname = get_name_in_parent(x), 
    .yname = get_name_in_parent(y)
){
    if (!(ok <- is_one_sided_formula(x, .xname = .xname)))                             return(ok)
    if (!(ok <- is_subset(all.vars(x), svars(y), .xname = .xname, .yname = .yname)))   return(ok)
    for (var in all.vars(x)){
        if (!(ok <- has_multiple_levels(y, var, .xname = .yname, .yname = .xname)))   return(ok)
    }
    TRUE
}

#' @rdname is_valid_formula
#' @export
assert_valid_formula <- function(
    x, y, .xname = get_name_in_parent(x), .yname = get_name_in_parent(y)
){
    assert_engine(is_valid_formula, x = x, y = y, .xname = .xname, .yname = .yname)
}

#---------------------

all_have_setidentical_colnames <- function(x, .xname = get_name_in_parent(x))
{
    assert_is_list(x)
    assert_all_are_true(sapply(x, is_data.frame))
    aicns <- lapply(x, names) %>%
      lapply(sort) %>%
      sapply(identical, .[[1]]) %>%
      all()
    if (! aicns)  return(false('Not all colnames in %s are setidentical', .xname))
    TRUE
}

assert_all_have_setidentical_colnames <- function(
    x, .xname = get_name_in_parent(x)
){
  assert_engine(all_have_setidentical_colnames, x = x, .xname = .xname)
}
