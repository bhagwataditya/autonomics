#=======================================================================
#
#                      assertive.properties
#
#=======================================================================


#----------------------------
#
#       assert_has_names
#
#-----------------------------

#' Does the input have names?
#'
#' Checks to see if the input has (row/column/dimension) names.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{has_names} returns \code{TRUE} if \code{names} is 
#' non-null. 
#' \code{has_rownames}, \code{has_colnames} and \code{has_dimnames} work
#' in a similar fashion, checking the corresponding attributes.
#' \code{assert_has_names} returns nothing but throws an error if 
#' \code{has_names} is not \code{TRUE}.
#' @note Empty names (i.e., \code{""}) are not allowed in R, and are 
#' not checked here.
#' @seealso \code{\link[base]{names}}, \code{\link[base]{rownames}}, 
#' \code{\link[base]{colnames}}, \code{\link[base]{dimnames}}.
#' @examples
#' assert_has_names(c(a = 1, 2))
#' dfr <- data.frame(x = 1:5)
#' assert_has_rownames(dfr)
#' assert_has_colnames(dfr)
#' assert_has_dimnames(dfr)
#' @author Richard Cotton
#' @noRd
has_names <- function(x, .xname = get_name_in_parent(x))
{
  namesx <- names(x)
  if(is.null(namesx)) 
  {
    return(false("The names of %s are NULL.", .xname))
  }
  if(!any(nzchar(namesx))) 
  {
    return(false("The names of %s are all empty.", .xname))
  }
  TRUE
} 

assert_has_names <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                            
  assert_engine(
    has_names, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


#-----------------------------------
#
#      assert_has_no_duplicates
#
#-----------------------------------

#' Does the input have duplicates?
#'
#' Checks to see if the input has duplicates.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{has_duplicates} returns \code{TRUE} if\code{anyDuplicated} 
#' is \code{TRUE}.  \code{assert_has_duplicates} returns nothing but 
#' throws an error if \code{has_duplicates} is not \code{TRUE}. 
#' \code{has_no_duplicates} is the negation of \code{has_duplicates}.
#' @seealso \code{\link{anyDuplicated}}.
#' @examples 
#' x <- sample(10, 100, replace = TRUE)
#' assert_has_duplicates(x)
#' has_no_duplicates(x)
#' @author Richard Cotton
#' @noRd
has_duplicates <- function(x, .xname = get_name_in_parent(x))
{
  if(!anyDuplicated(x)) 
  {
    return(false(gettext("%s has no duplicates."), .xname))
  }
  TRUE
}

has_no_duplicates <- function(x, .xname = get_name_in_parent(x))
{
  if(anyDuplicated(x)) 
  {
    dupe_indicies <- which(duplicated(x))
    return(
      false(
        ngettext(
          length(dupe_indicies),
          "%s has a duplicate at position %s.",
          "%s has duplicates at positions %s."
        ), 
        .xname, 
        toString(dupe_indicies, width = 100)
      )
    )
  }
  TRUE
}

assert_has_no_duplicates <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                             
  assert_engine(
    has_no_duplicates, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


#----------------------------
#
#      assert_is_null
#      assert_is_not_null
#
#----------------------------

#' Is the input (not) null?
#'
#' Checks to see if the input is (not) null.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_null} wraps \code{is.null}, providing more 
#' information on failure. \code{is_not_null} returns \code{TRUE} in
#' the opposite case.  The \code{assert_*} functions return nothing but
#' throw an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @seealso \code{\link[base]{is.null}}.
#' @examples
#' # Predicate for NULL. 
#' is_null(NULL)
#' is_null(c())
#' 
#' # Atomic vectors of length zero are not NULL!
#' is_null(numeric())
#' # ... and neither is NA
#' is_null(NA)
#' 
#' # The opposite check
#' is_not_null(NULL)
#' is_not_null(c())
#' is_not_null(numeric())
#' 
#' # These checks should pass
#' assert_is_null(NULL)
#' assert_is_null(c())
#' assert_is_not_null(NA)
#' 
#' # This should fail
#' assertive.base::dont_stop(assert_is_null(NaN))
#' @author Richard Cotton
#' @noRd
is_null <- function(x, .xname = get_name_in_parent(x))
{
  if(!is.null(x))
  {
    return(false("%s is not NULL; its value is %s.", .xname, assertive.base::safe_deparse(x)))
  }
  TRUE
}

is_not_null <- function(x, .xname = get_name_in_parent(x))
{
  if(is.null(x))
  {
    return(false("%s is NULL.", .xname))
  }
  TRUE
}

assert_is_not_null <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                      
  assert_engine(is_not_null, x, .xname = get_name_in_parent(x))   
}


#---------------------
#
#    assert_is_empty
#    assert_is_scalar
#
#---------------------

is_of_length <- function(x, n, .xname = get_name_in_parent(x))
{
  n <- assertive.base::use_first(n)
  check_n(n)
  length_x <- length(x)
  if(length_x != n)
  {
    return(false("%s has length %d, not %d.", .xname, length_x, n))
  }
  TRUE
}

#' Get the dimensions of an object
#' 
#' Get the dimensions of an object, retuning the length if that object has no
#' \code{dim} attribute.
#' @param x Any object.
#' @return A integer vector of non-negative values.
#' @seealso \code{\link[base]{NROW}}, \code{\link[base]{dim}}
#' @examples
#' # For data frames and matrices, DIM is the same as dim.
#' DIM(sleep) 
#' # For vectors (and other objects without a dim attribute), DIM is the 
#' # same as length.
#' DIM(1:10)
#' DIM(list(x = 1:10))
#' @noRd
DIM <- function(x)
{
  dim_x <- dim(x)
  if(is.null(dim_x)) length(x) else dim_x
}

n_elements <- function(x)
{
  if(is.recursive(x))
  {
    sum(vapply(x, n_elements, integer(1)))
  } else
  {
    as.integer(prod(DIM(x)))
  }  
}

has_elements <- function(x, n, .xname = get_name_in_parent(x))
{
  n <- assertive.base::use_first(n)
  check_n(n)
  n_elements_x <- n_elements(x)
  if(n_elements_x != n)
  {
    return(
      false(
        ngettext(
          n_elements_x, 
          "%s has %d element, not %d.", 
          "%s has %d elements, not %d."
        ),
        .xname, 
        n_elements_x,
        n
      )
    )
  }
  TRUE
}

check_n <- function(n)
{
  if(any(n < 0 | n != round(n)))
  {
    stop("n should be a non-negative integer vector.")
  }
}

get_metric <- function(metric)
{
  switch(
    metric,
    length   = is_of_length,
    elements = has_elements,
    stop("Bug in assertive; the metric", metric, "is not valid.", domain = NA)
  )
}

#' Is the input empty/scalar?
#'
#' Checks to see if the input has length zero/one.
#'
#' @param x Input to check.
#' @param n Non-negative integer(s) of the expected length/number of elements/
#' lengths of dimensions.  See note.
#' @param metric A string. Should be length or the number of elements be used to
#' determine if the object is empty/non-empty/scalar?
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_empty} returns \code{TRUE} if the input has length 
#' zero.  \code{is_scalar} returns \code{TRUE} if the input has length 
#' one.  The \code{assert_*} functions return nothing but throw an
#' error if the corresponding \code{is_*} function returns \code{FALSE}.
#' @note For \code{is_empty}, \code{is_non_empty} and \code{is_scalar}, \code{n}
#' should be an single integer representing either the expected length or the
#' expected number of elements in \code{x}.  For \code{is_of_dimension} \code{n}
#' should be a vector of integers representing the expected lengths of 
#' dimensions.
#' @seealso \code{\link{length}}.
#' @examples
#' # is_of_length returns TRUE if the length of an object
#' # matches a specified value.
#' is_of_length(1:5, 5)
#' assert_is_of_length(1:5, 5)
#' 
#' # has_elements returns TRUE if an object has a specified
#' # number of elements.  This is usually the same thing.
#' has_elements(1:5, 5)
#' assert_has_elements(1:5, 5)
#' 
#' # Data frames and lists behave differently for length
#' # and number of elements.
#' d <- data.frame(x = 1:5, y = letters[1:5])
#' assert_is_of_length(d, 2)
#' assert_has_elements(d, 10)
#' 
#' l <- list(a = 1:5, b = list(b.a = 1:3, b.b = 1:7))
#' assert_is_of_length(l, 2)
#' assert_has_elements(l, 15)
#' 
#' # Functions always have length one, but may have lots of 
#' # elements.
#' assert_is_of_length(var, 1)
#' assert_has_elements(var, 54)
#' 
#' # is_scalar is a shortcut for length one, or one elements.
#' assert_is_scalar(99)
#' assert_is_scalar("Multiple words in a single string are scalar.")
#' assert_is_scalar(NA)
#' 
#' # The two metrics can yield different results!
#' is_scalar(list(1:5))
#' is_scalar(list(1:5), "elements")
#' is_scalar(var)
#' is_scalar(var, "elements")
#' 
#' # Similarly, is_empty is a shortcut for length zero/zero elements.
#' assert_is_empty(NULL)
#' assert_is_empty(numeric())
#' assert_is_non_empty(1:10)
#' assert_is_non_empty(NA)
#' 
#' # is_of_dimension tests the lengths of all dimensions.
#' assert_is_of_dimension(d, c(5, 2))
#' assert_is_of_dimension(l, NULL)
#' @author Richard Cotton
#' @noRd
is_empty <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x))
{  
  metric <- match.arg(metric)
  metric_fn <- get_metric(metric)
  metric_fn(x, 0L, .xname)
}

is_scalar <- function(x, metric = c("length", "elements"), 
                      .xname = get_name_in_parent(x))
{
  metric <- match.arg(metric)
  metric_fn <- get_metric(metric)
  metric_fn(x, 1L, .xname)
}


assert_is_empty <- function(x, metric = c("length", "elements"), 
  severity = getOption("assertive.severity", "stop"))
{                             
  metric <- match.arg(metric)                             
  assert_engine(
    is_empty, 
    x, 
    metric = metric, 
    .xname = get_name_in_parent(x),
    severity = severity
  ) 
}


assert_is_scalar <- function(x, metric = c("length", "elements"), 
  severity = getOption("assertive.severity", "stop"))
{                                        
  metric <- match.arg(metric)
  assert_engine(
    is_scalar, 
    x, 
    metric = metric, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#=======================================================================
#
#                      assertive.types
#
#=======================================================================


#--------------------------
#
#     assert_is_all_of
#     assert_is_any_of
#
#--------------------------

#' Does x belong to these classes?
#' 
#' Checks to see if x belongs to any of the classes in classes.
#' @param x Input to check.
#' @param classes As for \code{class}. 
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return The functions return nothing but throw an error if
#' \code{x} does not have any/all of the class \code{classes}.
#' @seealso \code{\link[assertive.base]{is2}}
#' @examples 
#' assert_is_all_of(1:10, c("integer", "numeric"))
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_any_of(1:10, c("list", "data.frame")))
#' @author Richard Cotton
#' @noRd
assert_is_all_of <- function(x, classes, 
  severity = getOption("assertive.severity", "stop"))
{  
  msg <- gettextf(
    "%s is not in all of the classes %s.", 
    get_name_in_parent(x), 
    toString(sQuote(classes))
  )
  assert_engine(
    assertive.base::is2, 
    x, 
    class = classes, 
    msg = msg, 
    severity = severity
  )
}


assert_is_any_of <- function(x, classes, 
  severity = getOption("assertive.severity", "stop"))
{  
  msg <- gettextf(
    "%s is not in any of the classes %s.", 
    get_name_in_parent(x), 
    toString(sQuote(classes))
  )
  assert_engine(
    assertive.base::is2, 
    x, 
    class = classes, 
    msg = msg, 
    what = "any",
    severity = severity
  )
}


#------------------------
#
#    assert_is_a_bool
#
#------------------------

#' Is the input logical?
#'
#' Checks to see if the input is logical.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_logical} wraps \code{is.logical}, providing more 
#' information on failure. \code{is_a_bool} returns \code{TRUE} if the 
#' input is logical and scalar.  The \code{assert_*} functions return
#' nothing but throw an error if the corresponding \code{is_*} function
#' returns \code{FALSE}.
#' @seealso \code{\link[base]{is.logical}} and \code{\link{is_scalar}}.
#' @examples
#' assert_is_logical(runif(10) > 0.5)
#' assert_is_a_bool(TRUE)
#' assert_is_a_bool(NA)
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_logical(1))
#' assertive.base::dont_stop(assert_is_a_bool(c(TRUE, FALSE)))
#' assertive.base::dont_stop(assert_is_a_bool(logical()))
#' @author Richard Cotton
#' @noRd
is_logical <- function(x, .xname = get_name_in_parent(x))
{
  assertive.base::is2(x, "logical", .xname)
}       

# assertive.properties (!)
is_scalar <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x))
{
  metric <- match.arg(metric)
  metric_fn <- get_metric(metric)
  metric_fn(x, 1L, .xname)
}     

is_a_bool <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_logical(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
}

assert_is_a_bool <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{      
  assert_engine(
    is_a_bool, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  )
}

assert_is_logical <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_logical, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


#-----------------------------
#
#      assert_is_a_number
#
#------------------------------

#' Is the input numeric?
#'
#' Checks to see if the input is numeric.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_numeric} wraps \code{is.numeric}, providing more 
#' information on failure. \code{is_a_number} returns \code{TRUE} if the 
#' input is numeric and scalar.  The \code{assert_*} functions return nothing
#' but throw an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @note \code{numeric} means either double or integer, inc this case.
#' @seealso \code{\link{is_integer}}, \code{\link[base]{is.numeric}} and 
#' \code{\link{is_scalar}}.
#' @examples
#' # "numeric" fns work on double or integers; 
#' assert_is_numeric(1:10)
#' 
#' # Here we check for length 1 as well as type
#' assert_is_a_number(pi)
#' assert_is_a_number(1L)
#' assert_is_a_number(NA_real_)
#' 
#' # "double" fns fail for integers.
#' assert_is_a_double(pi)
#' 
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_numeric(c(TRUE, FALSE)))
#' assertive.base::dont_stop(assert_is_a_number(1:10))
#' assertive.base::dont_stop(assert_is_a_number(numeric()))
#' assertive.base::dont_stop(assert_is_double(1:10))
#' @author Richard Cotton
#' @noRd
is_numeric <- function(x, .xname = get_name_in_parent(x))
{
  assertive.base::is2(x, "numeric", .xname)
}


is_a_number <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_numeric(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
} 


assert_is_a_number <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                          
  assert_engine(
    is_a_number, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  )
}


assert_is_numeric <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_numeric, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


#-----------------------------
#
#      assert_is_a_string
#      assert_is_character
#
#------------------------------

#' Is the input of type character?
#'
#' Checks to see if the input is of type character.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_character} wraps \code{is.character}, providing more 
#' information on failure. \code{is_a_string} returns \code{TRUE} if the 
#' input is character and scalar.
#' The \code{assert_*} functions return nothing but throw an error if the
#' corresponding \code{is_*} function returns \code{FALSE}.
#' @seealso \code{\link[base]{is.character}} and \code{\link{is_scalar}}.
#' @examples
#' assert_is_character(letters)
#' assertive.base::dont_stop(assert_is_character(factor(letters)))
#' @author Richard Cotton
#' @noRd
is_character <- function(x, .xname = get_name_in_parent(x))
{
  assertive.base::is2(x, "character", .xname)
}


is_a_string <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_character(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
}


assert_is_a_string <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_a_string, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  )
}

assert_is_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_character, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


#-----------------------------
#
#      assert_is_data_frame
#
#------------------------------

#' Is the input is a data.frame?
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_data.frame} wraps \code{is.data.frame}, 
#' providing more information on failure.  \code{assert_is_data.frame} 
#' returns nothing but throws an error if \code{is_data.frame} 
#' returns \code{FALSE}.
#' @seealso \code{\link[base]{is.data.frame}}.
#' @examples
#' assert_is_data.frame(data.frame())
#' assert_is_data.frame(datasets::CO2)
#' @author Richard Cotton
#' @noRd
is_data.frame <- function(x, .xname = get_name_in_parent(x))
{
  assertive.base::is2(x, "data.frame", .xname)
}


assert_is_data.frame <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_data.frame, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


#-----------------------------
#
#      assert_is_function
#
#------------------------------

#' Is the input a function?
#'
#' Checks to see if the input is a function.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_function}, \code{is_primitive} and \code{is_stepfun}
#' wrap \code{is.function}, \code{is.primitive} and \code{is.stepfun} 
#' repsectively, providing more information on failure.  The 
#' \code{assert_*} functions return nothing but throw an error if the
#' corresponding \code{is_*} function returns \code{FALSE}.
#' @seealso \code{\link[base]{is.function}}.
#' @examples
#' assert_is_function(sqrt)
#' assert_is_function(function(){})
#' @author Richard Cotton
#' @noRd
is_function <- function(x, .xname = get_name_in_parent(x))
{
  assertive.base::is2(x, "function", .xname)
}

assert_is_function <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_function, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


#-----------------------------
#
#      assert_is_list
#
#------------------------------

#' Is the input a list?
#'
#' Checks to see if the input is a list.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_list} wraps \code{is.list}, providing more 
#' information on failure.
#' @seealso \code{\link[base]{is.list}}.
#' @examples
#' assert_is_list(list(1,2,3))
#' assert_is_pairlist(.Options)
#' #These examples should fail.
#' assertive.base::dont_stop({
#'   assert_is_list(1:10)
#'   assert_is_pairlist(options())
#' })
#' @author Richard Cotton
#' @noRd
is_list <- function(x, .xname = get_name_in_parent(x))
{
  assertive.base::is2(x, "list", .xname)
}

assert_is_list <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_list, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


#-----------------------------
#
#      assert_is_matrix
#
#------------------------------

#' Is the input an array or matrix?
#'
#' Checks to see if the input is an array or matrix.
#'
#' @param x Input to check.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return \code{is_array} and \code{is_matrix} wrap \code{is.array}, 
#' and \code{is.matrix} respectively, providing more information on
#' failure.  The \code{assert_*} functions return nothing but throw
#' an error if the corresponding \code{is_*} function returns
#' \code{FALSE}.
#' @examples
#' assert_is_array(array())
#' assert_is_array(matrix())
#' assert_is_matrix(matrix())
#' #These examples should fail.
#' assertive.base::dont_stop(assert_is_matrix(array()))
#' @author Richard Cotton
#' @noRd
is_array <- function(x, .xname = get_name_in_parent(x))
{
  assertive.base::is2(x, "array", .xname)
}

is_matrix <- function(x, .xname = get_name_in_parent(x))
{
  assertive.base::is2(x, "matrix", .xname)
}


assert_is_matrix <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_matrix, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


#-----------------------------
#
#      is_formula
#
#------------------------------

#' Is the input a formula?
#'
#' Checks to see if the input is a formula.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return The \code{is_*} functions return \code{TRUE} when the input is a 
#' formula.  The \code{assert_*} functions return nothing but throw an error 
#' if the corresponding \code{is_*} function returns \code{FALSE}.
#' @seealso \code{\link{is_environment}} and  \code{\link{is_language}}
#' @examples
#' is_one_sided_formula(~ x)
#' is_two_sided_formula(y ~ x)
#' @author Richard Cotton
#' @noRd
is_formula <- function(x, .xname = get_name_in_parent(x))
{
  assertive.base::is2(x, "formula", .xname)
}


#=======================================================================
#
#                      assertive.strings
#
#=======================================================================


#' Does the input contain empty or missing strings?
#' 
#' Checks for empty or missing strings.
#' @param x A character vector.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return The \code{is_*} functions return logical vectors for strings which
#' are (non) empty or missing, and the \code{assert_*} functions throw errors
#' on failure.
#' @note In R, \code{NA_character_} is considered to be a non-empty string
#' (at least by \code{\link[base]{nzchar}}), which is why many functions are
#' needed to to clarify the situation.
#' @seealso \code{is_character}, 
#' \code{\link[base]{nzchar}}
#' @examples 
#' # These functions return a vector:
#' x <- c("", "a", NA)
#' is_empty_character(x)
#' is_non_empty_character(x)
#' is_missing_or_empty_character(x)
#' is_non_missing_nor_empty_character(x)
#' 
#' # These functions return a single value:
#' is_an_empty_string("")
#' is_an_empty_string("a")
#' is_an_empty_string(NA_character_)
#' 
#' is_a_non_empty_string("")
#' is_a_non_empty_string("a")
#' is_a_non_empty_string(NA_character_)
#' 
#' is_a_missing_or_empty_string("")
#' is_a_missing_or_empty_string("a")
#' is_a_missing_or_empty_string(NA_character_)
#' 
#' is_a_non_missing_nor_empty_string("")
#' is_a_non_missing_nor_empty_string("a")
#' is_a_non_missing_nor_empty_string(NA_character_)
#' @author Richard Cotton
#' @noRd
is_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- assertive.base::coerce_to(x, "character", .xname)
  assertive.base::call_and_name(
    function(x)
    {
      ok <- !nzchar(x)
      assertive.base::set_cause(ok, ifelse(is.na(x), "missing", "nonempty"))
    },
    x
  )
}


is_non_missing_nor_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- assertive.base::coerce_to(x, "character", .xname)
  ok <- nzchar(x) & !assertive.base::is_na(x)
  assertive.base::set_cause(ok, ifelse(is.na(x), "missing", "empty"))
}


assert_all_are_non_missing_nor_empty_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are not all non-missing nor non-empty strings.", 
    .xname
  )
  assert_engine(
    is_non_missing_nor_empty_character, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}



#=======================================================================
#
#                      autonomics
#
#=======================================================================


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
#' file <- download_data('billing19.proteingroups.txt')
#' object <- read_proteingroups(file, plot=FALSE)
#' contains_ratios(object)
#'
#' file <- download_data('halama18.metabolon.xlsx')
#' object <- read_metabolon(file, plot=FALSE)
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
    if (!(ok <- assertive.base::is2(x, "SummarizedExperiment")))  return(ok)
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







