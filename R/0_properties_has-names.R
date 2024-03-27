has_colnames <- function(x, .xname = get_name_in_parent(x))
{
  colnamesx <- colnames(x)
  if(is.null(colnamesx)) 
  {
    return(false("The column names of %s are NULL.", .xname))
  }
  if(!any(nzchar(colnamesx))) 
  {
    return(false("The column names of %s are all empty.", .xname))
  }
  TRUE
} 

has_dimnames <- function(x, .xname = get_name_in_parent(x))
{
  dimnamesx <- dimnames(x)
  if(is.null(dimnamesx)) 
  {
    return(false("The dimension names of %s are NULL.", .xname))
  }
  if(!any(nzchar(unlist(dimnamesx, use.names = FALSE)))) 
  {
    return(false("The dimension names of %s are all empty.", .xname))
  }
  TRUE
} 

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

has_rownames <- function(x, .xname = get_name_in_parent(x))
{
  rownamesx <- rownames(x)
  if(is.null(rownamesx)) 
  {
    return(false("The row names of %s are NULL.", .xname))
  }
  if(!any(nzchar(rownamesx))) 
  {
    return(false("The row names of %s are all empty.", .xname))
  }
  TRUE
} 

