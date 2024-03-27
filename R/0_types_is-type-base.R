is_array <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "array", .xname)
}


is_call <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "call", .xname)
}


is_character <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "character", .xname)
}


is_complex <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "complex", .xname)
}       


is_data.frame <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "data.frame", .xname)
}


is_double <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "double", .xname)
}


is_environment <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "environment", .xname)
}

is_expression <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "expression", .xname)
}


is_externalptr <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "externalptr", .xname)
}


is_factor <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "factor", .xname)
}


is_function <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "function", .xname)
}


is_integer <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "integer", .xname)
}

is_language <- function(x, .xname = get_name_in_parent(x))
{
  if(!is.language(x)) 
  {
    return(
      false(
        gettext("%s is not a language object (name, call or expression)."), 
        .xname
      )
    )
  }
  TRUE
}



is_list <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "list", .xname)
}



is_logical <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "logical", .xname)
}       



is_matrix <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "matrix", .xname)
}



is_name <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "name", .xname)
}



is_numeric <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "numeric", .xname)
}


is_ordered <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_factor(x))) 
  {
    return(ok)
  }
  if(!is.ordered(x))
  {
    return(false(gettext("%s is not an ordered factor."), .xname))
  }
  TRUE
}


is_pairlist <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "pairlist", .xname)
}


is_primitive <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_function(x))) return(ok)
  if(!is.primitive(x))
  {
    return(false(gettext("%s is not a primitive function."), .xname))
  }
  TRUE
} 


is_qr <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "qr", .xname)
}


is_raw <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "raw", .xname)
}


is_s4 <- function(x, .xname = get_name_in_parent(x))
{
  if(!isS4(x))
  {
    return(false(gettext("%s is not an S4 object."), .xname))
  }
  TRUE
} 


is_S4 <- function(x, .xname = get_name_in_parent(x))
{
  .Deprecated("is_s4")
  is_s4(x, .xname)
}


is_symbol <- is_name


is_table <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "table", .xname)
}
