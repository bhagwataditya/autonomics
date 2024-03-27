is_a_bool <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_logical(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
}

is_a_complex <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_complex(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
}

is_a_double <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_double(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
} 

is_a_number <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_numeric(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
} 

is_a_raw <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_raw(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
} 

is_a_string <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_character(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
}

is_an_integer <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_integer(x, .xname))) return(ok)
  if(!(ok <- is_scalar(x, .xname = .xname))) return(ok)
  TRUE
} 

is_inherited_from <- function(x, classes, .xname = get_name_in_parent(x))
{
  ok <- bapply(classes, function(class) inherits(x, class))
  if(!any(ok)) 
  {
    msg <- ngettext(
      length(classes),
      "%s does not inherit from the class %s. It has class %s.",
      "%s does not inherit from any of the classes %s. It has class %s."
    )
    return(
      false(msg, .xname, toString(classes), toString(class(x)))
    )
  }
  TRUE
}
