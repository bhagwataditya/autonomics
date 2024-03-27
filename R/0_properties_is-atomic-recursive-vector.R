is_atomic <- function(x, .xname = get_name_in_parent(x))
{
  if(!is.atomic(x))
  {
    return(false(gettext("%s is not atomic."), .xname))
  }
  TRUE
}

is_nested <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_recursive(x, .xname)))
  {
    return(ok)
  }
  are_rec <- bapply(x, is.recursive)
  if(!any(are_rec))
  {
    return(false(gettext("%s has no recursive elements."), .xname))
  }
  TRUE
}

is_non_nested <- function(x, .xname = get_name_in_parent(x))
{
  are_rec <- bapply(x, is.recursive)
  if(any(are_rec))
  {
    msg <- ngettext(
      sum(are_rec),
      "Element %s of %s is recursive.",
      "Elements %s of %s are recursive."
    )
    return(false(msg, toString(which(are_rec), width = 50), .xname))
  }
  TRUE
}

is_recursive <- function(x, .xname = get_name_in_parent(x))
{
  if(!is.recursive(x))
  {
    return(false("%s is not recursive.", .xname))
  }
  TRUE
}

is_vector <- function(x, .xname = get_name_in_parent(x))
{
  if(!is.vector(x)) 
  {
    return(false("%s is not a vector.", .xname))
  }
  TRUE
}                
