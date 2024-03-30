
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

