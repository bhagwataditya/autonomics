is_formula <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "formula", .xname)
}

is_one_sided_formula <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_formula(x, .xname)))
  {
    return(ok)
  }
  if(!(ok <- is_of_length(x, 2L, .xname)))
  {
    return(ok)
  } 
  TRUE
}

is_two_sided_formula <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_formula(x, .xname)))
  {
    return(ok)
  }
  if(!(ok <- is_of_length(x, 3L, .xname)))
  {
    return(ok)
  } 
  TRUE
}
