
is_vector <- function(x, .xname = get_name_in_parent(x))
{
  if(!is.vector(x)) 
  {
    return(false("%s is not a vector.", .xname))
  }
  TRUE
}                
