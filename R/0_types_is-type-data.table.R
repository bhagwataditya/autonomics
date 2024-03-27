is_data.table <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_data.frame(x, .xname)))
  {
    return(ok)
  }
  is2(x, "data.table", .xname)
}
