

is_data.frame <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "data.frame", .xname)
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


is_list <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "list", .xname)
}



is_matrix <- function(x, .xname = get_name_in_parent(x))
{
  is2(x, "matrix", .xname)
}


