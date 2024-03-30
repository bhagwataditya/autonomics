are_same_length <- function(x, y, .xname = get_name_in_parent(x),
  .yname = get_name_in_parent(y))
{
  len_x <- length(x)
  len_y <- length(y)
  if(len_x != len_y)
  {
    return(
      false(
        gettext("%s has length %d but %s has length %d."),
        .xname,
        len_x,
        .yname,
        len_y
      )
    )
  }
  TRUE
}

