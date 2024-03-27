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

get_dim_string <- function(x)
{
  if(is.null(x)) "NULL" else toString(x)
}

have_same_dims <- function(x, y, .xname = get_name_in_parent(x),
  .yname = get_name_in_parent(y))
{
  dim_x <- dim(x)
  dim_y <- dim(y)
  if(!identical(dim_x, dim_y))
  {
    return(
      false(
        gettext("%s has dim %s but %s has dim %s."),
        .xname,
        get_dim_string(dim_x),
        .yname,
        get_dim_string(dim_y)
      )
    )
  }
  TRUE
}

are_same_length_legacy <- function(..., l = list())
{
  .Deprecated("are_same_length")
  envir <- parent.frame()
  inputs <- as.list(match.call())[-1]
  inputs_in_list <- as.list(inputs$l)[-1]
  inputs <- c(inputs[names(inputs) != "l"], inputs_in_list)
  input_pairs <- expand.grid(expr1 = inputs, expr2 = inputs)
  equality <- apply(
    input_pairs, 
    1, 
    function(row)
    {       
      with(
        row,         
        length(eval(expr1, envir = envir)) == length(eval(expr2, envir = envir))
      )
    }
  )
  matrix(
    equality,
    nrow     = length(inputs),
    dimnames = list(inputs, inputs) 
  )
}
