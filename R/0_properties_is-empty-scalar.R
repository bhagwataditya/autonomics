is_empty <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x))
{  
  metric <- match.arg(metric)
  metric_fn <- get_metric(metric)
  metric_fn(x, 0L, .xname)
}


is_non_empty <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x))
{
  metric <- match.arg(metric)
  metric_fn <- get_metric(metric)
  if(metric_fn(x, 0)) 
  {
    msg <- switch(
      metric,
      length = gettext("%s has length 0."),
      elements = gettext("%s has 0 elements.")
    )
    return(false(msg, .xname))
  }
  TRUE
}

is_non_scalar <- function(x, metric = c("length", "elements"), .xname = get_name_in_parent(x))
{
  metric <- match.arg(metric)
  metric_fn <- get_metric(metric)
  if(metric_fn(x, 1)) 
  {
    msg <- switch(
      metric,
      length = gettext("%s has length 1."),
      elements = gettext("%s has 1 element.")
    )
    return(false(msg, .xname))
  }
  TRUE
}

is_scalar <- function(x, metric = c("length", "elements"), 
                      .xname = get_name_in_parent(x))
{
  metric <- match.arg(metric)
  metric_fn <- get_metric(metric)
  metric_fn(x, 1L, .xname)
}     


has_elements <- function(x, n, .xname = get_name_in_parent(x))
{
  n <- use_first(n)
  check_n(n)
  n_elements_x <- n_elements(x)
  if(n_elements_x != n)
  {
    return(
      false(
        ngettext(
          n_elements_x, 
          "%s has %d element, not %d.", 
          "%s has %d elements, not %d."
        ),
        .xname, 
        n_elements_x,
        n
      )
    )
  }
  TRUE
}


is_of_length <- function(x, n, .xname = get_name_in_parent(x))
{
  n <- use_first(n)
  check_n(n)
  length_x <- length(x)
  if(length_x != n)
  {
    return(false("%s has length %d, not %d.", .xname, length_x, n))
  }
  TRUE
}

check_n <- function(n)
{
  if(any(n < 0 | n != round(n)))
  {
    stop("n should be a non-negative integer vector.")
  }
}

get_metric <- function(metric)
{
  switch(
    metric,
    length   = is_of_length,
    elements = has_elements,
    stop("Bug in assertive; the metric", metric, "is not valid.", domain = NA)
  )
}
