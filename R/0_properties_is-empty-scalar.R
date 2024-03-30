
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
