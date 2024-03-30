
assert_is_empty <- function(x, metric = c("length", "elements"), 
  severity = getOption("assertive.severity", "stop"))
{                             
  metric <- match.arg(metric)                             
  assert_engine(
    is_empty, 
    x, 
    metric = metric, 
    .xname = get_name_in_parent(x),
    severity = severity
  ) 
}

assert_is_non_empty <- function(x, metric = c("length", "elements"), 
  severity = getOption("assertive.severity", "stop"))
{                            
  metric <- match.arg(metric)                                 
  assert_engine(
    is_non_empty, 
    x, 
    metric = metric, 
    .xname = get_name_in_parent(x),
    severity = severity
  )  
}

