assert_has_cols <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                       
  assert_engine(
    has_cols, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

assert_has_dims <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                                
  assert_engine(
    has_dims, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

assert_has_rows <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                           
  assert_engine(
    has_rows, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

