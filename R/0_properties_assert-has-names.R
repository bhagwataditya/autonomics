assert_has_names <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                            
  assert_engine(
    has_names, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
