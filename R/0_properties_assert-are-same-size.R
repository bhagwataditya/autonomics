assert_are_same_length <- function(x, y, 
  severity = getOption("assertive.severity", "stop"))
{
  assert_engine(
    are_same_length,
    x, 
    y = y,
    .xname = get_name_in_parent(x),
    .yname = get_name_in_parent(y),
    severity = severity
  )
}
