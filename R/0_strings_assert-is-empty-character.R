

assert_all_are_non_missing_nor_empty_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are not all non-missing nor non-empty strings.", 
    .xname
  )
  assert_engine(
    is_non_missing_nor_empty_character, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}

