

#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_any_are_matching_regex <- function(x, pattern, opts_regex = NULL, 
                                          na_ignore = FALSE,
                                          severity = getOption("assertive.severity", "stop")){
  .xname <- get_name_in_parent(x)                                             
  .fixedname <- get_name_in_parent(pattern)
  msg <- sprintf(
    "%s does not match %s", 
    .xname, .fixedname
  )
  assert_engine(
    is_matching_regex, 
    x,
    pattern,
    opts_regex = opts_regex,
    .xname     = .xname,
    msg        = msg, 
    what       = 'any',
    na_ignore  = na_ignore,
    severity   = severity
  )
}

