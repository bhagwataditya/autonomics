
is_matching_fixed <- function(x, pattern, opts_fixed = NULL, 
                              .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- stringi::stri_detect_fixed(x, pattern, opts_fixed = opts_fixed)
      set_cause(ok, gettextf("does not match '%s'", pattern))
    },
    x
  )
}


is_not_matching_fixed <- function(x, pattern, opts_fixed = NULL, 
                              .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- !stringi::stri_detect_fixed(x, pattern, opts_fixed = opts_fixed)
      set_cause(ok, gettextf("matches '%s'", pattern))
    },
    x
  )
}



is_matching_regex <- function(x, pattern, opts_regex = NULL, 
                              .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- stringi::stri_detect_regex(x, pattern, opts_regex = opts_regex)
      set_cause(ok, gettextf("does not match '%s'", pattern))
    },
    x
  )
}


is_not_matching_regex <- function(x, pattern, opts_regex = NULL, 
                                  .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- !stringi::stri_detect_regex(x, pattern, opts_regex = opts_regex)
      set_cause(ok, gettextf("matches '%s'", pattern))
    },
    x
  )
}