

is_missing_or_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- coerce_to(x, "character", .xname)
  ok <- !nzchar(x) | is_na(x)
  set_cause(ok, "nonempty")
}

is_non_missing_nor_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- coerce_to(x, "character", .xname)
  ok <- nzchar(x) & !is_na(x)
  set_cause(ok, ifelse(is.na(x), "missing", "empty"))
}

