is_numeric_string <- function(x, .xname)
{
  x <- coerce_to(x, "character", .xname)
  ok <- call_and_name(
    function(x)
    {
      suppressWarnings(
        {
          numx <- as.numeric(x)
          is_not_na(numx)
        }
      )
    },
    x
  )
  set_cause(ok, ifelse(is.na(x), "missing", "bad format"))
}

