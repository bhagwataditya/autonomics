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

is_logical_string <- function(x, .xname)
{
  x <- coerce_to(x, "character", .xname)
  ok <- call_and_name(
    function(x)
    {
      suppressWarnings(
        {
          logx <- as.logical(x)
          is_not_na(logx)
        }
      )
    },
    x
  )
  set_cause(ok, ifelse(is.na(x), "missing", "bad format"))
}

is_single_character <- function(x, .xname)
{
  x <- coerce_to(x, "character", .xname)
  ok <- call_and_name(
    function(x)
    {
      nch_is_1 <- nchar(x) == 1
      is_na_x <- is.na(x)
      if(any(is_na_x))
      {
        message("New behaviour: NA inputs now return NA.")
        nch_is_1[is_na_x] <- NA
      }
      nch_is_1
    },
    x
  )
  set_cause(ok, ifelse(is.na(x), "missing", "bad format"))  
}
