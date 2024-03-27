has_duplicates <- function(x, .xname = get_name_in_parent(x))
{
  if(!anyDuplicated(x)) 
  {
    return(false(gettext("%s has no duplicates."), .xname))
  }
  TRUE
}

has_no_duplicates <- function(x, .xname = get_name_in_parent(x))
{
  if(anyDuplicated(x)) 
  {
    dupe_indicies <- which(duplicated(x))
    return(
      false(
        ngettext(
          length(dupe_indicies),
          "%s has a duplicate at position %s.",
          "%s has duplicates at positions %s."
        ), 
        .xname, 
        toString(dupe_indicies, width = 100)
      )
    )
  }
  TRUE
}

