

assert_is_data.frame <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_data.frame, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


assert_is_factor <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_factor, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


assert_is_function <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_function, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


assert_is_integer <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_integer, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


assert_is_list <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_list, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


assert_is_logical <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_logical, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


assert_is_matrix <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_matrix, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}


