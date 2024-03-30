assert_is_a_bool <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{      
  assert_engine(
    is_a_bool, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  )
}

assert_is_a_number <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                          
  assert_engine(
    is_a_number, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  )
}

assert_is_a_string <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_a_string, 
    x, 
    .xname = get_name_in_parent(x), 
    severity = severity
  )
}

