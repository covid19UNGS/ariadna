#
# Test S4 classes
# 


person <- setClass(
  # Set the name for the class
  "person",
  
  # Define the slots
  slots = list(
    age = "integer",
    actual_place   = "integer",
    state   = "integer",
    fb_region = 'integer'
  ),
  
  # Set the default values for the slots. (optional)
  prototype=list(
    age = 0,
    actual_place = 0,
    state   = 1,           # suceptible
    fb_region = c(0,0)
  ),
  
  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  validity=function(object)
  {
    if(object@age > 100)  {
      return("The age should be less than 100.")
    }
    return(TRUE)
  }
)

    
a <- person(10)
  a    
