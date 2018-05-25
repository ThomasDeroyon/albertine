

compare <- function(x,min,max){

  if (min == max){
    return(x == min)
  } else {
    return((min <= x) & (x < max))
  }
}
