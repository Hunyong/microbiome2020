# time calculator
tt <- function(s){
  if (s==1) {time.tmp <<- Sys.time() # record time
  } else if (s==2) { # calculate time
    return(data.frame(begin = time.tmp, end = Sys.time(), elapsed = Sys.time() - time.tmp))
  }
}

"%btw%" <- function(x, range) {(x >= range[1])&(x <= range[2])}
if (FALSE) {#example
  "%btw%"(1:5, c(2.4,4.5))
  1:5 %btw% c(2.4,4.5)  
}

