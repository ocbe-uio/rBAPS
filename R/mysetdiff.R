mysetdiff <- function(A, B) {
 # MYSETDIFF Set difference of two sets of positive integers (much faster than
 # built - in setdiff)
 # C <- mysetdiff(A, B)
 # C <- A \ B = { things in A that are not in B }
 # MATLAB Original by Kevin Murphy, modified by Leon Peshkin
  if (is.null(A)) {
    return(vector())
  } else if (is.null(B)) {
    return(A)
  } else { # both non-empty
    bits <- zeros(1, base::max(base::max(A), base::max(B)))
    bits[A] <- 1
    bits[B] <- 0
    C <- A[as.logical(bits[A])]
  }
    return(C)
}
