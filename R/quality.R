
#' @export

calcule_erreurs <- function(r, phat){

  return(
    list(
      mse         = sqrt(sum((r-phat)**2) / length(r)),
      mad         = sum(abs(r-phat)) / length(r),
      vp          = mean(r*phat),
      fp          = mean((1-r)*phat),
      fn          = mean(r*(1-phat)),
      vn          = mean((1-r)*(1-phat)),
      sensibilite = mean(r*phat) / (mean(r*phat)+mean(r*(1-phat))),
      specificite = mean((1-r)*(1-phat)) / (mean((1-r)*(1-phat)) + mean((1-r)*phat)),
      vpp         = mean(r*phat) / (mean(r*phat) + mean((1-r)*phat)),
      vpn         = mean((1-r)*(1-phat)) / (mean((1-r)*(1-phat)) + mean(r*(1-phat)))
    )
  )

}

#' @export

calcule_variances <- function(phat, x, p){

  list(v1 = sum((1-phat)/phat**2*x**2/p),
       v2 = sum((1-p)/p**2*x**2/phat**2)
  )

}
