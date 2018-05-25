


hb <- function(phat, parameters, leaf = 100, print_rhg = T){

  NbGRH <- 2

  sorted_parameters <- sort(parameters)

  param_max  <- max(parameters)
  minsizeQ   <- 1000
  minsizeKM   <- 1000
  r2Q        <- 0
  r2KM       <- 0
  v_r2Q      <- c()
  v_r2KM     <- c()
  list_grhQ  <- list()
  list_grhKM <- list()
  last_grhQ <- rep(1, length(phat))
  last_grhKM <- rep(1, length(phat))

  while(((r2Q < param_max) & (minsizeQ >= leaf)) | ((r2KM < param_max) & (minsizeKM >= leaf)) ){

    groupsQ <-
      cut(phat,unique(quantile(phat,seq(0,1,by=1/NbGRH))),labels=FALSE,include.lowest=TRUE,right=TRUE)

    t <- data.table(
      phat = phat,
      groupsQ = groupsQ
    )

    kernel     <- t[, list(p = mean(phat)), by = groupsQ]
    KM         <- suppressWarnings(stats::kmeans(phat, kernel$p, iter.max = 30))
    if(KM$ifault == 4)
      KM <- stats::kmeans(phat, kernel$p, algorithm = "MacQueen", iter.max = 30)
    t$groupsKM <- KM$cluster

    r2Q <-
      1 - sum(t[,list(sqr = (phat - mean(phat))**2),by = list(groupsQ)]$sqr) /
      sum((t$phat - mean(t$phat))**2)

    r2KM <-
      1 - sum(t[,list(sqr = (phat - mean(phat))**2),by = list(groupsKM)]$sqr) /
      sum((t$phat - mean(t$phat))**2)

    minsizeKM <- min(table(KM$cluster))
    minsizeQ <- min(table(groupsQ))

    if (minsizeQ >= leaf){

      v_r2Q[(NbGRH-1)]       <- r2Q
      list_grhQ[[(NbGRH-1)]] <- t$groupsQ
      last_grhQ <- t$groupsQ

    } else {

      v_r2Q[(NbGRH-1)]       <- param_max
      list_grhQ[[(NbGRH-1)]] <- last_grhQ

    }

    if (minsizeKM >= leaf){

      v_r2KM[(NbGRH-1)]       <- r2KM
      list_grhKM[[(NbGRH-1)]] <- t$groupsKM
      last_grhKM <- t$groupsKM

    } else {

      v_r2KM[(NbGRH-1)]       <- param_max
      list_grhKM[[(NbGRH-1)]] <- last_grhKM

    }

    if (print_rhg == T){
      message("Nombre de GRH  : ",NbGRH," R2 Quantiles : ", round(r2Q, 6),
              " Minsize Quantiles : ", minsizeQ,
              " R2 Centres Mobiles : ",round(r2KM, 6),
              " Minsize Centres Mobiles : ", minsizeKM
      )
    }
    NbGRH <- NbGRH + 1

  }

  if (length(parameters) > 1){

    indicesQ <- apply(sapply(v_r2Q,FUN=function(x){(x >= sorted_parameters)}),
                      1,
                      FUN=function(x){min(which(x))}
    )

    indicesKM <- apply(sapply(v_r2KM,FUN=function(x){(x >= sorted_parameters)}),
                       1,
                       FUN=function(x){min(which(x))}
    )
  } else {

    indicesQ <- min(which(sapply(v_r2Q,FUN=function(x){(x >= sorted_parameters)})))
    indicesKM <- min(which(sapply(v_r2KM,FUN=function(x){(x >= sorted_parameters)})))

  }

  results <-
    list(q=list(nb_rhg = sapply(indicesQ,FUN=function(x){x+1}),
                parameters = sorted_parameters,
                rhg = list_grhQ[indicesQ]),
         km=list(nb_rhg = sapply(indicesKM,FUN=function(x){x+1}),
                 parameters = sorted_parameters,
                 rhg = list_grhKM[indicesKM])
    )

  return(results)

}
