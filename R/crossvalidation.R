

affect_to_rhg <- function(p, phatmax, phat){

  temp_phatmax <- phatmax[order(phatmax)]
  temp_phat    <- phat[order(phatmax)]

  indice <- 1:length(temp_phatmax)

  temp_phatmin <- c(-0.1,temp_phatmax[1:(length(temp_phatmax)-1)])

  temp_phatmax[indice == length(temp_phatmax)] <- 1.1

  pouet <- mapply(min=temp_phatmin,max=temp_phatmax,FUN=compare,MoreArgs=list(x=p))

  return(temp_phat[unlist(apply(pouet,1,FUN=which))])

}


compute_error_rhg <- function(ly, lphat, lrhg, lweight, ty, tphat, method = 'mse'){

  temp <- data.table(phat = lphat,
                     weight = lweight,
                     rhg    = lrhg)

  temp_rhg <-
    temp[,list(phatmax_rhg = max(phat),
               phat_rhg = stats::weighted.mean(phat, weights = weight)),
         by = rhg]

  phathat <- affect_to_rhg(tphat, phatmax = temp_rhg$phatmax_rhg, phat = temp_rhg$phat_rhg)

  if (method == 'mse') {
    error <- sqrt(sum((ty - phathat)**2) / length(tphat))
  } else if (method == 'mad') {
    error <- sum(abs(ty - phathat)) / length(tphat)
  } else if (method == 'logloss') {
    error <- - sum(ty*log(phathat) + (1-ty)*log(1-phathat)) / length(tphat)
  }
  error
}

# compute_error_hb #

compute_error_hb <- function(y, phat, weight, learning_identifier, parameters,
                             method = 'mse',
                             output = 'matrix',
                             print_rhg = T){

  rhg_hb <- hb(phat[learning_identifier], parameters, print_rhg = F)

  if (output == 'matrix'){

    results <-
      matrix(c(sapply(rhg_hb$q$rhg,
                      FUN = function(r) compute_error_rhg(r,
                                                          ly      = y[learning_identifier],
                                                          lphat   = phat[learning_identifier],
                                                          lweight = weight[learning_identifier],
                                                          ty      = y[!learning_identifier],
                                                          tphat   = phat[!learning_identifier],
                                                          method  = method)),
               sapply(rhg_hb$km$rhg,
                      FUN     = function(r) compute_error_rhg(r,
                                                              ly      = y[learning_identifier],
                                                              lphat   = phat[learning_identifier],
                                                              lweight = weight[learning_identifier],
                                                              ty      = y[!learning_identifier],
                                                              tphat   = phat[!learning_identifier],
                                                              method  = method))),
             ncol = 2)

    rownames(results) <- rhg_hb$q$parameters
    colnames(results) <- c('q', 'km')

  } else {

    results <- c(sapply(rhg_hb$q$rhg,
                        FUN = function(r) compute_error_rhg(r,
                                                            ly      = y[learning_identifier],
                                                            lphat   = phat[learning_identifier],
                                                            lweight = weight[learning_identifier],
                                                            ty      = y[!learning_identifier],
                                                            tphat   = phat[!learning_identifier],
                                                            method  = method)),
                 sapply(rhg_hb$km$rhg,
                        FUN     = function(r) compute_error_rhg(r,
                                                                ly      = y[learning_identifier],
                                                                lphat   = phat[learning_identifier],
                                                                lweight = weight[learning_identifier],
                                                                ty      = y[!learning_identifier],
                                                                tphat   = phat[!learning_identifier],
                                                                method  = method)))

    names(results) <- c(paste0('q_',as.character(rhg_hb$q$parameters)),
                        paste0('km_',as.character(rhg_hb$km$parameters)))

  }

  results

}


select_hb <- function(y, phat, weight = rep(1, length(y)), K = 10, parameters, method = 'mse',
                      print_rhg = F){

  # Génération des échantillons pour la validation croisée

  which <- rep(seq_len(K), length.out = length(y))
  samples_cv <- lapply(1:K, FUN = function(x){ which != x})

  # Application de la validation croisée

  cv <- apply(
    sapply(samples_cv,
           FUN = function(l) compute_error_hb(learning_identifier = l,
                                              y = y,
                                              phat = phat,
                                              weight = weight,
                                              parameters = parameters,
                                              method = method,
                                              output = 'vector',
                                              print_rhg = print_rhg)),
    MARGIN = 1,
    FUN = mean)

  type <- rep(c('q','km'), each = length(parameters))

  sorted_parameters <- parameters[order(parameters)]

  c(q = sorted_parameters[which.min(cv[type == 'q'])],
    km = sorted_parameters[which.min(cv[type == 'km'])])

}

################################# Fonction choose_best_hb #################################


choose_best_hb <- function(y, phat, weight = rep(1, length(y)), u, parameters,
                           method = 'mse', K = 10, print_rhg = F){

  optimal_parameters <- select_hb(y = y[u],
                                  phat = phat[u],
                                  weight = weight[u],
                                  K = K,
                                  parameters = parameters,
                                  method = method,
                                  print_rhg = print_rhg)

  temp <-
    data.table(rhg_q = hb(phat[u], parameters = optimal_parameters['q'], print_rhg = F)$q$rhg[[1]],
               rhg_km = hb(phat[u], parameters = optimal_parameters['km'], print_rhg = F)$km$rhg[[1]],
               phat = phat[u],
               weight = weight[u])

  temp_rhg_q <-
    temp[,list(phatmax_rhg = max(phat),phat_rhg = stats::weighted.mean(phat, weights = weight)),
         by = rhg_q]

  temp_rhg_km <-
    temp[,list(phatmax_rhg = max(phat),phat_rhg = stats::weighted.mean(phat, weights = weight)),
         by = rhg_km]

  list(
    q = list(param_opt = optimal_parameters['q'],
             phat = affect_to_rhg(phat[!u],
                                  phatmax = temp_rhg_q$phatmax_rhg,
                                  phat = temp_rhg_q$phat_rhg)),
    km = list(param_opt = optimal_parameters['km'],
              phat = affect_to_rhg(phat[!u],
                                   phatmax = temp_rhg_q$phatmax_rhg,
                                   phat = temp_rhg_q$phat_rhg))
  )

}

