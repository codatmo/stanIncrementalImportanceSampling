#' Calculate Importance Weights
#'
#' This function calculates importance weights for a previously fitted Stan
#' model given a new set of input data. Given a stanfit object with MCMC samples
#' and a new set of input data for the stan model, calculate importance sampling
#' weights. For example, a time series model predicting into the future can be
#' updated with new data post-hoc when data becomes available.
#'
#' @param originalFit stanfit containing MCMC samples
#' @param oldData list() containing data used to fit originalFit
#' @param newData list() containing the new input data
#' @param model_code string containing Stan model code for originalFit (optional)
#'
#' @return list containing weights, neff, lpOriginal and lpNew
#' @export
importance_sampling <- function(originalFit, oldData, newData, model_code = originalFit@stanmodel@model_code[1]){

  #Compile model with new data
  newFit <- rstan::stan(model_code = model_code, chains = 0, data = newData)

  #Have to recompile the model to use the log_prob function
  originalFitTmp <- rstan::stan(model_code = model_code, chains = 0, data = oldData)

  sim <- originalFit@sim
  originalLogPosterior <- unlist(lapply(seq(1,sim$chains), function(j){ sim$samples[[j]]$lp__ }))

  #Get the parameter names from the original model (recompiled)
  parnames <- originalFitTmp@model_pars
  pardims <- originalFitTmp@par_dims
  num_upars <- rstan::get_num_upars(originalFitTmp)

  #MCMC samples from the original fit
  samples <- rstan::extract(originalFit)
  ndraws <- nrow(samples[[1]])

  #For each iteration, get the log posterior of the MCMC samples
  #under the original fit and under the new data
  lps <- lapply(seq(1,ndraws), function(iter){
    pars <- lapply(parnames, function(par){

      #This check for the dimensions of the parameters could be neater
      if (is.na(ncol(samples[[par]]))){
        x <- samples[[par]][iter]
      } else {
        if (length(dim(samples[[par]])) > 2){
          x <-  samples[[par]][iter,,]
        } else {
          x <- samples[[par]][iter,]
        }
      }

      x
    })
    names(pars) <- parnames

    #Parameters in unconstrained space
    upars <- rstan::unconstrain_pars(newFit, pars)

    #Original logposterior
    lpOriginal <- rstan::log_prob(originalFitTmp, upars, adjust_transform = T, gradient = T)

    #Logposterior in new model
    lpNew <- rstan::log_prob(newFit, upars, adjust_transform = T, gradient = T)

    list(lpOriginal, lpNew)
  })

  lpOriginal <- unlist(lapply(lps, function(lp){ lp[[1]] }))
  lpNew <- unlist(lapply(lps, function(lp){ lp[[2]] }))

  #logRatio of posteriors
  lpRatio <- lpNew - lpOriginal
  lpRatio <- lpRatio - matrixStats::logSumExp(lpRatio)

  #Importance weights
  weights <- exp(lpRatio)

  #Effective sample size
  neff <- (sum(weights))^2/(sum(weights^2))

  list(weights = weights, neff = neff, lpOriginal = lpOriginal, lpNew = lpNew)
}
