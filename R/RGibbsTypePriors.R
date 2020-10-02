#' Computes the prior number of clusters for a Normalised Generalised Gamma process
#'
#' @param	k			The value for k
#' @param	n			The value for n
#' @param	beta		The value for beta
#' @param	sigma		The value for sigma
#' @param	prec		The number of bits of precision used during internal computation
#' @return	The density in k for the prior distribution on the number of clusters for a Normalised Generalised Gamma process
#' @export
#' @useDynLib	RGibbsTypePriors
pkn_ngg <- function(k, n, beta, sigma, prec = 53) {
  res <- .C("pkn_ngg_api",
            k = as.integer(k),
            n = as.integer(n),
            beta = as.double(beta),
            sigma = as.double(sigma),
            prec = as.integer(prec),
            pkn = as.double(0),
            radius = as.double(0),
            accu = as.integer(0)
  )
  return(list(pkn = res$pkn, radius = res$radius, accuracy_bits = res$accu))
}

#' Computes the weights for a Normalised Generalised Gamma process
#'
#' @param	n			The value for n
#' @param	k			The value for k
#' @param	beta		The value for beta
#' @param	sigma		The value for sigma
#' @param	prec		The number of bits of precision used during internal computation
#' @return	The Vnk weights for a Normalised Generalised Gamma process
#' @export
#' @useDynLib	RGibbsTypePriors
vnk_ngg <- function(n, k, beta, sigma, prec = 53) {
  res <- .C("vnk_ngg_api",
            k = as.integer(k),
            n = as.integer(n),
            beta = as.double(beta),
            sigma = as.double(sigma),
            prec = as.integer(prec),
            vnk = as.double(0),
            radius = as.double(0),
            accu = as.integer(0)
  )
  return(list(vnk = res$vnk, radius = res$radius, accuracy_bits = res$accu))
}

#' Computes the generalized binomial coefficients
#'
#' @param	n			The value for n
#' @param	k			The value for k
#' @param	sigma		The value for sigma
#' @param	prec		The number of bits of precision used during internal computation
#' @return	The value of the generalized binomial coefficient
#' @export
#' @useDynLib	RGibbsTypePriors
cnk_ngg <- function(n, k, sigma, prec = 53) {
  res <- .C("cnk_ngg_api",
            k = as.integer(k),
            n = as.integer(n),
            sigma = as.double(sigma),
            prec = as.integer(prec),
            cnk = as.double(0),
            radius = as.double(0),
            accu = as.integer(0)
  )
  return(list(cnk = res$cnk, radius = res$radius, accuracy_bits = res$accu))
}
