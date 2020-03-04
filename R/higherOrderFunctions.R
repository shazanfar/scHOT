
##############################################

#' the weightedZIKendall function calculates weighted Tau*,
#' where Tau* is described in Pimentel et al (2015)
#' doi:10.1016/j.spl.2014.09.002. This association measure
#' is defined for zero-inflated, non-negative random variables.
#'
#' @title weightedZIKendall
#' @param x x and y are non-negative data vectors
#' @param y x and y are non-negative data vectors
#' @param w weight vector, values should be between 0 and 1
#' @return \code{numeric} weighted Tau* association value between x and y

#' @examples
#'
#' x = pmax(0,rnorm(100))
#' y = pmax(0,rnorm(100))
#' w = runif(100)
#' weightedZIKendall(x,y,w)
#'
#' @export

weightedZIKendall <- function(x, y, w = 1) {

  if (any(x < 0 | y < 0)) stop("x and/or y values have negative values")

  if (length(x) != length(y)) stop("x and y should have the same length")

  if (length(w) == 1) w <- rep(w, length(x))

  posx = x > 0
  posy = y > 0

  pospos = posx & posy

  x_diff = outer(x[pospos],x[pospos], FUN = "-")
  y_diff = outer(y[pospos],y[pospos], FUN = "-")
  w_diff = outer(w[pospos],w[pospos])
  w_diff[lower.tri(w_diff, diag = TRUE)] <- 0

  p_conc = sum((x_diff*y_diff*upper.tri(x_diff) > 0)*w_diff)/sum(w_diff)
  p_disc = sum((x_diff*y_diff*upper.tri(x_diff) < 0)*w_diff)/sum(w_diff)

  tau_11 = p_conc - p_disc

  p_11 = sum(w*pospos)/sum(w)

  p_00 = sum(w*(!posx & !posy))/sum(w)

  p_01 = sum(w*(!posx & posy))/sum(w)

  p_10 = sum(w*(posx & !posy))/sum(w)

  # define p_1 and p_2
  # p_1 (weighted) proportion of times the x-values with y = 0 is greater than
  # the x-values with y > 0
  x_10 = x[!posy]
  x_11 = x[posy]

  w_x_10 = w[!posy]
  w_x_11 = w[posy]

  w_x_outer = outer(w_x_10, w_x_11)

  if (sum(w_x_outer) == 0) {
    p_1 <- 0
  } else {
    p_1 = sum(outer(x_10, x_11, FUN = ">")*w_x_outer) / sum(w_x_outer)
  }

  # analogous for y
  y_10 = y[!posx]
  y_11 = y[posx]

  w_y_10 = w[!posx]
  w_y_11 = w[posx]

  w_y_outer = outer(w_y_10, w_y_11)

  if (sum(w_y_outer) == 0) {
    p_2 <- 0
  } else {
    p_2 = sum(outer(y_10, y_11, FUN = ">")*w_y_outer) / sum(w_y_outer)
  }

  tauStar = (p_11^2)*tau_11 +
    2*(p_00*p_11 - p_01*p_10) +
    2*p_11*(p_10*(1 - 2*p_1) + p_01*(1 - 2*p_2))

  return(tauStar)
}



##############################################

#' the weightedPearson function
#'
#' @title weightedPearson
#' @param w weight vector, values should be between 0 and 1
#' @param x x and y are data vectors
#' @param y x and y are data vectors
#' @return \code{numeric} weighted correlation value between x and y

#' @examples
#'
#' x = rnorm(100)
#' y = rnorm(100)
#' w = runif(100)
#' weightedPearson(x,y,w)
#'
#' @export

weightedPearson = function(x, y, w = 1) {

  if (length(x) != length(y)) stop("data must be the same length")

  if (length(w) == 1) {
    w <- rep(w, length(x))
  }

  nw = sum(w)
  wssx = nw * sum(w * (x^2)) - sum(w * x)^2
  wssy = nw * sum(w * (y^2)) - sum(w * y)^2
  wssxy = nw * sum(w * x * y) - sum(w * x) * sum(w * y)
  wcor = wssxy/sqrt(wssx * wssy)
  return(wcor)
}



##############################################

#' the weightedSpearman function
#'
#' @title weightedSpearman
#' @param w weight vector, values should be between 0 and 1
#' @param x x and y are data vectors
#' @param y x and y are data vectors
#' @return \code{numeric} weighted correlation value between x and y

#' @examples
#'
#' x = rnorm(100)
#' y = rnorm(100)
#' w = runif(100)
#' weightedSpearman(x,y,w)
#'
#' @export

weightedSpearman = function(x, y, w = 1) {

  if (length(x) != length(y)) {
    stop("x and y should have the same length")
  }
  if (length(w) == 1) {
    w <- rep(w, length(x))
  }

  keep = w > 0

  xr = rank(x[keep])
  yr = rank(y[keep])
  return(weightedPearson(x = xr, y = yr, w =w[keep]))
}






##############################################

#' the weightedVariance function
#'
#' @title weightedVariance
#' @param w weight vector, values should be between 0 and 1
#' @param x x is a data vector
#' @param y default to NULL, if given it is ignored
#' @return \code{numeric} weighted variance value for x

#' @examples
#'
#' x = rnorm(100)
#' w = runif(100)
#' weightedVariance(x,w = w)
#'
#' @export

weightedVariance = function(x, y = NULL, w) {
  # args x,w (if y given, it is ignored)

  w <- w/sum(w)
  nw <- sum(w)
  n <- sum(w != 0)
  wssx <- nw * sum(w * (x^2)) - sum(w * x)^2
  wvar <- (n/(n - 1))*wssx
  return(wvar)

}




##############################################

#' the weightedZISpearman function calculates weighted rho\*,
#' where rho\* is described in Pimentel et al (2009).
#' This association measure is defined for zero-inflated,
#' non-negative random variables.
#'
#' @title weightedZISpearman
#' @param w weight vector, values should be between 0 and 1
#' @param x x and y are non-negative data vectors
#' @param y x and y are non-negative data vectors
#' @return \code{numeric} weighted rho* association value between x and y
#'
#'  Pimentel, Ronald Silva, "Kendall's Tau and Spearman's Rho for
#'  Zero-Inflated Data" (2009). Dissertations. 721.
#'  https://scholarworks.wmich.edu/dissertations/721

#' @examples
#'
#' x = pmax(0,rnorm(100))
#' y = pmax(0,rnorm(100))
#' w = runif(100)
#' weightedZISpearman(x,y,w)
#'
#' @export

weightedZISpearman <- function(x, y, w = 1) {

  # needs the original values, not the ranks

  if (any(x < 0 | y < 0)) {
    stop("x and/or y values have negative values")
  }
  if (length(x) != length(y)) {
    stop("x and y should have the same length")
  }
  if (length(w) == 1) {
    w <- rep(w, length(x))
  }

  posx = x > 0
  posy = y > 0
  pospos = posx & posy

  p_11 = sum(w * pospos)/sum(w)
  p_00 = sum(w * (!posx & !posy))/sum(w)
  p_01 = sum(w * (!posx & posy))/sum(w)
  p_10 = sum(w * (posx & !posy))/sum(w)
  rho_11 = weightedSpearman(x = x[pospos], y = y[pospos], w = w[pospos])
  rho_star = p_11 * (p_01 + p_11) * (p_10 + p_11) * rho_11 +
    3*(p_00 * p_11 - p_10 * p_01)

  if (is.na(rho_star)) {
    print("Zero inflated Spearman correlation is undefined,
          returning Spearman correlation")
    rho = weightedSpearman(x = x, y = y, w = w)
    return(rho)
  }

  return(rho_star)
}

