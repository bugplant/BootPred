#' Generate random deviates from beta-binomial distribution
#'
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param size number of trials (zero or more).
#' @param mu mean, defined by shape1/(shape1 + shape2) in beta distribution
#' @param phi dispersion patameter, defined by shape1 + shape2 in beta distribution
#'
#' @importFrom stats rbinom rbeta
#'
#' @examples rbetabinom(10, 5, 1/3, 7)
#' @examples rbetabinom(c(10,13), c(4,9), c(1/3, 1/5), c(7, 15))
#'
#' @export
#'
rbetabinom <- function(n, size, mu, phi){
  rbinom(n = n, size = size, prob = rbeta(n = n, shape1 = mu*phi, shape2 = (1-mu)*phi))
}

#' Give the density of beta-binomial distribution
#'
#' @param x vector of quantiles.
#' @param size number of trials (zero or more).
#' @param mu mean, defined by shape1/(shape1 + shape2) in beta distribution
#' @param phi dispersion patameter, defined by shape1 + shape2 in beta distribution
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @examples dbetabinom(0:10, 10, 1/3, 7)
#'
#' @export
#'
dbetabinom <- function(x, size, mu, phi, log = F){
  f <- function(x, size, mu, phi){
    choose(size, x)*(gamma(mu*phi + x)/gamma(mu*phi))*(gamma((1-mu)*phi + size - x)/gamma(((1-mu)*phi)))*(gamma(phi)/gamma(phi+size))
  }
  d <- mapply(f, x = x, size = size, mu = mu, phi = phi)

  if(log == T) d <- log(d)

  return(d)
}

#' Give the distribution function of beta-binomial distribution
#'
#' @param q vector of quantiles.
#' @param size number of trials (zero or more).
#' @param mu mean, defined by shape1/(shape1 + shape2) in beta distribution
#' @param phi dispersion patameter, defined by shape1 + shape2 in beta distribution
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ≤ x], otherwise, P[X > x].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @examples pbetabinom(5, 10, 1/3, 7)
#'
#' @export
#'
pbetabinom <- function(q, size, mu, phi, lower.tail = T, log.p = F){
  f <- function(q, size, mu, phi){
    tmp <- sum(dbetabinom(x = 0:floor(q), size, mu, phi))
    tmp[q == size] <- 1
    return(tmp)
  }
  p <- mapply(f, q = q, size = size, mu = mu, phi = phi)
  if(lower.tail == F) p <- 1 - p
  if(log.p == T) p <- log(p)

  return(p)
}

#' Give the distribution function of beta-binomial distribution
#'
#' @param p vector of probabilities.
#' @param size number of trials (zero or more).
#' @param mu mean, defined by shape1/(shape1 + shape2) in beta distribution
#' @param phi dispersion patameter, defined by shape1 + shape2 in beta distribution
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ≤ x], otherwise, P[X > x].
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @examples qbetabinom(0.2, 10, 1/3, 7)
#'
#' @export
#'
qbetabinom <- function(p, size, mu, phi, lower.tail = T, log.p = F){
  if(log.p == T) p <- exp(p)

  f <- function(p, size, mu, phi, lower.tail){
    tmp <- pbetabinom(0:size, size, mu, phi, lower.tail)
    if(lower.tail == T){
      sum(tmp <= p)
    } else {
      size + 1 - sum(tmp <= p)
    }
  }
  q <- mapply(f, p = p, size = size, mu = mu, phi = phi, lower.tail = lower.tail)
  q[q > size] <- size
  return(q)
}

#' Generate random deviates from hurdle Gamma distribution
#'
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param shape shape parameters. Must be positive.
#' @param rate  rate parameters. Must be positive.
#' @param ziprob probability of zero(0 < ziprob < 1)
#'
#' @importFrom stats rgamma
#'
#' @examples r_hgamma(10, 5, 1, 1/3)
#'
#' @export
#'
r_hgamma <- function(n, shape, rate = 1, ziprob){
  tmp <- rbinom(n = n, size = 1, prob = ziprob)
  v <- mapply(function(x, shape, rate){
    if(x == 1) {
      rgamma(n = 1, shape, rate)
    } else {0}
  },
  x = tmp, shape = shape, rate = rate
  )
  unlist(v)
}

#' Generate random deviates from hurdle beta distribution
#'
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param shape1 shape parameters. Must be positive.
#' @param shape2 shape parameters. Must be positive.
#' @param ziprob probability of zero(0 < ziprob < 1)
#'
#' @importFrom stats rbeta
#'
#' @examples r_hbeta(10, 5, 3, 1/3)
#'
#' @export
#'
r_hbeta <- function(n, shape1, shape2, ziprob){
  tmp <- rbinom(n = n, size = 1, prob = ziprob)
  v <- mapply(function(x, shape1, shape2){
    if(x == 1) {
      rbeta(n = 1, shape1, shape2)
    } else {0}
  },
  x = tmp, shape1 = shape1, shape2 = shape2
  )
  unlist(v)
}

#' Generate random deviates from hurdle poisson distribution
#'
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param lambda vector of (non-negative) means.
#' @param ziprob probability of zero(0 < ziprob < 1)
#'
#' @importFrom stats rbinom
#' @importFrom actuar rztpois
#'
#' @seealso [actuar] <https://cran.r-project.org/web/packages/actuar/index.html>
#'
#' @examples r_hpois(10, 5, 1/3)
#'
#' @export
#'
r_hpois <- function(n, lambda, ziprob){
  tmp <- rbinom(n = n, size = 1, prob = ziprob)
  v <- mapply(function(x, lambda){
    if(x == 1) {
      rztpois(n = 1, lambda)
    } else {0}
  },
  x = tmp, lambda = lambda
  )
  unlist(v)
}

#' Generate random deviates from zero-inflated poisson distribution
#'
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param lambda vector of (non-negative) means.
#' @param ziprob probability of zero(0 < ziprob < 1) by Bernoulli distribution
#'
#' @importFrom stats rbinom rpois
#'
#' @examples r_zipois(10, 5, 1/3)
#'
#' @export
#'
r_zipois <- function(n, lambda, ziprob){
  tmp <- rbinom(n = n, size = 1, prob = ziprob)
  v <- mapply(function(x, lambda){
    if(x == 1) {
      rpois(n = 1, lambda)
    } else {0}
  },
  x = tmp, lambda = lambda
  )
  unlist(v)
}

#' Generate random deviates from hurdle negative binomial distribution
#'
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success in each trial. 0 < prob <= 1.
#' @param ziprob probability of zero(0 < ziprob < 1)
#'
#' @importFrom stats rbinom
#' @importFrom actuar rztnbinom
#'
#' @seealso [actuar] <https://cran.r-project.org/web/packages/actuar/index.html>
#'
#' @examples r_hnbinom(10, 10, 1/3, 1/5)
#'
#' @export
#'
r_hnbinom <- function(n, size, prob, ziprob){
  tmp <- rbinom(n = n, size = 1, prob = ziprob)
  v <- mapply(function(x, size, prob){
    if(x == 1) {
      rztnbinom(n = 1, size, prob)
    } else {0}
  },
  x = tmp, size = size, prob = prob
  )
  unlist(v)
}

#' Generate random deviates from zero-inflated negative binomial distribution
#'
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success in each trial. 0 < prob <= 1.
#' @param ziprob probability of zero(0 < ziprob < 1) by Bernoulli distribution
#'
#' @import stats
#'
#' @examples r_zinbinom(10, 10, 1/3, 1/5)
#'
#' @export
#'
r_zinbinom <- function(n, size, prob, ziprob){
  tmp <- rbinom(n = n, size = 1, prob = ziprob)
  v <- mapply(function(x, size, prob){
    if(x == 1) {
      rnbinom(n = 1, size, prob)
    } else {0}
  },
  x = tmp, size = size, prob = prob
  )
  unlist(v)
}

#' Model predictions by bootstrap
#'
#' @param object a model object for which prediction is desired. "lm", "glm" and "glmmTMB" are supported.
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param interval Type of interval calculation. "prediction" and "confidence" are supported.
#' @param level Tolerance/confidence level.
#' @param bootstrap.n number of simulation
#'
#' @examples x <- rnorm(15)
#' @examples y <- x + rnorm(15)
#' @examples new <- data.frame(x = seq(-3, 3, 0.5))
#' @examples predict(lm(y ~ x), new, interval = "prediction")
#'
#' @export
#'
bootstrap_predict <- function(object, newdata, interval = NULL,
                              level = 0.95, bootstrap.n = 10000){
  if(!(class(object)[1] %in% c("lm", "glm", "glmmTMB")))
    stop(sprintf("No methods for %s", class(object)[1]))

  if(class(object)[1] == "lm") {
    result <- bootstrap_predict.lm(object = object,
                                   newdata = newdata,
                                   interval = interval,
                                   level = level,
                                   bootstrap.n = bootstrap.n)
    return(result)
  }
  if(class(object)[1] == "glm") {
    result <- bootstrap_predict.glm(object = object,
                                    newdata = newdata,
                                    interval = interval,
                                    level = level,
                                    bootstrap.n = bootstrap.n)
    return(result)
  }

  if(class(object)[1] == "glmmTMB") {
    result <- bootstrap_predict.glmmTMB(object = object,
                                        newdata = newdata,
                                        interval = interval,
                                        level = level,
                                        bootstrap.n = bootstrap.n)
    return(result)
  }
}

#' Model predictions by bootstrap for lm object
#'
#' @param object a model object for which prediction is desired.
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param interval Type of interval calculation. "prediction" and "confidence" are supported.
#' @param level Tolerance/confidence level.
#' @param bootstrap.n number of simulation
#'
#' @importFrom stats coef formula model.matrix quantile rbeta rbinom rgamma rnbinom rnorm rpois vcov
#' @importFrom MASS mvrnorm
#'
bootstrap_predict.lm <- function(object, newdata, interval = NULL,
                                 level = 0.95, bootstrap.n = 10000){
  crit <- (1 - level)/2
  beta_hat <- coef(object)
  vcov_hat <- vcov(object)
  sigma_hat <- summary(object)$sigma
  beta_sim <- mvrnorm(bootstrap.n, beta_hat, vcov_hat)
  X <- model.matrix(formula(object)[-2], data = newdata)
  m <- beta_hat %*% t(X)
  m_sim <- beta_sim %*% t(X)

  if(interval == "prediction"){
    y_sim <- apply(m_sim, 1, function(m) rnorm(length(m), mean = m, sd = sigma_hat))
    PI <- t(apply(y_sim, 1, quantile, probs = c(crit, 1-crit)))
    PI <- cbind(t(m), PI)
    colnames(PI) <- c("fit", "lwr", "upr")
    return(PI)

  }
  if(interval == "confidence"){
    CI <- t(apply(m_sim, 2, quantile, probs = c(crit, 1-crit)))
    CI <- cbind(t(m), CI)
    colnames(CI) <- c("fit", "lwr", "upr")
    return(CI)
  }
  if(!interval %in% c("prediction", "confidence"))
    stop(sprintf("No methods for %s", interval))
}

#' Model predictions by bootstrap for glm object
#'
#' @param object a model object for which prediction is desired.
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param interval Type of interval calculation. "prediction" and "confidence" are supported.
#' @param level Tolerance/confidence level.
#' @param bootstrap.n number of simulation
#'
#' @importFrom stats coef formula model.matrix quantile rbeta rbinom rgamma rnbinom rnorm rpois vcov
#' @importFrom MASS mvrnorm
#'
bootstrap_predict.glm <- function(object, newdata, interval = NULL,
                                  level = 0.95, bootstrap.n = 10000){
  if(!object$family$family %in% c("binomial", "Gamma", "gaussian", "poisson"))
    stop("Unsupported family")
  crit <- (1 - level)/2

  beta_hat <- coef(object)
  vcov_hat <- vcov(object)
  sigma_hat <- sqrt(summary(object)$dispersion)
  shape_hat <- 1/(sigma_hat)^2
  size_hat <- round(1/mean(1/rowSums(object$model[,1, drop = F])))

  beta_sim <- mvrnorm(bootstrap.n, beta_hat, vcov_hat)
  X <- model.matrix(formula(object)[-2], data = newdata)
  m <- object$family$linkinv(beta_hat %*% t(X))
  m_sim <- object$family$linkinv(beta_sim %*% t(X))
  if(interval == "prediction"){
    if(object$family$family == "binomial") {
      y_sim <- apply(m_sim, 1, function(m) rbinom(length(m), size = size_hat,
                                                  prob = m)/size_hat)
    }

    if(object$family$family == "Gamma") {
      y_sim <- apply(m_sim, 1, function(m) rgamma(length(m),
                                                  shape = shape_hat,
                                                  rate = shape_hat/m))
    }

    if(object$family$family == "gaussian") {
      y_sim <- apply(m_sim, 1, function(m) rnorm(length(m), mean = m, sd = sigma_hat))
    }

    if(object$family$family == "poisson") {
      y_sim <- apply(m_sim, 1, function(m) rpois(length(m), lambda = m))
    }

    PI <- t(apply(y_sim, 1, quantile, probs = c(crit, 1-crit)))
    PI <- cbind(t(m), PI)
    colnames(PI) <- c("fit", "lwr", "upr")
    return(PI)

  } else {

    CI <- t(apply(m_sim, 2, quantile, probs = c(crit, 1-crit)))
    CI <- cbind(t(m), CI)
    colnames(CI) <- c("fit", "lwr", "upr")
    return(CI)

  }

  if(!interval %in% c("prediction", "confidence"))
    stop(sprintf("No methods for %s", interval))
}

#' Model predictions by bootstrap for glmmTMB object
#'
#' @param object a model object for which prediction is desired.
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param interval Type of interval calculation. "prediction" and "confidence" are supported.
#' @param level Tolerance/confidence level.
#' @param bootstrap.n number of simulation
#'
#' @importFrom stats coef formula model.matrix quantile rbeta rbinom rgamma rnbinom rnorm rpois vcov
#' @importFrom MASS mvrnorm
#' @import glmmTMB
#'
bootstrap_predict.glmmTMB <- function(object, newdata, interval = NULL,
                                      level = 0.95, bootstrap.n = 10000){
  if(!object$modelInfo$family$family %in% c("Gamma", "poisson", "truncated_poisson",
                                            "nbinom1", "nbinom2", "truncated_nbinom2",
                                            "binomial", "beta", "betabinomial"))
    stop("Unsupported family")
  crit <- (1 - level)/2
  zi <- length(fixef(object)$zi) > 0

  beta_hat <- fixef(object)$cond
  vcov_hat <- vcov(object)$cond

  disp_hat <- summary(object)$sigma
  shape_hat <- 1/((disp_hat)^2)
  size_hat <- round(1/mean(1/rowSums(object$frame[,1, drop = F])))

  beta_sim <- mvrnorm(bootstrap.n, beta_hat, vcov_hat)

  X <- model.matrix(formula(object)[-2], data = newdata)
  m_sim <- object$modelInfo$family$linkinv(beta_sim %*% t(X))

  if(zi){
    zibeta_hat <- fixef(object)$zi
    zivcov_hat <- vcov(object)$zi
    zibeta_sim <- mvrnorm(bootstrap.n, zibeta_hat, zivcov_hat)
    ziX <- model.matrix(object$call$ziformula, data = newdata)
    zim_sim <- 1/(1+exp(zibeta_sim %*% t(ziX)))
  }

  m <- object$modelInfo$family$linkinv(beta_hat %*% t(X))
  if(zi) m <- m * 1/(1+exp(zibeta_hat %*% t(ziX)))

  if(interval == "prediction"){
    if(object$modelInfo$family$family == "Gamma") {
      if(zi){
        y_sim <- r_hgamma(length(as.vector(t(m_sim))),
                          shape = shape_hat,
                          rate = shape_hat/as.vector(t(m_sim)),
                          ziprob = as.vector(t(zim_sim)))
        y_sim <- matrix(y_sim, ncol = bootstrap.n)
      } else {
        y_sim <- apply(m_sim, 1, function(m) rgamma(length(m),
                                                    shape = shape_hat,
                                                    rate = shape_hat/m))
      }

    }

    if(object$modelInfo$family$family == "poisson") {
      if(zi){
        y_sim <- r_zipois(length(as.vector(t(m_sim))),
                          lambda = as.vector(t(m_sim)),
                          ziprob = as.vector(t(zim_sim)))
        y_sim <- matrix(y_sim, ncol = bootstrap.n)
      } else {
        y_sim <- apply(m_sim, 1, function(m) rpois(length(m), lambda = m))
      }
    }

    if(object$modelInfo$family$family == "truncated_poisson") {
      y_sim <- r_hpois(length(as.vector(t(m_sim))),
                       lambda = as.vector(t(m_sim)),
                       ziprob = as.vector(t(zim_sim)))
      y_sim <- matrix(y_sim, ncol = bootstrap.n)
    }

    if(object$modelInfo$family$family == "nbinom1") {
      if(zi){
        y_sim <- r_zinbinom(length(as.vector(t(m_sim))),
                            size = m_sim/disp_hat,
                            prob = 1/(1+disp_hat),
                            ziprob = as.vector(t(zim_sim)))
        y_sim <- matrix(y_sim, ncol = bootstrap.n)
      } else {
        y_sim <- apply(m_sim, 1, function(m) rnbinom(length(m), size = m/disp_hat,
                                                     prob = 1/(1+disp_hat)))
      }
    }

    if(object$modelInfo$family$family == "nbinom2") {
      if(zi){
        y_sim <- r_zinbinom(length(as.vector(t(m_sim))),
                            size = disp_hat,
                            prob = disp_hat/(disp_hat + as.vector(t(m_sim))),
                            ziprob = as.vector(t(zim_sim)))
        y_sim <- matrix(y_sim, ncol = bootstrap.n)
      } else {
        y_sim <- apply(m_sim, 1, function(m) rnbinom(length(m), size = disp_hat,
                                                     prob = disp_hat/(disp_hat + m)))
      }
    }

    if(object$modelInfo$family$family == "truncated_nbinom2") {
      y_sim <- r_hnbinom(length(as.vector(t(m_sim))),
                         size = disp_hat,
                         prob = disp_hat/(disp_hat + as.vector(t(m_sim))),
                         ziprob = as.vector(t(zim_sim)))
      y_sim <- matrix(y_sim, ncol = bootstrap.n)
    }


    if(object$modelInfo$family$family == "binomial") {
      y_sim <- apply(m_sim, 1, function(m) rbinom(length(m), size = size_hat,
                                                  prob = m)/size_hat)
    }

    if(object$modelInfo$family$family == "beta") {
      if(zi){
        y_sim <- r_hbeta(length(as.vector(t(m_sim))),
                         shape1 = disp_hat*as.vector(t(m_sim)),
                         shape2 = disp_hat*(1-as.vector(t(m_sim))),
                         ziprob = as.vector(t(zim_sim)))
        y_sim <- matrix(y_sim, ncol = bootstrap.n)
      } else {
        y_sim <- apply(m_sim, 1, function(m) rbeta(length(m),
                                                   shape1 = disp_hat*m,
                                                   shape2 = disp_hat*(1-m)))
      }

    }

    if(object$modelInfo$family$family == "betabinomial") {
      y_sim <- apply(m_sim, 1, function(m) rbetabinom(length(m), size = size_hat,
                                                      mu = m,
                                                      phi = disp_hat))/size_hat
    }

    PI <- t(apply(y_sim, 1, quantile, probs = c(crit, 1-crit)))
    PI <- cbind(t(m), PI)
    colnames(PI) <- c("fit", "lwr", "upr")
    return(PI)

  }else{

    if(zi) m_sim <- m_sim * zim_sim

    CI <- t(apply(m_sim, 2, quantile, probs = c(crit, 1-crit)))
    CI <- cbind(t(m), CI)
    colnames(CI) <- c("fit", "lwr", "upr")
    return(CI)
  }

  if(!interval %in% c("prediction", "confidence"))
    stop(sprintf("No methods for %s", interval))
}
