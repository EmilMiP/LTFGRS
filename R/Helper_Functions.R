utils::globalVariables("cir")
utils::globalVariables("thresh")
utils::globalVariables("K_i")
utils::globalVariables("K_pop")

#' Positive definite matrices
#'
#' \code{correct_positive_definite} verifies that a given covariance matrix
#' is indeed positive definite by checking that all eigenvalues are positive.
#' If the given covariance matrix is not positive definite,
#' \code{correct_positive_definite} tries to modify the underlying correlation matrices
#' genetic_corrmat and full_corrmat in order to obtain a positive definite
#' covariance matrix.
#'
#' This function can be used to verify that a given covariance matrix
#' is positive definite. It calculates all eigenvalues in order to
#' investigate whether they are all positive. This property is necessary
#' for the covariance matrix to be used as a Gaussian covariance matrix.
#' It is especially useful to check whether any covariance matrix obtained
#' by \code{\link{construct_covmat_multi}} is positive definite.
#' If the given covariance matrix is not positive definite, \code{correct_positive_definite}
#' tries to modify the underlying correlation matrices (called \code{genetic_corrmat} and
#' \code{full_corrmat} in \code{\link{construct_covmat}} or \code{\link{construct_covmat_multi}}) by
#' multiplying all off-diagonal entries in the correlation matrices by a given number.
#'
#' @param covmat A symmetric and numeric matrix. If the covariance matrix
#' should be corrected, it must have a number of attributes, such as
#' \code{attr(covmat,"fam_vec")}, \code{attr(covmat,"n_fam")},
#' \code{attr(covmat,"add_ind")}, \code{attr(covmat,"h2")},
#' \code{attr(covmat,"genetic_corrmat")}, \code{attr(covmat,"full_corrmat")}
#' and \code{attr(covmat,"phenotype_names")}. Any covariance matrix
#' obtained by \code{\link{construct_covmat}}, \code{\link{construct_covmat_single}}
#' or \code{\link{construct_covmat_multi}} will have these attributes by default.
#' @param correction_val A positive number representing the amount by which
#' \code{genetic_corrmat} and \code{full_corrmat} will be changed, if some
#' eigenvalues are non-positive. That is, correction_val is the number that will be
#' multiplied to all off_diagonal entries in \code{genetic_corrmat} and \code{full_corrmat}.
#' Defaults to 0.99.
#' @param correction_limit A positive integer representing the upper limit for the correction
#' procedure. Defaults to 100.
#'
#' @return If \code{covmat} is a symmetric and numeric matrix and all eigenvalues are
#' positive, \code{correct_positive_definite} simply returns \code{covmat}. If some
#' eigenvalues are not positive and \code{correction_val} is a positive number,
#' \code{correct_positive_definite} tries to convert \code{covmat} into a positive definite
#' matrix. If \code{covmat} has attributes \code{add_ind}, \code{h2},
#' \code{genetic_corrmat}, \code{full_corrmat} and \code{phenotype_names},
#' \code{correct_positive_definite} computes a new covariance matrix using slightly
#' modified correlation matrices \code{genetic_corrmat} and \code{full_corrmat}.
#' If the correction is performed successfully, i.e. if the new covariance matrix
#' is positive definite,the new covariance matrix is returned.
#' Otherwise, \code{correct_positive_definite} returns the original covariance matrix.
#'
#' @examples
#' ntrait <- 2
#' genetic_corrmat <- matrix(0.6, ncol = ntrait, nrow = ntrait)
#' diag(genetic_corrmat) <- 1
#' full_corrmat <- matrix(-0.25, ncol = ntrait, nrow = ntrait)
#' diag(full_corrmat) <- 1
#' h2_vec <- rep(0.6, ntrait)
#' cov <- construct_covmat(fam_vec = c("m", "f"),
#'   genetic_corrmat = genetic_corrmat,
#'   h2 = h2_vec,
#'   full_corrmat = full_corrmat)
#' cov
#' correct_positive_definite(cov)
#'
#' @seealso \code{\link{construct_covmat}}, \code{\link{construct_covmat_single}} and
#' \code{\link{construct_covmat_multi}}.
#'
#' @export
correct_positive_definite = function(covmat, correction_val = .99, correction_limit = 100) {

  # Checking that covmat is symmetric
  if (!isSymmetric.matrix(covmat)) stop("The covariance matrix covmat must be symmetric!")
  # and numeric
  if (!is.numeric(covmat)) stop("The covariance matrix covmat must be numeric!")

  # Checking whether all eigenvalues are positive
  if (any(eigen(covmat)$values < 0)) {

    message("The specified covariance matrix is not positive definite. \n")
  }else{

    return(covmat)
  }

  # If some eigenvalues are negative, correction_val must be specified,
  # it must be numeric and positive in order to correct the covariance matrix.
  if (!is.numeric(correction_val) && !is.integer(correction_val)) stop("correction_val must be numeric!")
  if (correction_val <= 0) stop("correction_val must be positive!")

  # In addition, covmat must have several attributes holding the
  # family members (fam_vec or n_fam), a logical add_ind as well as
  # a numeric value or numeric matrix h2 and numeric matrices
  # genetic_corrmat and full_corrmat in order to change
  # the correlation matrix.
  if (is.null(attr(covmat,"add_ind")) || is.null(attr(covmat,"h2")) ||
     is.null(attr(covmat,"genetic_corrmat")) || is.null(attr(covmat,"full_corrmat"))) {

    warning("The required attributes are missing... The covariance matrix could not be corrected!")
    return(covmat)
  }

  # If the covariance matrix is for a single phenotype, it is not
  # possible to correct the covariance matrix.
  if (length(attr(covmat,"h2")) == 1) {
    warning("The covariance matrix cannot be corrected...")
    return(covmat)
  }

  # Furthermore, correction_limit must be a positive number
  if (!is.numeric(correction_limit)) stop("correction_limit must be numeric!")
  if (correction_limit <= 0) stop("Correction limit must be positive!")

  message("Trying to correct the covariance matrix...\n")

  # The covariance matrix will be modified at most correction_limit times.
  n <- 0
  # Storing the old covariance matrix
  old_covmat <- covmat

  # We also extract the vectors holding the family members
  fam_vec <- setdiff(attr(covmat,"fam_vec"), c("g","o"))
  #n_fam <- attr(covmat,"n_fam")[stringr::str_detect(names(attr(covmat,"n_fam")), "^[^go]")]

  while (any(eigen(covmat)$values < 0) && n <= correction_limit) {

    # Changing the correlation matrices slightly by
    # multiplying all entries by correction_val.
    genetic_corrmat <- attr(covmat,"genetic_corrmat")*correction_val
    diag(genetic_corrmat) <- 1
    full_corrmat <- attr(covmat,"full_corrmat")*correction_val
    diag(full_corrmat) <- 1
    # Computing a new covariance matrix
    covmat <- construct_covmat(fam_vec = fam_vec,
                               n_fam = NULL,
                               add_ind = attr(covmat,"add_ind"),
                               genetic_corrmat = genetic_corrmat,
                               full_corrmat = full_corrmat,
                               h2 = attr(covmat,"h2"),
                               phen_names = attr(covmat,"phenotype_names"))
    # Updating
    n <- n + 1
  }

  # If the matrix has been modified correction_limit times,
  # the correction step is aborted
  if (n > correction_limit) {
    message("The covariance matrix could not be corrected. Consider revisiting it.\n")
    return(old_covmat)
  }

  message(paste0("The correction was performed successfully! All off-diagonal entries are corrected by", round(correction_val^n,digits = 3),".\n"))

  return(covmat)
}

#' A simplied version of correct_positive_definite
#'
#' multiplies off-diagonal elements of a matrix with a correction value until all eigen values are positive or correction limit is reached.
#'
#' @param covmat Covariance matrix that needs to be positive definite
#' @param correction_limit The maximum number of times to correct off-diagonal elements
#' @param correction_val The value to multiply on the off-diagnoal elements.
#'
#' @return Returns list with the corrected covmat and number of iterations used to achieve positive definite.
#' @noRd

correct_positive_definite_simplified = function(covmat, correction_limit = 100, correction_val = 0.99) {
  eigen_val = eigen(covmat)

  if ( any(eigen_val$values < 0) ) {
    og_diag = diag(covmat)
    n = 0

    while (any(eigen(covmat)$values < 0) & n <= correction_limit ) {
      covmat = covmat*correction_val
      diag(covmat) = og_diag
      n = n + 1
    }

    if (any(eigen_val$values < 0)) stop("Unable to enforce a positive definite covariance matrix.")
    return(list(covmat = covmat, nitr = n))
  } else {
    return(list(covmat = covmat, nitr = 0))
  }
}



#' CDF for truncated normal distribution.
#'
#' \code{truncated_normal_cdf} computes the cumulative density
#' function for a truncated normal distribution.
#'
#' This function can be used to compute the value of the cumulative
#' density function for a truncated normal distribution given an
#' individual's true underlying liability.
#'
#' @param liability  A number representing the individual's
#' true underlying liability.
#' @param lower A number representing the lower cutoff point for the
#' truncated normal distribution. Defaults to 1.645
#' (stats::qnorm(0.05, lower.tail = FALSE)).
#' @param upper A number representing the upper cutoff point of the
#' truncated normal distribution. Must be greater or equal to lower.
#' Defaults to Inf.
#'
#' @return If liability is a number and the lower and upper cutoff points
#' are numbers satisfying lower <= upper, then \code{truncated_normal_cdf}
#' returns the probability that the liability will take on a value less than
#' or equal to \code{liability}.
#'
#' @examples
#' curve(sapply(liability, truncated_normal_cdf), from = qnorm(0.05, lower.tail = FALSE), to = 3.5,
#'  xname = "liability")
#'
#' @export
truncated_normal_cdf = function(liability, lower = stats::qnorm(0.05, lower.tail = FALSE), upper = Inf) {

  # Checking that the liability is valid
  if (!is.numeric(liability) && !is.integer(liability)) stop("The liability must be numeric!")
  # Checking that the lower and upper cutoff points are valid
  if (!is.numeric(lower) && !is.integer(lower)) stop("The lower cutoff point must be numeric!")
  if (!is.numeric(upper) && !is.integer(upper)) stop("The upper cutoff point must be numeric!")
  if (upper < lower) {
    warning("The upper cutoff point is below the lower cutoff point! \n
The upper and lower cutoff points will be swapped...")

    temp <- lower
    lower <- upper
    upper <- temp
  }

  return(stats::pnorm(liability) - stats::pnorm(lower)) / (stats::pnorm(upper) - stats::pnorm(lower))
}


#' Convert age to cumulative incidence rate
#'
#' \code{convert_age_to_cir} computes the cumulative incidence
#' rate from a person's age.
#'
#' Given a person's age, \code{convert_age_to_cir} can be used
#' to compute the cumulative incidence rate (cir), which is given
#' by the formula
#' \deqn{pop\_ prev / (1 + exp((mid\_ point - age) * slope))}
#'
#' @param age A non-negative number representing the individual's age.
#' @param pop_prev A positive number representing the overall
#' population prevalence. Must be at most 1. Defaults to 0.1.
#' @param mid_point A positive number representing the mid point
#' logistic function. Defaults to 60.
#' @param slope A number holding the rate of increase.
#' Defaults to 1/8.
#'
#' @return If age and mid_point are positive numbers, if pop_prev
#' is a positive number between 0 and 1 and if slope is a valid number,
#' then \code{convert_age_to_cir} returns a number, which is equal to
#' the cumulative incidence rate.
#'
#' @examples
#' curve(sapply(age, convert_age_to_cir), from = 10, to = 110, xname = "age")
#' @export
convert_age_to_cir = function(age, pop_prev = .1, mid_point = 60, slope = 1/8) {

  # Checking that age is valid
  if (!is.numeric(age) && !is.integer(age)) stop("The age must be numeric!")
  if (age < 0) stop("The age must be non-negative!")
  if (age >= 150) warning("At this point, it is unrealistic to be of age 150 or older...")

  # Checking that pop_prev is valid
  if (!is.numeric(pop_prev) && !is.integer(pop_prev)) stop("The population prevalence pop_prev must be numeric!")
  if (pop_prev <= 0) stop("The population prevalence pop_prev must be positive!")
  if (pop_prev > 1) stop("The population prevalence pop_prev must be smaller or equal to 1!")

  # Checking that mid_point is valid
  if (!is.numeric(mid_point) && !is.integer(mid_point)) stop("The mid point mid_point must be numeric!")
  if (mid_point <= 0) stop("The mid point mid_point must be positive!")

  # Checking that slope is valid
  if (!is.numeric(slope) && !is.integer(slope)) stop("The slope must be numeric!")

  return(pop_prev / (1 + exp((mid_point - age) * slope)))
}

#' Convert age to threshold
#'
#' \code{convert_age_to_thresh} computes the threshold
#' from a person's age using either the logistic function
#' or the truncated normal distribution
#'
#' Given a person's age, \code{convert_age_to_thresh} can be used
#' to first compute the cumulative incidence rate (cir), which is
#' then used to compute the threshold using either the
#' logistic function or the truncated normal distribution.
#' Under the logistic function, the formula used to compute
#' the threshold from an individual's age is given by
#' \deqn{qnorm(pop\_ prev / (1 + exp((mid\_ point - age) * slope)), lower.tail = F)},
#' while it is given by
#' \deqn{qnorm((1 - (age-min\_ age)/max\_ age) * (pnorm(upper) - pnorm(lower)) + pnorm(lower))}
#' under the truncated normal distribution.
#'
#' @param age A non-negative number representing the individual's age.
#' @param dist A string indicating which distribution to use.
#' If dist = "logistic", the logistic function will be used to
#' compute the age of onset.
#' If dist = "normal", the truncated normal distribution will be used instead.
#' Defaults to "logistic".
#' @param pop_prev Only necessary if dist = "logistic". A positive number representing the overall
#' population prevalence. Must be at most 1. Defaults to 0.1.
#' @param mid_point Only necessary if dist = "logistic". A positive number representing the mid point
#' logistic function. Defaults to 60.
#' @param slope Only necessary if dist = "logistic". A number holding the rate of increase.
#' Defaults to 1/8.
#' @param min_age Only necessary if dist = "normal". A positive number representing the individual's earliest age.
#' Defaults to 10.
#' @param max_age Only necessary if dist = "normal". A positive number representing the individual's latest age.
#' Must be greater than min_aoo. Defaults to 90.
#' @param lower Only necessary if dist = "normal". A number representing the lower cutoff point for the
#' truncated normal distribution. Defaults to 1.645
#' (stats::qnorm(0.05, lower.tail = FALSE)).
#' @param upper Only necessary if dist = "normal". A number representing the upper cutoff point of the
#' truncated normal distribution. Must be greater or equal to lower.
#' Defaults to Inf.
#'
#' @return If age is a positive number and all other necessary arguments are valid,
#' then \code{convert_age_to_thresh} returns a number, which is equal to
#' the threshold.
#'
#' @examples
#' curve(sapply(age, convert_age_to_thresh), from = 10, to = 110, xname = "age")
#' @export
convert_age_to_thresh = function(age, dist = "logistic", pop_prev = .1, mid_point = 60, slope = 1/8,
                                 min_age = 10, max_age = 90, lower = stats::qnorm(0.05, lower.tail = FALSE), upper = Inf) {
  # Checking that age is valid
  if (!is.numeric(age) && !is.integer(age)) stop("The age must be numeric!")
  if (age < 0) stop("The age must be non-negative!")

  # Checking that dist is either logistic or normal.
  if (is.character(dist)) {

    valid_is <- grep("logistic|normal", dist)

  }else{
    stop("dist must be a string!")
  }

  # Checking whether dist is empty or of length >1.
  if (length(valid_is) == 0) {

    warning("dist is not of the required format! \n The function will use the logistic function to compute the age of onset!")
    dist <- "logistic"
  } else {
    if (length(valid_is) > 1) warning("dist is not of the required format! \n The function will use the first valid entry to compute the age of onset!")
    dist <- dist[valid_is[1]]
  }

  # If dist = logistic, the logistic function will be used to compute the age of onset
  if (dist == "logistic") {

    # Checking that pop_prev is valid
    if (!is.numeric(pop_prev) && !is.integer(pop_prev)) stop("The population prevalence pop_prev must be numeric!")
    if (pop_prev <= 0) stop("The population prevalence pop_prev must be positive!")
    if (pop_prev > 1) stop("The population prevalence pop_prev must be smaller or equal to 1!")

    # Checking that mid_point is valid
    if (!is.numeric(mid_point) && !is.integer(mid_point)) stop("The mid point mid_point must be numeric!")
    if (mid_point <= 0) stop("The mid point mid_point must be positive!")

    # Checking that slope is valid
    if (!is.numeric(slope) && !is.integer(slope)) stop("The slope must be numeric!")

    # Computing the threshold
    return(stats::qnorm(pop_prev / (1 + exp((mid_point - age) * slope)), lower.tail = F))

  }

  # if dist = normal, the truncated normal distribution will be used to compute the age of onset
  if (dist == "normal") {

    # Checking that min_age and max_age are valid.
    if (!is.numeric(min_age) && !is.integer(min_age)) stop("The earliest age min_age must be numeric!")
    if (!is.numeric(max_age) && !is.integer(max_age)) stop("The latest age max_age must be numeric!")
    if (min_age <= 0) stop("The earliest age min_age must be positive!")
    if (max_age <= 0) stop("The latest age max_age must be positive!")
    if (min_age > max_age) {
      warning("The latest age max_age is below the earliest age min_age! \n
The earliest and latest age will be swapped...")

      min_age <- min_age + max_age
      max_age <- min_age - max_age
      min_age <- min_age - max_age
    }

    # Checking that the lower and upper cutoff points are valid
    if (!is.numeric(lower) && !is.integer(lower)) stop("The lower cutoff point must be numeric!")
    if (!is.numeric(upper) && !is.integer(upper)) stop("The upper cutoff point must be numeric!")
    if (upper < lower) {
      warning("The upper cutoff point is below the lower cutoff point! \n
The upper and lower cutoff points will be swapped...")

      temp <- lower
      lower <- upper
      upper <- temp
    }
    # Computing the threshold
    return(stats::qnorm((1 - (age - min_age)/max_age) * (stats::pnorm(upper) - stats::pnorm(lower)) + stats::pnorm(lower)))
  }
}


#' Convert cumulative incidence rate to age
#'
#' \code{convert_cir_to_age} computes the age
#' from a person's cumulative incidence rate.
#'
#' Given a person's cumulative incidence rate (cir), \code{convert_cir_to_age}
#' can be used to compute the corresponding age, which is given by
#' \deqn{mid\_ point - \log(pop\_ prev/cir - 1) * 1/slope}
#'
#' @param cir A positive number representing the individual's cumulative
#' incidence rate.
#' @param pop_prev A positive number representing the overall
#' population prevalence. Must be at most 1 and must be larger than
#' cir. Defaults to 0.1.
#' @param mid_point A positive number representing the mid point
#' logistic function. Defaults to 60.
#' @param slope A number holding the rate of increase.
#' Defaults to 1/8.
#'
#' @return If cir and mid_point are positive numbers, if pop_prev
#' is a positive number between 0 and 1 and if slope is a valid number,
#' then \code{convert_cir_to_age} returns a number, which is equal to
#' the current age.
#'
#' @examples
#' curve(sapply(cir, convert_cir_to_age), from = 0.001, to = 0.099, xname = "cir")
#' @export
convert_cir_to_age = function(cir, pop_prev = .1, mid_point = 60, slope = 1/8) {

  # Checking that age is valid
  if (!is.numeric(cir) && !is.integer(cir)) stop("The cumulative incidence rate cir must be numeric!")
  if (cir <= 0) stop("The cumulative incidence rate cir must be positive!")

  # Checking that pop_prev is valid
  if (!is.numeric(pop_prev) && !is.integer(pop_prev)) stop("The population prevalence pop_prev must be numeric!")
  if (pop_prev <= 0) stop("The population prevalence pop_prev must be positive!")
  if (pop_prev > 1) stop("The population prevalence pop_prev must be smaller or equal to 1!")

  # Checking that mid_point is valid
  if (!is.numeric(mid_point) && !is.integer(mid_point)) stop("The mid point mid_point must be numeric!")
  if (mid_point <= 0) stop("The mid point mid_point must be positive!")

  # Checking that slope is valid
  if (!is.numeric(slope) && !is.integer(slope)) stop("The slope must be numeric!")

  if (cir >= pop_prev) {
    return(NA)
  }else{

    res <- mid_point - log(pop_prev/cir - 1) * 1/slope

    if (res > 0) return(res)
    if (res <= 0) return(0)
  }
}


#' Convert liability to age of onset
#'
#' \code{convert_liability_to_aoo} computes the age
#' of onset from an individual's true underlying liability using
#' either the logistic function or the truncated normal distribution.
#'
#' Given a person's cumulative incidence rate (cir), \code{convert_liability_to_aoo}
#' can be used to compute the corresponding age. Under the logistic function,
#' the age is given by
#' \deqn{mid\_ point - log(pop\_ prev/cir - 1) * 1/slope},
#' while it is given by
#' \deqn{(1 - truncated\_ normal\_ cdf(liability = liability, lower = lower , upper = upper)) * max\_ aoo + min\_ aoo}
#' under the truncated normal distribution.
#'
#' @param liability A number representing the individual's
#' true underlying liability.
#' @param dist A string indicating which distribution to use.
#' If dist = "logistic", the logistic function will be used to
#' compute the age of onset.
#' If dist = "normal", the truncated normal distribution will be used instead.
#' Defaults to "logistic".
#' @param pop_prev Only necessary if dist = "logistic". A positive number representing the overall
#' population prevalence. Must be at most 1. Defaults to 0.1.
#' @param mid_point Only necessary if dist = "logistic". A positive number representing the mid point
#' logistic function. Defaults to 60.
#' @param slope Only necessary if dist = "logistic". A number holding the rate of increase.
#' Defaults to 1/8.
#' @param min_aoo Only necessary if dist = "normal". A positive number representing the individual's earliest age of onset.
#' Defaults to 10.
#' @param max_aoo Only necessary if dist = "normal". A positive number representing the individual's latest age of onset.
#' Must be greater than min_aoo. Defaults to 90.
#' @param lower Only necessary if dist = "normal". A number representing the lower cutoff point for the
#' truncated normal distribution. Defaults to 1.645
#' (stats::qnorm(0.05, lower.tail = FALSE)).
#' @param upper Only necessary if dist = "normal". A number representing the upper cutoff point of the
#' truncated normal distribution. Must be greater or equal to lower.
#' Defaults to Inf.
#'
#' @return If liability is a number and all other necessary arguments are valid,
#' then \code{convert_liability_to_aoo} returns a positive number, which is equal to
#' the age of onset.
#'
#' @examples
#' curve(sapply(liability, convert_liability_to_aoo), from = 1.3, to = 3.5, xname = "liability")
#' curve(sapply(liability, convert_liability_to_aoo, dist = "normal"),
#'  from = qnorm(0.05, lower.tail = FALSE), to = 3.5, xname = "liability")
#'
#' @export
convert_liability_to_aoo = function(liability, dist = "logistic", pop_prev = .1, mid_point = 60, slope = 1/8,
                                    min_aoo = 10, max_aoo = 90, lower = stats::qnorm(0.05, lower.tail = FALSE), upper = Inf ) {

  # Checking that liability is valid
  if (!is.numeric(liability) && !is.integer(liability)) stop("The liability must be numeric!")

  # Checking that dist is either logistic or normal.
  if (is.character(dist)) {

    valid_is <- grep("logistic|normal", dist)

  } else {
    stop("dist must be a string!")
  }

  # Checking whether dist is empty or of length >1.
  if (length(valid_is) == 0) {

    warning("dist is not of the required format! \n The function will use the logistic function to compute the age of onset!")
    dist <- "logistic"
  } else {

    if (length(valid_is) > 1) warning("dist is not of the required format! \n The function will use the first valid entry to compute the age of onset!")
    dist <- dist[valid_is[1]]
  }

  # If dist = logistic, the logistic function will be used to compute the age of onset
  if (dist == "logistic") {

    # Checking that pop_prev is valid
    if (!is.numeric(pop_prev) && !is.integer(pop_prev)) stop("The population prevalence pop_prev must be numeric!")
    if (pop_prev <= 0) stop("The population prevalence pop_prev must be positive!")
    if (pop_prev > 1) stop("The population prevalence pop_prev must be smaller or equal to 1!")

    # Checking that mid_point is valid
    if (!is.numeric(mid_point) && !is.integer(mid_point)) stop("The mid point mid_point must be numeric!")
    if (mid_point <= 0) stop("The mid point mid_point must be positive!")

    # Checking that slope is valid
    if (!is.numeric(slope) && !is.integer(slope)) stop("The slope must be numeric!")

    # Computing the age of onset
    if (stats::pnorm(liability, lower.tail = F) >= pop_prev) {
      return(NA)
    } else {

      res <- mid_point - log(pop_prev/stats::pnorm(liability, lower.tail = F) - 1) * 1/slope

      if (res > 0) return(res)
      if (res <= 0) return(0)
    }
  }

  # if dist = normal, the truncated normal distribution will be used to compute the age of onset
  if (dist == "normal") {

    # Checking that min_aoo and max_aoo are valid.
    if (!is.numeric(min_aoo) && !is.integer(min_aoo)) stop("The earliest age of onset min_aoo must be numeric!")
    if (!is.numeric(max_aoo) && !is.integer(max_aoo)) stop("The latest age of onset max_aoo must be numeric!")
    if (min_aoo <= 0) stop("The earliest age of onset min_aoo must be positive!")
    if (max_aoo <= 0) stop("The latest age of onset max_aoo must be positive!")
    if (min_aoo > max_aoo) {
      warning("The latest age of onset max_aoo is below the earliest age of onset min_aoo! \n
The earliest and latest age of onset will be swapped...")

      min_aoo <- min_aoo + max_aoo
      max_aoo <- min_aoo - max_aoo
      min_aoo <- min_aoo - max_aoo
    }

    # Checking that the lower and upper cutoff points are valid
    if (!is.numeric(lower) && !is.integer(lower)) stop("The lower cutoff point must be numeric!")
    if (!is.numeric(upper) && !is.integer(upper)) stop("The upper cutoff point must be numeric!")
    if (upper < lower) {
      warning("The upper cutoff point is below the lower cutoff point! \n
The upper and lower cutoff points will be swapped...")

      temp <- lower
      lower <- upper
      upper <- temp
    }

    # Computing the age of onset
    res <- (1 - truncated_normal_cdf(liability = liability, lower = lower , upper = upper)) * max_aoo + min_aoo

    return(res)
  }
}

#' Convert the heritability on the observed scale to that on the liability scale
#'
#' \code{convert_observed_to_liability_scale} transforms the heritability on the
#' observed scale to the heritability on the liability scale.
#'
#' This function can be used to transform the heritability on the observed
#' scale to that on the liability scale. \code{convert_observed_to_liability_scale}
#' uses either Equation 17 (if prop_cases = NULL) or Equation 23 from
#' Sang Hong Lee, Naomi R. Wray, Michael E. Goddard and Peter M. Visscher, "Estimating
#' Missing Heritability for Diseases from Genome-wide Association Studies",
#' The American Journal of Human Genetics, Volume 88, Issue 3, 2011, pp. 294-305,
#' \doi{10.1016/j.ajhg.2011.02.002} to transform the heritability on the observed
#' scale to the heritability on the liability scale.
#'
#' @param obs_h2 A number or numeric vector representing the liability-scale
#' heritability(ies)on the observed scale. Must be non-negative and at most 1.
#' Defaults to 0.5
#' @param pop_prev A number or numeric vector representing the population prevalence(s). All
#' entries must be non-negative and at most one.
#' If it is a vector, it must have the same length as obs_h2. Defaults to 0.05.
#' @param prop_cases Either NULL or a number or a numeric vector representing the proportion
#' of cases in the sample. All entries must be non-negative and at most one.
#' If it is a vector, it must have the same length as obs_h2. Defaults to 0.5.
#'
#' @return If \code{obs_h2}, \code{pop_prev} and \code{prop_cases} are non-negative numbers
#' that are at most one, the function returns the heritability on the liability
#' scale using Equation 23 from
#' Sang Hong Lee, Naomi R. Wray, Michael E. Goddard and Peter M. Visscher, "Estimating
#' Missing Heritability for Diseases from Genome-wide Association Studies",
#' The American Journal of Human Genetics, Volume 88, Issue 3, 2011, pp. 294-305,
#' \doi{10.1016/j.ajhg.2011.02.002}.
#' If \code{obs_h2}, \code{pop_prev} and \code{prop_cases} are non-negative numeric
#' vectors where all entries are at most one, the function returns a vector of the same
#' length as obs_h2. Each entry holds to the heritability on the liability
#' scale which was obtained from the corresponding entry in obs_h2 using Equation 23.
#' If \code{obs_h2} and \code{pop_prev} are non-negative numbers that are at most
#' one and \code{prop_cases} is \code{NULL}, the function returns the heritability
#' on the liability scale using Equation 17 from
#' Sang Hong Lee, Naomi R. Wray, Michael E. Goddard and Peter M. Visscher, "Estimating
#' Missing Heritability for Diseases from Genome-wide Association Studies",
#' The American Journal of Human Genetics, Volume 88, Issue 3, 2011, pp. 294-305,
#' \doi{10.1016/j.ajhg.2011.02.002}.
#' If \code{obs_h2} and \code{pop_prev} are non-negative numeric vectors such that
#' all entries are at most one, while \code{prop_cases} is \code{NULL},
#' \code{convert_observed_to_liability_scale} returns a vector of the same
#' length as obq_h2. Each entry holds to the liability-scale heritability that
#' was obtained from the corresponding entry in obs_h2 using Equation 17.
#'
#' @examples
#' convert_observed_to_liability_scale()
#' convert_observed_to_liability_scale(prop_cases=NULL)
#' convert_observed_to_liability_scale(obs_h2 = 0.8, pop_prev = 1/44,
#'                                     prop_cases = NULL)
#' convert_observed_to_liability_scale(obs_h2 = c(0.5,0.8),
#'                                     pop_prev = c(0.05, 1/44),
#'                                     prop_cases = NULL)
#'
#' @references
#' Sang Hong Lee, Naomi R. Wray, Michael E. Goddard, Peter M. Visscher (2011, March). Estimating
#' Missing Heritability for Diseases from Genome-wide Association Studies. In The American Journal
#' of Human Genetics (Vol. 88, Issue 3, pp. 294-305). \doi{10.1016/j.ajhg.2011.02.002}
#'
#' @export
convert_observed_to_liability_scale <- function(obs_h2 = 0.5, pop_prev = 0.05, prop_cases = 0.5) {

  # Checking that the observed heritabilities are valid
  if (!is.numeric(obs_h2) && !is.integer(obs_h2)) stop("The observed heritability(ies) must be numeric!")
  if (any(obs_h2 < 0)) stop("The observed heritability(ies) must be non-negative!")
  if (any(obs_h2 > 1)) stop("The observed heritability(ies) must be smaller than or equal to one!")
  # Checking that the population prevalences are valid
  if (!is.numeric(pop_prev) && !is.integer(pop_prev)) stop("The population prevalence(s) must be numeric!")
  if (any(pop_prev < 0)) stop("The population prevalence(s) must be non-negative!")
  if (any(pop_prev > 1)) stop("The population prevalence(s) must be smaller than or equal to one!")

  # Defining the variable z, which is the height of the truncated
  # normal curve at the point t, and where t is the truncated point,
  # such that the fraction of observations larger than t is equal to
  # the population prevalence pop_prev. That is
  # t = qnorm(pop_prev, lower.tail = FALSE)
  z <- stats::dnorm(stats::qnorm(pop_prev, lower.tail = FALSE))

  # Using one of the two possible transformations depending on whether
  # prop_cases is NULL or not.
  if (is.null(prop_cases)) {

    return(obs_h2 * (pop_prev*(1 - pop_prev))/(z^2))
  }else{

    # Checking that the proportions of cases are valid
    if (!is.numeric(prop_cases) && !is.integer(prop_cases)) stop("The proportion(s) of cases must be numeric!")
    if (any(prop_cases < 0)) stop("The proportion(s) of cases must be non-negative!")
    if (any(prop_cases > 1)) stop("The proportion(s) of cases must be smaller than or equal to one!")

    return(obs_h2 * (pop_prev*(1 - pop_prev))/(z^2) * (pop_prev*(1 - pop_prev))/(prop_cases*(1 - prop_cases)))
  }
}



# PA-FGRS related helper functions ----------------------------------------



#' Title: Calculate the mean of the truncated normal distribution
#'
#' @param mu mean value of normal distribution
#' @param sigma standard deviation of normal distribution
#' @param lower lower threshold
#' @param upper upper threshold
#'
#' @returns mean value of the truncated normal distribution
#' @export
#'
#' @importFrom stats dnorm pnorm
#' @examples
#' tnorm_mean()
tnorm_mean = function(mu = 0, sigma = 1, lower = -Inf, upper = Inf) {
  # calculate the mean of the truncated normal distribution
  # using the formula from https://en.wikipedia.org/wiki/Truncated_normal_distribution

  if (upper == -Inf | lower == Inf) stop("tnorm_mean: you may have switched the lower and upper bounds!")

  # returns normal distribution mean of thresholds are -Inf and Inf
  if (lower == -Inf & upper == Inf) {
    return(mu)
  }
  # if lower and upper bounds are the same, we will interpret it as a normal distribution with 0 variance
  # and mean equal to the limits, which can be interpreted as a direc delta function on the limit's values
  if (lower == upper) {
    return(lower)
  }
  # Otherwise, we will use the general formula
  alpha = (lower - mu) / sigma
  beta = (upper - mu) / sigma
  return(mu - sigma * (dnorm(beta) - dnorm(alpha)) / (pnorm(beta) - pnorm(alpha)))
}



#' Title: Calculate the variance of the truncated normal distribution
#'
#' @param mu mean value of normal distribution
#' @param sigma standard deviation of normal distribution
#' @param lower lower threshold
#' @param upper upper threshold
#'
#' @returns mean value of the truncated normal distribution
#' @export
#'
#' @importFrom stats dnorm pnorm
#' @examples
#' tnorm_var()
tnorm_var = function(mu = 0, sigma = 1, lower = -Inf, upper = Inf) {
  # calculate the variance of the truncated normal distribution
  # using the formula from https://en.wikipedia.org/wiki/Truncated_normal_distribution

  if (upper == -Inf | lower == Inf) stop("tnorm_mean: you may have switched the lower and upper bounds!")

  # if the bounds are infinite, the variance of the normal distribution is returned
  if (lower == -Inf & upper == Inf) {
    return(sigma^2)
  }
  # if lower and upper bounds are the same, we will interpret it as a normal
  # distribution with 0 variance (direc distribution on the mean)
  if (lower == upper) {
    return(0)
  }
  # Otherwise, we will use the general formula(s)
  alpha = (lower - mu) / sigma
  beta = (upper - mu) / sigma

  # we will use the special cases to handle infinite bounds, since they may results in NaNs by, e.g. Inf * 0
  if ( beta == Inf ) {
    return(sigma^2 * (1 + (alpha * dnorm(alpha)) / (1 - pnorm(alpha)) - (dnorm(alpha))^2 / (1 - pnorm(alpha))^2))
  }
  if ( alpha == -Inf ) {
    return(sigma^2 * (1 - (beta * dnorm(beta)) / pnorm(beta) - (dnorm(beta))^2 / pnorm(beta)^2))
  }
  # if both alpha and beta are finite we do not experience any NaNs with the general formula
  return(sigma^2 * (1 + (alpha * dnorm(alpha) - beta * dnorm(beta)) / (pnorm(beta) - pnorm(alpha)) - (dnorm(alpha) - dnorm(beta))^2 / (pnorm(beta) - pnorm(alpha))^2))
}




#' Title: Calculates mean and variance of mixture of two truncated normal distributions
#'
#' @param mu Mean value of normal distribution.
#' @param var Variance of normal distribution.
#' @param lower Lower threshold (can be -Inf).
#' @param upper Upper threshold (can be Inf).
#' @param K_i (Stratified) cumulative incidence proportion for the individual.
#' @param K_pop Population prevalence (cumulative incidence proportion).
#'
#' @returns mean and variance of mixture distribution between two truncated normal distributions
#' @export
#'
#' @importFrom dplyr case_when
#' @examples
#' tnorm_mixture_conditional(mu = 0, var = 1, lower = -Inf, upper = Inf, K_i = 0, K_pop = 0.01)
#' tnorm_mixture_conditional(mu = 0, var = 1, lower = -Inf, upper = 2, K_i = .01, K_pop = 0.05)
tnorm_mixture_conditional = function(mu, var, lower, upper, K_i, K_pop) {
  # converting to sd for computations
  cur_sigma = sqrt(var)

  # are calculations required for mixture probabilities?
  if ( (!is.na(K_pop) | !is.na(K_i)) & (upper != Inf | lower == upper) ) {
    # enter here if K_i and K_pop are provided and individual is NOT a case

    # T and cdf values needed for mixture prob:
    thr_pop = qnorm(K_pop, lower.tail = FALSE)
    cdf_pop = pnorm((thr_pop - mu) / cur_sigma)

    # mixture prob - eq s3 supp notes of PA-FGRS (slightly modified)
    # the supp notes indicate that 1 - mixture prob should have been given as below
    # However, that seemed to lead to incorrect estimates, hence the change.
    mixture_prob = cdf_pop / (cdf_pop + (1 - cdf_pop) * (K_pop - K_i) / K_pop)

  } else {
    # mean and var is updated only through m0 and sd0
    mixture_prob = rep(1, length(K_i))
  }

  # # calculating mixture probabilities
  # w_below = case_when(
  #   K_i == 0 ~ pnorm(upper, mean = mu, sd = cur_sigma),
  #   upper == Inf | upper == lower | is.na(K_i) ~ 1,
  #   TRUE ~ pnorm(upper, mean = mu, sd = cur_sigma) / (1 - pnorm(upper, mean = mu, sd = cur_sigma, lower.tail = FALSE) * K_i / pnorm(upper, lower.tail = FALSE))
  #   )
  # w_above = 1 - w_below

  # new mean - eq s4 supp notes of PA-FGRS
  m0 = tnorm_mean(mu = mu, sigma = cur_sigma, lower = lower, upper = upper)
  m1 = ifelse(upper == Inf | lower == upper,
              0, # if case
              tnorm_mean(mu = mu, sigma = cur_sigma, lower = upper, upper = Inf)) # if control
  new_mean = mixture_prob * m0 + (1 - mixture_prob) * m1

  # new variance - eq S5 supp notes of PA-FGRS
  sd0 = tnorm_var(mu = mu, sigma = cur_sigma, lower = lower, upper = upper)
  sd1 = ifelse(upper == Inf | lower == upper,
               0, # if case
               tnorm_var(mu = mu, sigma = cur_sigma, lower = upper, upper = Inf)) # if control
  new_var = mixture_prob * (m0^2 + sd0) + (1 - mixture_prob) * (m1^2 + sd1) - new_mean^2

  return(list(mean = new_mean, var = new_var))
}


