utils::globalVariables(c("S", "varZ", "mZ", "sum_r", "cur_fid", "k_temp"))

#' Title Helper function for Kendler's FGRS
#'
#' @param tbl tibble with columns K_i, lower, upper, and pid (the personal identifier column).
#' @param cov Kinship matrix with proband as first row and column
#' @param pid column name of personal identifier
#' @param cur_dad_id ID of father (not column name, but the actual ID)
#' @param cur_mom_id ID of mother (not column name, but the actual ID)
#' @param env_cor_sib Cohabitation effect, i.e. Factor by which the siblings are weighted. Defaults to 1.
#' @param env_cor_f  Cohabitation effect, i.e. Factor by which the father is weighted. Defaults to 1.
#' @param env_cor_m  Cohabitation effect, i.e. Factor by which the mother is weighted. Defaults to 1.
#'
#' @returns A tibble with family specific values required for Kendler's FGRS calculation.
#'
#' @importFrom stats var pnorm
#' @importFrom dplyr pull
#'
#' @export
#'
#' @examples
#' # See Vignettes.
kendler_family_calculations = function(tbl,
                                       cov,
                                       pid,
                                       cur_dad_id,
                                       cur_mom_id,
                                       env_cor_sib = 1,
                                       env_cor_f = 1,
                                       env_cor_m = 1) {
  #### this function is kept mostly as-is from Morten's code (PAFGRS package) ####
  # assumes proband's full liab is first row in cov and tbl
  k_temp = cov[-1,1]
  dad_ind = names(k_temp) == cur_dad_id
  mom_ind = names(k_temp) == cur_mom_id
  if (is.na(cur_dad_id)) dad_ind = rep(0,length(k_temp))
  if (is.na(cur_mom_id)) mom_ind = rep(0,length(k_temp))
  # we only use kinship matrix, i.e. h2 = 1
  sib_ind = k_temp == 0.5 & !dad_ind & !mom_ind

  sum_r = sum(k_temp)
  k_temp = k_temp -
    (k_temp * sib_ind * rep(env_cor_sib,length(k_temp))) -
    (k_temp * dad_ind * rep(env_cor_f,length(k_temp))) -
    (k_temp * mom_ind * rep(env_cor_m,length(k_temp)))
  w = pull(tbl, K_i)[-1]/pnorm(-pull(tbl, upper)[-1])# assuming proband is always first
  w[!is.finite(w)] = 1
  lower = pull(tbl, lower)[-1] # assuming proband is always first
  upper = pull(tbl, upper)[-1] # assuming proband is always first
  m_above = function(t) sapply(t, function(t) f = tnorm_mean(mu = 0, sigma = 1, lower =    t, upper = Inf))
  m_below = function(t) sapply(t, function(t) f = tnorm_mean(mu = 0, sigma = 1, lower = -Inf, upper = t))
  z = ifelse(is.finite(lower),
             m_above(lower),
             m_below(upper))

  S <- sum(k_temp * z * w)
  S <- S/sum(k_temp > 0)
  vZ = var(z, na.rm = TRUE)
  tibble(
    !!as.symbol(pid) := tbl[[pid]][[1]],
    S = S,
    varZ = vZ,
    mZ = mean(z),
    sum_r = sum_r)
}


#' (Simplified) Kendler's FGRS
#'
#' Function to calculate the simplified version of Kendler's FGRS based on family data.
#'
#' @param .tbl A matrix, list or data frame that can be converted into a tibble.
#' Must have at least five columns that hold the family identifier, the personal
#' identifier, the role and the lower and upper thresholds. Note that the
#' role must be one of the following abbreviations
#' \itemize{
#'  \item \code{g} (Genetic component of full liability)
#'  \item \code{o} (Full liability)
#'  \item \code{m} (Mother)
#'  \item \code{f} (Father)
#'  \item \code{c[0-9]*.[0-9]*} (Children)
#'  \item \code{mgm} (Maternal grandmother)
#'  \item \code{mgf} (Maternal grandfather)
#'  \item \code{pgm} (Paternal grandmother)
#'  \item \code{pgf} (Paternal grandfather)
#'  \item \code{s[0-9]*} (Full siblings)
#'  \item \code{mhs[0-9]*} (Half-siblings - maternal side)
#'  \item \code{phs[0-9]*} (Half-siblings - paternal side)
#'  \item \code{mau[0-9]*} (Aunts/Uncles - maternal side)
#'  \item \code{pau[0-9]*} (Aunts/Uncles - paternal side).
#' }
#' Defaults to \code{NULL}. If \code{.tbl} is provided, \code{family_graphs} must be \code{NULL}.
#' @param role A string holding the name of the column in \code{.tbl} that
#' holds the role. Each role must be chosen from the following list of abbreviations
#' \itemize{
#'  \item \code{g} (Genetic component of full liability)
#'  \item \code{o} (Full liability)
#'  \item \code{m} (Mother)
#'  \item \code{f} (Father)
#'  \item \code{c[0-9]*.[0-9]*} (Children)
#'  \item \code{mgm} (Maternal grandmother)
#'  \item \code{mgf} (Maternal grandfather)
#'  \item \code{pgm} (Paternal grandmother)
#'  \item \code{pgf} (Paternal grandfather)
#'  \item \code{s[0-9]*} (Full siblings)
#'  \item \code{mhs[0-9]*} (Half-siblings - maternal side)
#'  \item \code{phs[0-9]*} (Half-siblings - paternal side)
#'  \item \code{mau[0-9]*} (Aunts/Uncles - maternal side)
#'  \item \code{pau[0-9]*} (Aunts/Uncles - paternal side).
#' }
#' Defaults to "role".
#' @param family_graphs A tibble with columns pid and family_graph_col, dadcol, and momcol.
#' See prepare_graph for construction of the graphs. The family graphs Defaults to NULL.
#' @param  pid A string holding the name of the column in \code{.tbl} (or \code{family} and
#' \code{threshs}) that hold the personal identifier(s). Defaults to "PID".
#' @param fid A string holding the name of the column in \code{.tbl} or \code{family} that
#' holds the family identifier. Defaults to "fid".
#' @param family_graphs_col Name of column with family graphs in family_graphs. Defaults to "fam_graph".
#' @param dadcol column name of father in family_graphs or .tbl.
#' @param momcol column name of mother in family_graphs or .tbl.
#' @param env_cor_sib Cohabitation effect, i.e. Factor by which the siblings are weighted. Defaults to 1.
#' @param env_cor_f  Cohabitation effect, i.e. Factor by which the father is weighted. Defaults to 1.
#' @param env_cor_m  Cohabitation effect, i.e. Factor by which the mother is weighted. Defaults to 1.
#'
#' @details
#' The coding of the cohabitation effects differ slightly from the one suggested by Kendler et al.
#' Here, it is coded as env_cor_\* = env_eff_\* / (gen_eff_\* + env_eff_\*), while the original implementation
#' suggestes coding it as gen_eff_\*  / (gen_eff_\* + env_eff_\*), i.e. the two are related as 1 - env_cor_\*.
#'
#' @returns A tibble with summary values used to calculate the simplified kendler FGRS and the FGRS itself.
#'
#' @export
#'
#' @examples
#' # See Vignettes.
kendler_simplified = function(.tbl = NULL,
                              family_graphs = NULL,
                              family_graphs_col = "fam_graph",
                              pid = "pid",
                              fid = "fid",
                              role = NULL,
                              dadcol,
                              momcol,
                              env_cor_sib = 0,
                              env_cor_f = 0,
                              env_cor_m = 0) {
  # only one type of input must be used:
  if ( !is.null(.tbl) & !is.null(family_graphs) ) {
    stop("Only one type of input must be used: either .tbl or family_graphs.")
  }

  # Validating input specific variables -------------------------------------

  if ( !is.null(.tbl) ) { #### .tbl input ####
    # Turning .tbl into a tibble
    # if it is not of class tbl
    if (!is.null(.tbl) && !tibble::is_tibble(.tbl))  .tbl <- tibble::as_tibble(.tbl)

    # role (as string) must be supplied
    if (is.null(role)) stop("role must be specified.")
    # if role is supplied, convert to string
    if (!is.null(role)) role <- as.character(role)

    # Checking that .tbl has three columns named pid, fid and role
    if (!(pid %in% colnames(.tbl))) stop(paste0("The column ", pid," does not exist in the tibble .tbl..."))
    if (!(fid %in% colnames(.tbl))) stop(paste0("The column ", fid," does not exist in the tibble .tbl..."))
    if (!(role %in% colnames(.tbl))) stop(paste0("The column ", role," does not exist in the tibble .tbl..."))
    if (!("K_i" %in% colnames(.tbl))) stop("The column K_i does not exist in the tibble .tbl...")

    # In addition, we check that two columns named lower and upper are present
    if (any(!c("lower","upper") %in% colnames(.tbl))) stop("The tibble .tbl must include two columns named 'lower' and 'upper'!")

    # Extracting the (unique) family identifiers
    fam_list <- unique(pull(.tbl, !!as.symbol(fid)))

  } else if ( !is.null(family_graphs) ) { #### Graph input ####

    # validating graph input (nothing is returned)
    return_catch <- validating_graph_input(family_graphs = family_graphs, fid = fid, family_graphs_col = family_graphs_col, useMixture = FALSE)

    # Extracting the (unique) family identifiers
    fam_list <- unique(pull(family_graphs, !!as.symbol(fid)))


  } else ( stop("no valid input used.") )


  # Performing Kendler's calculations per family ----------------------------

  res = lapply(1:length(fam_list), function(i){
    cur_fid = fam_list[i]
    if ( !is.null(family_graphs) ) { # family_graph based covariance construction
      info = extract_estimation_info_graph(
        # no need to pass family_Graphs, family_graphs_Col, or 'i' this way:
        cur_fam_graph = family_graphs[[family_graphs_col]][[i]],
        cur_fid = cur_fid,
        h2 = 1,
        pid = pid,
        useMixture = FALSE,
        add_ind = FALSE)

      # getting parental ids:
      cur_dad_id = family_graphs[[dadcol]][i]
      cur_mom_id = family_graphs[[momcol]][i]

    } else { # role based covariance construction
      # saving to place holder object
      info = extract_estimation_info_tbl(
        .tbl = .tbl,
        pid = pid,
        fid = fid,
        cur_fid = cur_fid,
        role = role,
        useMixture = FALSE,
        h2 = 1,
        add_ind = TRUE)
      # removing genetic liability from both cov mat and tbl
      info$tbl = info$tbl %>% filter(!!as.symbol(role) != "g")
      info$cov = info$cov[rownames(info$cov) != "g", colnames(info$cov) != "g"]

      # getting parental ids:
      cur_dad_id = info$tbl %>% filter(!!as.symbol(role) == "f") %>% pull(!!as.symbol(pid))
      cur_mom_id = info$tbl %>% filter(!!as.symbol(role) == "m") %>% pull(!!as.symbol(pid))
    }

    # prepare and pass input
    kendler_family_calculations(tbl = info$tbl,
                                cov = info$cov,
                                cur_dad_id = cur_dad_id,
                                cur_mom_id = cur_mom_id,
                                env_cor_sib = env_cor_sib,
                                env_cor_f = env_cor_f,
                                env_cor_m = env_cor_m,
                                pid = pid)
  }) %>% do.call("bind_rows", .)

  # finalising calculations that require full population
  if (nrow(res) == 1) {
    vs = 1
    vmz = 1
  } else {
    vs = var(res$S)
    vmz = var(res$mZ)
  }
  total_varZ = mean(res$varZ) + vmz

  # FGRS calculation and return
  res %>% mutate(FGRS = S * vs/(vs + total_varZ/sum_r))
}
