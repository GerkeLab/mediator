#' Causal Mediation Analysis
#'
#' The `mediator` R function conducts mediation analysis under
#' the counterfactual framework assuming interation between the exposure
#' and mediator. Currently the function works for binary and continuous
#' outcomes and mediators.
#'
#' @param data Data set to use for analysis
#' @param out.model A fitted model object for the outcome.
#'   Can be of class 'glm','lm', or 'coxph'.
#' @param med.model A fitted model object for the mediator.
#'   Can be of class 'glm','lm'.
#' @param treat A character string indicating the name of the
#'   treatment/exposure variable used.
#' @param a A numeric value indicating the exposure level. Default = 1
#' @param a_star A numeric value indicating the compared exposure level.
#'   Default = 0.
#' @param m A numeric value indicating the level of the mediator. Default = 1
#' @param boot_rep A numeric value indicating the number of repetitions
#'   to use when utalizing bootstrap to calculate confidence intervals.
#'   When `boot_rep` = 0, the Delta method for calculating confidence
#'   intervals is used. Default = 0.
#'
#' @return Tibble containing point estimates and 95 percent CI for the
#'   CDE, NDE, NIE and TE and the point estimate for the proportion mediated.
#'

mediator <- function(...) {
  UseMethod("mediator")
}

#'
#' @export
mediator.default <- function(data,out.model, med.model, treat, a = 1, a_star = 0,
                     m = 0, boot_rep = 0){

  # identifying mediator variable
  mediator_name <- stringr::str_trim(gsub("~.*","",as.character(med.model$call)[2]))

  # identifying out model type
  out.reg <- if (class(out.model)[1] == "coxph") {
    "coxph"
  } else if (class(out.model)[1] == "lm") {
      "linear"
  } else if (class(out.model)[1] == "glm") {
      if (out.model$family$family == "binomial") "logistic" else
        "linear"
  }

  med.reg <- if (class(med.model)[1] == "lm") {
    "linear"
  } else if (class(med.model)[1] == "glm") {
    if (med.model$family$family == "binomial") "logistic" else
      "linear"
  }


  # subset data to the set of variables from data which are relevant
  out_vars <- if (out.reg=="coxph") names(attr(out.model$terms,"dataClasses"))[-1] else
    names(attr(out.model$terms,"dataClasses"))

  var_set <- unique(c(out_vars,
                      names(attr(med.model$terms,"dataClasses"))))
  data <- data %>% dplyr::select(var_set)

  betas <- stats::coef(med.model) # coefficients from mediation model
  beta_info <- cov_pred(treat, mediator_name, med.model, data)
  betasum <- sum(beta_info$betasum, na.rm=TRUE)
  betameans <- beta_info$betamean

  # get covariate names
  cnames <- names(betameans)

  # Covariance matrix for standar errors
  sigmaV <- stats::sigma(med.model)^2
  Sigma <- comb_sigma(med.model, out.model, treat, mediator_name,
                                  out.reg, cnames, med.reg)

  # setting coefficients for no interaction = 0 -------------------------------
  if(is.na(out.model$coefficients[paste0(treat, ":", mediator_name)])){
    out.model$coefficients[paste0(treat, ":", mediator_name)] <- 0
  } else {
    out.model$coefficients[paste0(treat, ":", mediator_name)] <-
      out.model$coefficients[paste0(treat, ":", mediator_name)]
  }

  # pulling coefficients from models
  theta1 <- out.model$coefficients[treat]
  theta2 <- out.model$coefficients[mediator_name]
  theta3 <- out.model$coefficients[paste0(treat, ":", mediator_name)]

  beta0 <- med.model$coefficients["(Intercept)"]
  beta1 <- med.model$coefficients[treat]

  arg_list <- list(theta1 = theta1, theta2 = theta2, theta3 = theta3,
                   beta0 = beta0, beta1 = beta1,
                   betasum = betasum, betameans = betameans,
                   a = a, a_star = a_star, m = m, out.reg = out.reg,
                   med.reg = med.reg,
                   sigmaV = sigmaV)

  ##### ----------------------------------------------------------------- #####
  # Calculate effect estimates and confidence intervals (delta method) --------
  ##### ----------------------------------------------------------------- #####

    # calculate effect estimates
    ## controlled direct effect
    CDE <- do.call(controlled_direct_effect, arg_list)

    ## natural direct effect
    NDE <- do.call(natural_direct_effect, arg_list)

    ## natural indirect effect
    NIE <- do.call(natural_indirect_effect, arg_list)

    ## total effect
    TE <- total_effect(NDE, NIE, out.reg)

    ## proportion mediated
    PM <- prop_mediated(NDE, NIE, out.reg, TE)

    # calculate GAMMA for SE
    if (boot_rep == 0) {
      ## controlled direct effect
      gCDE <- do.call(gamma_cde, arg_list)

      ## natural direct effect
      gNDE <- do.call(gamma_nde, arg_list)

      ## natural indirect effect
      gNIE <- do.call(gamma_nie, arg_list)

      ## total effect
      gTE <- do.call(gamma_te, c(arg_list, list("gNDE" = gNDE, "gNIE" = gNIE)))

      # delta method of calculating confidence intervals

      ## controlled direct effect
      CI_CDE <- delta_cde_nde(gCDE, Sigma, CDE, a, a_star, out.reg, med.reg)

      ## natural direct effect
      CI_NDE <- delta_cde_nde(gNDE, Sigma, NDE, a, a_star, out.reg, med.reg)

      ## natural indirect effect
      CI_NIE <- delta_nie(gNIE, Sigma, NIE, a, a_star, out.reg, med.reg)

      ## total effect
      CI_TE <- delta_te(gTE, Sigma, TE, a, a_star, out.reg, med.reg)

    }

    ##### ------------------------------------------------------------------------------ #####
    # CI using boostrapping
    ##### ------------------------------------------------------------------------------ #####

    if (boot_rep > 0) {

      CIs <- function(data , indices, ...) {

        d <- data[indices, ]

        out <- stats::update(out.model, data = d)
        med <- stats::update(med.model, data = d)

        betas <- stats::coef(med) # coefficients from mediation model
        beta_info <- cov_pred(treat, mediator_name, med, d)
        betasum <- sum(beta_info$betasum, na.rm=TRUE)
        betameans <- beta_info$betamean
        cnames <- names(betameans)

        # Covariance matrix for standar errors
        sigmaV <- stats::sigma(med.model)^2
        Sigma <- comb_sigma(med.model, out.model, treat, mediator_name,
                            out.reg, cnames, med.reg)

        # setting coefficients for no interaction = 0 -------------------------------
        if(is.na(out$coefficients[paste0(treat, ":", mediator_name)])){
          out$coefficients[paste0(treat, ":", mediator_name)] <- 0
        } else {
          out$coefficients[paste0(treat, ":", mediator_name)] <-
            out$coefficients[paste0(treat, ":", mediator_name)]
        }

        # pulling coefficients from models
        theta1 <- out$coefficients[treat]
        theta2 <- out$coefficients[mediator_name]
        theta3 <- out$coefficients[paste0(treat, ":", mediator_name)]
        thetas <- list(theta1 = unname(theta1), theta2 = unname(theta2),
                       theta3 = unname(theta3))

        beta0 <- med$coefficients["(Intercept)"]
        beta1 <- med$coefficients[treat]

        arg_list <- list(theta1, theta2, theta3, beta0, beta1, a, a_star, m,
                         betasum, sigmaV, out.reg, med.reg)

          # calculate effect estimates

          ## controlled direct effect
          CDE <- do.call(controlled_direct_effect, arg_list)

          ## natural direct effect
          NDE <- do.call(natural_direct_effect, arg_list)

          ## natural indirect effect
          NIE <- do.call(natural_indirect_effect, arg_list)

          ## total effect
          TE <-  total_effect(NDE, NIE, out.reg)

        val <- c(CDE, NDE, NIE, TE)

        return(val)
      }

      boot_results <- boot::boot(data = data, statistic = CIs, R = boot_rep,
                                 parallel = "multicore",
                                 ncpus = parallel::detectCores(logical = FALSE))

      CI_CDE <- c(boot::boot.ci(boot_results, index = 1, type = "bca")$bca[[4]],
                  boot::boot.ci(boot_results, index = 1, type = "bca")$bca[[5]])
      CI_NDE <- c(boot::boot.ci(boot_results, index = 2, type = "bca")$bca[[4]],
                  boot::boot.ci(boot_results, index = 2, type = "bca")$bca[[5]])
      CI_NIE <- c(boot::boot.ci(boot_results, index = 3, type = "bca")$bca[[4]],
                  boot::boot.ci(boot_results, index = 3, type = "bca")$bca[[5]])
      CI_TE  <- c(boot::boot.ci(boot_results, index = 4, type = "bca")$bca[[4]],
                  boot::boot.ci(boot_results, index = 4, type = "bca")$bca[[5]])

    }

  # }

  estimates <- tibble::tibble(Effect = c("CDE", "NDE", "NIE",
                               "Total Effect", "Proportion Mediated"),
                   Estimate = c(CDE, NDE, NIE, TE, PM),
                   `Lower 95% CI` = c(CI_CDE[[1]],
                                      CI_NDE[[1]],
                                      CI_NIE[[1]],
                                      CI_TE[[1]],NA),
                   `Upper 95% CI` = c(CI_CDE[[2]],
                                      CI_NDE[[2]],
                                      CI_NIE[[2]],
                                      CI_TE[[2]], NA))

  # create mediator results object
  mediator <- estimates
  class(mediator) <- c("mediator", class(mediator))

  attr(mediator,"arguments") <- list(treat = treat, mediator = mediator_name,
                                    a = a, a_star = a_star, m = m)
  attr(mediator,"gammas") <- c(list(CDE = gCDE), list(NDE = gNDE),
                               list(NIE = gNIE), list(TE = gTE))
  attr(mediator,"sigma") <- Sigma

  return(mediator)
}
