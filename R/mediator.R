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
#' @param m A numeric value indicating the level of the mediator. Default = 0.
#' @param boot_rep A numeric value indicating the number of repetitions
#'   to use when utalizing bootstrap to calculate confidence intervals.
#'   When `boot_rep` = 0, the Delta method for calculating confidence
#'   intervals is used. Default = 0.
#' @param pm_ci A logical indicator for calculating the CI for the proportion
#'   mediated. Default = FALSE. Currently, the CI can only be determined using
#'   boostrapping. If `pm_ci` = TRUE and `boot_rep` = 0 then 100 replicated
#'   are automatically used.
#'
#' @return Tibble containing point estimates and 95 percent CI for the
#'   CDE, NDE, NIE and TE and the point estimate for the proportion mediated.
#'
#' @importFrom rlang .data
#'
#' @export
mediator <- function(data, out.model, med.model, treat, a = 1, a_star = 0,
                     m = 0, boot_rep = 0, pm_ci = FALSE){

  # identifying mediator variable
  mediator <- stringr::str_trim(gsub("~.*","",as.character(med.model$call)[2]))

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
  out_vars <- if (out.reg=="coxph") {
    unlist(stringr::str_extract_all(names(attr(out.model$terms,"dataClasses")),
                                    stringr::boundary("word")))[-1]
  } else {
      names(attr(out.model$terms,"dataClasses"))
    }

  var_set <- unique(c(out_vars,
                      names(attr(med.model$terms,"dataClasses"))))
  data <- data %>% dplyr::select(var_set)

  betas <- stats::coef(med.model) # coefficients from mediation model
  # beta_info <- cov_pred(cmeans, cmodes, treat, mediator, med.model, data)
  beta_info <- cov_pred(treat, mediator, med.model, data)
  betasum <- sum(beta_info$betasum, na.rm=TRUE)
  betameans <- beta_info$betamean

  # get covariate names
  cnames <- names(betameans)

  # Covariance matrix for standar errors
  sigmaV <- stats::sigma(med.model)^2
  Sigma <- comb_sigma(med.model, out.model, treat, mediator,
                                  out.reg, cnames, med.reg)

  # setting coefficients for no interaction = 0 -------------------------------
  if(is.na(out.model$coefficients[paste0(treat, ":", mediator)])){
    out.model$coefficients[paste0(treat, ":", mediator)] <- 0
  } else {
    out.model$coefficients[paste0(treat, ":", mediator)] <-
      out.model$coefficients[paste0(treat, ":", mediator)]
  }

  # pulling coefficients from models
  theta1 <- out.model$coefficients[treat]
  theta2 <- out.model$coefficients[mediator]
  theta3 <- out.model$coefficients[paste0(treat, ":", mediator)]

  beta0 <- med.model$coefficients["(Intercept)"]
  beta1 <- med.model$coefficients[treat]

  arg_list <- list(theta1 = theta1, theta2 = theta2, theta3 = theta3,
                   beta0 = beta0, beta1 = beta1,
                   betasum = betasum, betameans = betameans,
                   a = a, a_star = a_star, m = m, out.reg = out.reg,
                   med.reg = med.reg,
                   sigmaV = sigmaV)

  boot_rep <- dplyr::case_when(
    boot_rep == 0 & pm_ci == FALSE ~ 0,
    boot_rep == 0 & pm_ci == TRUE ~ 100,
    boot_rep != 0  ~ abs(boot_rep),
    TRUE ~ boot_rep
  )

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

      ## proportion mediated
      CI_PM <- c(NA, NA)

    }

    ##### ------------------------------------------------------------------------------ #####
    # CI using boostrapping
    ##### ------------------------------------------------------------------------------ #####

    if (boot_rep > 0) {

      pb <- progress::progress_bar$new(total = boot_rep  + (boot_rep + 1) * 15)

      boot_dat <-
        rsample::bootstraps(data, times = boot_rep, apparent = TRUE) %>%
        dplyr::mutate(out = furrr::future_map(.data$splits, function(split, pb) {
          x <- stats::update(out.model, data = rsample::analysis(split))
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(med = furrr::future_map(.data$splits, function(split, pb) {
          x <- stats::update(med.model, data = rsample::analysis(split))
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(beta_info = furrr::future_map2(.data$splits, .data$med,
                                                     function(split, med, pb) {
          x <- cov_pred(treat = treat, mediator = mediator, med.model = med,
                        data = rsample::analysis(split))
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(betasum = furrr::future_map(beta_info, function(betainfo, pb) {
          x <- sum(betainfo$betasum, na.rm=TRUE)
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(sigmaV = furrr::future_map(.data$med, function(med, pb) {
          x <- stats::sigma(med)^2
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(theta1 = furrr::future_map(.data$out, function(out, pb) {
          x <- out$coefficients[treat]
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(theta2 = furrr::future_map(.data$out, function(out, pb) {
          x <- out$coefficients[mediator]
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(theta3 = furrr::future_map(.data$out, function(out, pb){
          x <- out$coefficients[paste0(treat, ":", mediator)]
          x <- dplyr::case_when(is.na(x) ~ 0, TRUE ~ as.numeric(x))
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(beta0 = furrr::future_map(.data$med, function(med, pb) {
          x <- med$coefficients["(Intercept)"]
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(beta1 = furrr::future_map(.data$med, function(med, pb) {
          x <- med$coefficients[treat]
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(CDE = furrr::future_map2(theta1, theta3, function(theta1, theta3, pb){
          x <- controlled_direct_effect(theta1 = theta1, a = a, a_star = a_star,
                                        theta3 = theta3, m = m, out.reg = out.reg)
          pb$tick()
          return(x)
        }, pb = pb))%>%
        dplyr::mutate(NDE = furrr::future_pmap(list(theta1, theta2, theta3, beta0,
                               beta1, betasum, sigmaV),
                          function(theta1, theta2, theta3, beta0,
                                   beta1, betasum, sigmaV, pb){
          x <- natural_direct_effect(theta1 = theta1, a = a, theta2 = theta2,
                                       theta3 = theta3, beta0 = beta0,
                                       beta1 = beta1, a_star = a_star,
                                       betasum = betasum, sigmaV = sigmaV,
                                       out.reg = out.reg, med.reg = med.reg)
          pb$tick()
          return(x)
        }, pb)) %>%
        dplyr::mutate(NIE = furrr::future_pmap(list(beta0, beta1, betasum, theta2, theta3),
                          function(beta0, beta1, betasum, theta2, theta3, pb){
          x <- natural_indirect_effect(out.reg = out.reg, med.reg = med.reg,
                                         beta0 = beta0, beta1 = beta1,
                                         a_star = a_star, betasum = betasum,
                                         theta2 = theta2, theta3 = theta3, a = a)
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(TE = furrr::future_map2(NDE, NIE, function(NDE, NIE, pb){
          x <- total_effect(NDE = NDE, NIE = NIE, out.reg = out.reg)
          pb$tick()
          return(x)
        }, pb = pb)) %>%
        dplyr::mutate(PM = furrr::future_pmap(list(NDE, NIE, TE), function(NDE, NIE, TE, pb){
          x <- prop_mediated(NDE = NDE, NIE = NIE, out.reg = out.reg, TE = TE)
          pb$tick()
          return(x)
        }, pb = pb))

      CI_CDE <- boot_dat %>%
        dplyr::pull(CDE) %>%
        unlist() %>%
        stats::quantile(probs = c(0.05 / 2, 1 - 0.05 / 2), na.rm = TRUE)
      CI_NDE <- boot_dat %>%
        dplyr::pull(NDE) %>%
        unlist() %>%
        stats::quantile(probs = c(0.05 / 2, 1 - 0.05 / 2), na.rm = TRUE)
      CI_NIE <- boot_dat %>%
        dplyr::pull(NIE) %>%
        unlist() %>%
        stats::quantile(probs = c(0.05 / 2, 1 - 0.05 / 2), na.rm = TRUE)
      CI_TE  <- boot_dat %>%
        dplyr::pull(TE) %>%
        unlist() %>%
        stats::quantile(probs = c(0.05 / 2, 1 - 0.05 / 2), na.rm = TRUE)
      CI_PM  <- boot_dat %>%
        dplyr::pull(PM) %>%
        unlist() %>%
        stats::quantile(probs = c(0.05 / 2, 1 - 0.05 / 2), na.rm = TRUE)




      # CIs <- function(data , indices, ...) {
      #
      #   d <- data[indices,]
      #   # pb$tick()
      #
      #   out <- stats::update(out.model, data = d)
      #   med <- stats::update(med.model, data = d)
      #
      #   betas <- stats::coef(med) # coefficients from mediation model
      #   beta_info <- cov_pred(treat, mediator, med, d)
      #   betasum <- sum(beta_info$betasum, na.rm=TRUE)
      #   betameans <- beta_info$betamean
      #   cnames <- names(betameans)
      #
      #   # Covariance matrix for standar errors
      #   sigmaV <- stats::sigma(med.model)^2
      #
      #   # setting coefficients for no interaction = 0 -------------------------------
      #   if(is.na(out$coefficients[paste0(treat, ":", mediator)])){
      #     out$coefficients[paste0(treat, ":", mediator)] <- 0
      #   } else {
      #     out$coefficients[paste0(treat, ":", mediator)] <-
      #       out$coefficients[paste0(treat, ":", mediator)]
      #   }
      #
      #   # pulling coefficients from models
      #   theta1 <- out$coefficients[treat]
      #   theta2 <- out$coefficients[mediator]
      #   theta3 <- out$coefficients[paste0(treat, ":", mediator)]
      #
      #   beta0 <- med$coefficients["(Intercept)"]
      #   beta1 <- med$coefficients[treat]
      #
      #   arg_list <- list(theta1 = theta1, theta2 = theta2,
      #                    theta3 = theta3, beta0 = beta0,
      #                    beta1 = beta1, a = a, a_star = a_star,
      #                    m = m, betasum = betasum, sigmaV = sigmaV,
      #                    out.reg = out.reg, med.reg = med.reg)
      #
      #     # calculate effect estimates
      #
      #     ## controlled direct effect
      #     CDE <- do.call(controlled_direct_effect, arg_list)
      #
      #     ## natural direct effect
      #     NDE <- do.call(natural_direct_effect, arg_list)
      #
      #     ## natural indirect effect
      #     NIE <- do.call(natural_indirect_effect, arg_list)
      #
      #     ## total effect
      #     TE <-  total_effect(NDE, NIE, out.reg)
      #
      #     ## total effect
      #     PM <- prop_mediated(NDE, NIE, out.reg, TE)
      #
      #   val <- c(CDE, NDE, NIE, TE, PM)
      #
      #   return(val)
      # }
      #
      # boot_results <- boot::boot(data = data, statistic = CIs, R = boot_rep,
      #                            parallel = "multicore",
      #                            ncpus = parallel::detectCores(logical = FALSE))
      #
      # # boot_results <- boot::boot(data = data, statistic = CIs, R = boot_rep)
      #
      # CI_CDE <- c(boot::boot.ci(boot_results, index = 1, type = "bca")$bca[[4]],
      #             boot::boot.ci(boot_results, index = 1, type = "bca")$bca[[5]])
      # CI_NDE <- c(boot::boot.ci(boot_results, index = 2, type = "bca")$bca[[4]],
      #             boot::boot.ci(boot_results, index = 2, type = "bca")$bca[[5]])
      # CI_NIE <- c(boot::boot.ci(boot_results, index = 3, type = "bca")$bca[[4]],
      #             boot::boot.ci(boot_results, index = 3, type = "bca")$bca[[5]])
      # CI_TE  <- c(boot::boot.ci(boot_results, index = 4, type = "bca")$bca[[4]],
      #             boot::boot.ci(boot_results, index = 4, type = "bca")$bca[[5]])
      # CI_PM  <- c(boot::boot.ci(boot_results, index = 5, type = "bca")$bca[[4]],
      #             boot::boot.ci(boot_results, index = 5, type = "bca")$bca[[5]])

    }


  # }

  output <- tibble::tibble(Effect = c("CDE", "NDE", "NIE",
                               "Total Effect", "Proportion Mediated"),
                   Estimate = c(round(CDE, 5), round(NDE, 5),
                                round(NIE, 5), round(TE, 5),
                                round(PM, 5)),
                   `Lower 95% CI` = c(round(CI_CDE[[1]], 5),
                                      round(CI_NDE[[1]], 5),
                                      round(CI_NIE[[1]], 5),
                                      round(CI_TE[[1]], 5),
                                      round(CI_PM[[1]], 5)),
                   `Upper 95% CI` = c(round(CI_CDE[[2]], 5),
                                      round(CI_NDE[[2]], 5),
                                      round(CI_NIE[[2]], 5),
                                      round(CI_TE[[2]], 5),
                                      round(CI_PM[[2]], 5)))
  # rownames(output) <- NULL

  return(output)
}
