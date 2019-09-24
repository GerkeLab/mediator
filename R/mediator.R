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
#' @param mediator A character string indicating the name of the
#'   mediator variable used.
#' @param out.reg A character string indicating the type of
#'   regression used in the outcome model. Can be either "logisitic",
#'   "linear", or "coxph".
#' @param med.reg A character string indicating the type of
#'   regression used in the mediator model. Can be either "logisitic" or
#'   "linear".
#' @param a A numeric value indicating the exposure level. Default = 1
#' @param a_star A numeric value indicating the compared exposure level.
#'   Default = 0.
#' @param m A numeric value indicating the level of the mediator. Default = 1
#' @param boot_rep A numeric value indicating the number of repetitions
#'   to use when utalizing bootstrap to calculate confidence intervals.
#'   When `boot_rep` = 0, the Delta method for calculating confidence
#'   intervals is used. Default = 0.
#'
#' @return Data frame containing point estimates and 95 percent CI for the
#'   CDE, NDE, NIE and TE and the point estimate for the proportion mediated.
#'
#' @export
mediator <- function(data,
                     out.model,
                     med.model,
                     treat,
                     mediator,
                     out.reg = "logistic",
                     med.reg = "logistic",
                     a = 1,
                     a_star = 0,
                     m = 0,
                     boot_rep = 0){

  # subset data to the set of variables from data which are relevant
  var_set <- unique(c(names(attr(out.model$terms,"dataClasses")),
                      names(attr(med.model$terms,"dataClasses"))))
  data <- data %>% dplyr::select(var_set)

  # calculate mean of numeric values, mode of characters/factors
  # cmeans <- data %>%
  #   dplyr::select_if(is.numeric) %>%
  #   purrr::map_dbl(~mean(.x, na.rm = TRUE)) # mean value for all numeric values
  # cmodes <- data %>%
  #   dplyr::select_if(purrr::negate(is.numeric)) %>%
  #   purrr::map_chr(~{
  #     ux <- unique(.x)
  #     ux[which.max(tabulate(match(.x, ux)))]
  #     })

  betas <- stats::coef(med.model) # coefficients from mediation model
  # beta_info <- cov_pred(cmeans, cmodes, treat, mediator, med.model, data)
  beta_info <- cov_pred(treat, mediator, med.model, data)
  betasum <- sum(beta_info$betasum, na.rm=TRUE)
  betameans <- beta_info$betamean

  # get covariate names
  cnames <- names(betameans)

  # Covariance matrix for standar errors
  SigmaB <- stats::vcov(med.model)
  SigmaT <- stats::vcov(out.model)
  # including 0 for variance for interaction if missing
  if(is.na(out.model$coefficients[paste0(treat, ":", mediator)])){
    SigmaT <- rbind(cbind(SigmaT,rep(0,nrow(SigmaT))),rep(0,nrow(SigmaT)))
    dimnames(SigmaT)[[1]][nrow(SigmaT)] <- paste0(treat, ":", mediator)
    dimnames(SigmaT)[[2]][nrow(SigmaT)] <- paste0(treat, ":", mediator)
  } else if (out.reg == "coxph") {
    SigmaT <- rbind(cbind(rep(0,nrow(SigmaT)),SigmaT),rep(0,nrow(SigmaT)))
    dimnames(SigmaT)[[1]][nrow(SigmaT)] <- "(Intercept)"
    dimnames(SigmaT)[[2]][1] <- "(Intercept)"
  } else {
    SigmaT <- SigmaT
  }
  SigmaT <- SigmaT[c("(Intercept)", treat, mediator, paste0(treat, ":", mediator), cnames),
                   c("(Intercept)", treat, mediator, paste0(treat, ":", mediator), cnames)]
  Sigma <- rbind(cbind(SigmaB, matrix(0, ncol = ncol(SigmaT), nrow = nrow(SigmaB))),
                 cbind(matrix(0, ncol = ncol(SigmaB), nrow = nrow(SigmaT)), SigmaT))
  # Sigma includes standard error only for logistic/linear and no others
  if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {
     sigmaV <- stats::sigma(med.model)^2
     Sigma <- rbind(cbind(Sigma, rep(0, nrow(Sigma))),
                    c(rep(0, ncol(Sigma)), sigmaV))
  } else {
     Sigma <- Sigma
  }

  rm(SigmaB, SigmaT)

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

  ##### ----------------------------------------------------------------- #####
  # Calculate effect estimates and confidence intervals (delta method) --------
  ##### ----------------------------------------------------------------- #####
  if (out.reg %in% c("logistic","coxph") & med.reg == "logistic") {

    # calculate effect estimates
    ## controlled direct effect
    CDE <- controlled_direct_effect(theta1, a, a_star, theta3, m, out.reg)

    ## natural direct effect
    NDE <- natural_direct_effect(theta1, a, theta2, theta3, beta0,
                                 beta1, a_star, betasum, sigmaV,
                                 out.reg, med.reg)

    ## natural indirect effect
    NIE <- natural_indirect_effect(out.reg, med.reg, beta0, beta1, a_star,
                                   betasum, theta2, theta3, a)

    ## total effect
    TE <- total_effect(NDE, NIE, out.reg)

    ## proportion mediated
    PM <- prop_mediated(NDE, NIE, out.reg, TE)

    # calculate GAMMA for SE
    if (boot_rep == 0) {
      ## controlled direct effect
      gCDE <- gamma_cde(betameans, a, a_star, m, out.reg, med.reg)

      ## natural direct effect
      gNDE <- gamma_nde(out.reg, med.reg, theta2, theta3, a, beta0,
                         beta1, a_star, betasum, betameans, sigmaV)

      ## natural indirect effect
      gNIE <- gamma_nie(beta0, beta1, a, betasum, a_star, out.reg,
                        med.reg, theta2, theta3, betameans)

      ## total effect
      gTE <- gamma_te(gNDE, gNIE, beta0, beta1, a, betasum, a_star,
                      theta3, theta2, betameans, out.reg, med.reg, sigmaV)

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

  } else if (out.reg == "linear" & med.reg == "logistic") {

    # calculate effect estimates

    ## controlled direct effect
    CDE <- controlled_direct_effect(theta1, a, a_star, theta3, m, out.reg)

    ## natural direct effect
    NDE <- natural_direct_effect(theta1, a, theta2, theta3, beta0,
                                 beta1, a_star, betasum, sigmaV,
                                 out.reg, med.reg)

    ## natural indirect effect
    NIE <- natural_indirect_effect(out.reg, med.reg, beta0, beta1, a_star,
                                   betasum, theta2, theta3, a)

    ## total effect
    TE <- total_effect(NDE, NIE, out.reg)

    ## proportion mediated
    PM <- prop_mediated(NDE, NIE, out.reg, TE)

    # calculate GAMMA for SE
    if (boot_rep == 0) {

      ## controlled direct effect
      gCDE <- gamma_cde(betameans, a, a_star, m, out.reg, med.reg)

      ## natural direct effect
      gNDE <- gamma_nde(out.reg, med.reg, theta2, theta3, a, beta0,
                        beta1, a_star, betasum, betameans, sigmaV)

      ## natural indirect effect
      gNIE <- gamma_nie(beta0, beta1, a, betasum, a_star, out.reg,
                        med.reg, theta2, theta3, betameans)

      ## total effect - exchange A for Q 07/01/2019
      gTE <- gamma_te(gNDE, gNIE, beta0, beta1, a, betasum, a_star,
                      theta3, theta2, betameans, out.reg, med.reg, sigmaV)

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

  } else if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {

    # calculate effect estimates

    ## controlled direct effect
    CDE <- controlled_direct_effect(theta1, a, a_star, theta3, m, out.reg)

    ## natural direct effect
    NDE <- natural_direct_effect(theta1, a, theta2, theta3, beta0,
                                 beta1, a_star, betasum, sigmaV,
                                 out.reg, med.reg)

    ## natural indirect effect
    NIE <- natural_indirect_effect(out.reg, med.reg, beta0, beta1, a_star,
                                   betasum, theta2, theta3, a)

    ## total effect
    TE <-  total_effect(NDE, NIE, out.reg)

    ## proportion mediated
    PM <- prop_mediated(NDE, NIE, out.reg, TE)

    # calculate GAMMA for SE
    if (boot_rep == 0) {
      ## controlled direct effect
      gCDE <- gamma_cde(betameans, a, a_star, m, out.reg, med.reg)

      ## natural direct effect
      gNDE <- gamma_nde(out.reg, med.reg, theta2, theta3, a, beta0,
                        beta1, a_star, betasum, betameans, sigmaV)

      ## natural indirect effect
      gNIE <- gamma_nie(beta0, beta1, a, betasum, a_star, out.reg,
                        med.reg, theta2, theta3, betameans)

      ## total effect
      gTE <- gamma_te(gNDE, gNIE, beta0, beta1, a, betasum, a_star,
                      theta3, theta2, betameans, out.reg, med.reg, sigmaV)

      # calculate CI using delta method

      ## controlled direct effect
      CI_CDE <- delta_cde_nde(gCDE, Sigma, CDE, a, a_star, out.reg, med.reg)

      ## natural direct effect
      CI_NDE <- delta_cde_nde(gNDE, Sigma, NDE, a, a_star, out.reg, med.reg)

      ## natural indirect effect
      CI_NIE <- delta_nie(gNIE, Sigma, NIE, a, a_star, out.reg, med.reg)

      ## total effect
      CI_TE <- delta_te(gTE, Sigma, TE, a, a_star, out.reg, med.reg)

    }

  } else {

    # calculate effect estimates

    ## controlled direct effect
    CDE <- controlled_direct_effect(theta1, a, a_star, theta3, m, out.reg)

    ## natural direct effect
    NDE <- natural_direct_effect(theta1, a, theta2, theta3, beta0,
                                 beta1, a_star, betasum, sigmaV,
                                 out.reg, med.reg)

    ## natural indirect effect
    NIE <- natural_indirect_effect(out.reg, med.reg, beta0, beta1, a_star,
                                   betasum, theta2, theta3, a)

    ## total effect
    TE <-  total_effect(NDE, NIE, out.reg)

    ## proportion mediated
    PM <- prop_mediated(NDE, NIE, out.reg, TE)

    # calculate GAMMA for SE
    if (boot_rep == 0) {

      ## controlled direct effect
      gCDE <- gamma_cde(betameans, a, a_star, m, out.reg, med.reg)

      ## natural direct effect
      gNDE <- gamma_nde(out.reg, med.reg, theta2, theta3, a, beta0,
                        beta1, a_star, betasum, betameans, sigmaV)

      ## pure natural indirect effect
      gNIE <- gamma_nie(beta0, beta1, a, betasum, a_star, out.reg,
                        med.reg, theta2, theta3, betameans)

      ## total effects
      gTE <- gamma_te(gNDE, gNIE, beta0, beta1, a, betasum, a_star,
                      theta3, theta2, betameans, out.reg, med.reg, sigmaV)

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

  }
    ##### ------------------------------------------------------------------------------ #####
    # CI using boostrapping
    ##### ------------------------------------------------------------------------------ #####

    if (boot_rep > 0) {

      CIs <- function(data , indices, ...) {

        d <- data[indices, ]

        out <- stats::update(out.model, data = d)
        med <- stats::update(med.model, data = d)

        # cmeans <- d %>%
        #   dplyr::select_if(is.numeric) %>%
        #   purrr::map_dbl(~mean(.x, na.rm = TRUE)) # mean value for all numeric values
        # cmodes <- d %>%
        #   dplyr::select_if(purrr::negate(is.numeric)) %>%
        #   purrr::map_chr(~{
        #     ux <- unique(.x)
        #     ux[which.max(tabulate(match(.x, ux)))]
        #   })

        betas <- stats::coef(med) # coefficients from mediation model
        # beta_info <- cov_pred(cmeans, cmodes, treat, mediator, med, d)
        beta_info <- cov_pred(treat, mediator, med, d)
        betasum <- sum(beta_info$betasum, na.rm=TRUE)
        betameans <- beta_info$betamean
        cnames <- names(betameans)

        # Covariance matrix for standar errors
        SigmaB <- stats::vcov(med)
        SigmaT <- stats::vcov(out)
        # including 0 for variance for interaction if missing
        if(is.na(out$coefficients[paste0(treat, ":", mediator)])){
          SigmaT <- rbind(cbind(SigmaT,rep(0,nrow(SigmaT))),rep(0,nrow(SigmaT)))
          dimnames(SigmaT)[[1]][nrow(SigmaT)] <- paste0(treat, ":", mediator)
          dimnames(SigmaT)[[2]][nrow(SigmaT)] <- paste0(treat, ":", mediator)
        } else if (out.reg == "coxph") {
          SigmaT <- rbind(cbind(rep(0,nrow(SigmaT)),SigmaT),rep(0,nrow(SigmaT)))
          dimnames(SigmaT)[[1]][nrow(SigmaT)] <- "(Intercept)"
          dimnames(SigmaT)[[2]][1] <- "(Intercept)"
        } else {
          SigmaT <- SigmaT
        }
        SigmaT <- SigmaT[c("(Intercept)", treat, mediator, paste0(treat, ":", mediator), cnames),
                         c("(Intercept)", treat, mediator, paste0(treat, ":", mediator), cnames)]
        Sigma <- rbind(cbind(SigmaB, matrix(0, ncol = ncol(SigmaT), nrow = nrow(SigmaB))),
                       cbind(matrix(0, ncol = ncol(SigmaB), nrow = nrow(SigmaT)), SigmaT))
        # Sigma includes standard error only for logistic/linear and no others
        if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {
          sigmaV <- stats::sigma(med)^2
          Sigma <- rbind(cbind(Sigma, rep(0, nrow(Sigma))),
                         c(rep(0, ncol(Sigma)), sigmaV))
        } else {
          Sigma <- Sigma
        }

        rm(SigmaB, SigmaT)

        # setting coefficients for no interaction = 0 -------------------------------
        if(is.na(out$coefficients[paste0(treat, ":", mediator)])){
          out$coefficients[paste0(treat, ":", mediator)] <- 0
        } else {
          out$coefficients[paste0(treat, ":", mediator)] <-
            out$coefficients[paste0(treat, ":", mediator)]
        }

        # pulling coefficients from models
        theta1 <- out$coefficients[treat]
        theta2 <- out$coefficients[mediator]
        theta3 <- out$coefficients[paste0(treat, ":", mediator)]

        beta0 <- med$coefficients["(Intercept)"]
        beta1 <- med$coefficients[treat]

        if (out.reg == "logistic" & med.reg == "logistic") {

          # calculate effect estimates

          ## controlled direct effect
          CDE <- controlled_direct_effect(theta1, a, a_star, theta3, m, out.reg)

          ## natural direct effect
          NDE <- natural_direct_effect(theta1, a, theta2, theta3, beta0,
                                       beta1, a_star, betasum, sigmaV,
                                       out.reg, med.reg)

          ## natural indirect effect
          NIE <- natural_indirect_effect(out.reg, med.reg, beta0, beta1, a_star,
                                         betasum, theta2, theta3, a)

          ## total effect
          TE <-  total_effect(NDE, NIE, out.reg)

        } else if (out.reg == "linear" & med.reg == "logistic") {

          # calculate effect estimates

          ## controlled direct effect
          CDE <-  controlled_direct_effect(theta1, a, a_star, theta3, m, out.reg)

          ## natural direct effect
          NDE <- natural_direct_effect(theta1, a, theta2, theta3, beta0,
                                       beta1, a_star, betasum, sigmaV,
                                       out.reg, med.reg)

          ## natural indirect effect
          NIE <- natural_indirect_effect(out.reg, med.reg, beta0, beta1, a_star,
                                         betasum, theta2, theta3, a)

          ## total effect
          TE <-  total_effect(NDE, NIE, out.reg)

        } else if (out.reg == "logistic" & med.reg == "linear") {

          # calculate effect estimates

          ## controlled direct effect
          CDE <-  controlled_direct_effect(theta1, a, a_star, theta3, m, out.reg)

          ## natural direct effect
          NDE <- natural_direct_effect(theta1, a, theta2, theta3, beta0,
                                       beta1, a_star, betasum, sigmaV,
                                       out.reg, med.reg)

          ## natural indirect effect
          NIE <- natural_indirect_effect(out.reg, med.reg, beta0, beta1, a_star,
                                         betasum, theta2, theta3, a)

          ## total effect
          TE <-  total_effect(NDE, NIE, out.reg)

        } else {

          # calculate effect estimates

          ## controlled direct effect
          CDE <- controlled_direct_effect(theta1, a, a_star, theta3, m, out.reg)

          ## natural direct effect
          NDE <- natural_direct_effect(theta1, a, theta2, theta3, beta0,
                                       beta1, a_star, betasum, sigmaV,
                                       out.reg, med.reg)

          ## natural indirect effect
          NIE <- natural_indirect_effect(out.reg, med.reg, beta0, beta1, a_star,
                                         betasum, theta2, theta3, a)

          ## total effect
          TE <-  total_effect(NDE, NIE, out.reg)

        }

        val <- c(CDE, NDE, NIE, TE)

        return(val)
      }


      # library(boot)
      boot_results <- boot::boot(data = data, statistic = CIs, R = boot_rep,
                                 parallel = "multicore",
                                 ncpus = parallel::detectCores(logical = FALSE))

      CI_CDE <- c(boot::boot.ci(boot_results, index = 1, type = "bca")$bca[[4]],
                  boot::boot.ci(boot_results, index = 1, type = "bca")$bca[[5]])
      CI_NDE <- c(boot::boot.ci(boot_results, index = 2, type = "bca")$bca[[4]],
                  boot::boot.ci(boot_results, index = 2, type = "bca")$bca[[5]])
      CI_NIE <- c(boot::boot.ci(boot_results, index = 3, type = "bca")$bca[[4]],
                  boot::boot.ci(boot_results, index = 3, type = "bca")$bca[[5]])
      CI_TE <- c(boot::boot.ci(boot_results, index = 4, type = "bca")$bca[[4]],
                 boot::boot.ci(boot_results, index = 4, type = "bca")$bca[[5]])

    }

  # }

  output <- as.data.frame(cbind(Effect = c("CDE", "NDE", "NIE",
                                           "Total Effect", "Proportion Mediated"),
                                Estimate = c(round(CDE, 5), round(NDE, 5),
                                             round(NIE, 5), round(TE, 5),
                                             round(PM, 5)),
                                `Lower 95% CI` = c(round(CI_CDE[[1]], 5),
                                                   round(CI_NDE[[1]], 5),
                                                   round(CI_NIE[[1]], 5),
                                                   round(CI_TE[[1]], 5), NA),
                                `Upper 95% CI` = c(round(CI_CDE[[2]], 5),
                                                   round(CI_NDE[[2]], 5),
                                                   round(CI_NIE[[2]], 5),
                                                   round(CI_TE[[2]], 5), NA)))
  rownames(output) <- NULL

  return(print(output))
}
