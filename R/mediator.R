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
  cmeans <- data %>%
    dplyr::select_if(is.numeric) %>%
    purrr::map_dbl(~mean(.x, na.rm = TRUE)) # mean value for all numeric values
  cmodes <- data %>%
    dplyr::select_if(purrr::negate(is.numeric)) %>%
    purrr::map_chr(~{
      ux <- unique(.x)
      ux[which.max(tabulate(match(.x, ux)))]
      })

  betas <- stats::coef(med.model) # coefficients from mediation model
  # betameans <- cmeans[which(names(cmeans) %in%
  #                             names(betas)[!(names(betas) %in%
  #                                              c("(Intercept)", treat))])] # subset to only covariates
  # betameans <- betameans[match(names(betas)[!(names(betas) %in% c("(Intercept)", treat))],
  #                              names(betameans))] # put in order to match coefficients
  # betasum <- betameans %*% betas[names(betas)[!(names(betas) %in% c("(Intercept)", treat))]] # mean * coefficient from model
  beta_info <- cov_pred(cmeans, cmodes, treat, mediator, med.model, data)
  betasum <- sum(beta_info$betasum, na.rm=TRUE)
  betameans <- beta_info$betamean

  # get covariate names
  # covmeans <- cmeans[which(names(cmeans) %in% names(betas)[!(names(betas) %in% c("(Intercept)", treat))])]
  # cnames <- names(covmeans)
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
    CDE <- as.numeric(theta1 * (a - a_star) + theta3 * m * (a - a_star))
    CDE <- exp(CDE)
    ## natural direct effect
    NDEnum <- exp(theta1 * a) * exp_sum_1(theta2, theta3 * a, beta0,
                                          beta1 * a_star, betasum)
    NDEden <- exp(theta1 * a_star) * exp_sum_1(theta2, theta3 * a_star,
                                               beta0, beta1 * a_star,
                                               betasum)
    NDE <- NDEnum / NDEden
    rm(NDEnum, NDEden)
    ## natural indirect effect
    NIEnum <- exp_sum_1(beta0, beta1 * a_star, betasum) *
      exp_sum_1(theta2, theta3 * a, beta0, beta1 * a, betasum)
    NIEden <- exp_sum_1(beta0, beta1 * a, betasum) *
      exp_sum_1(theta2, theta3 * a, beta0, beta1 * a_star, betasum)
    NIE <- as.vector(NIEnum / NIEden)
    rm(NIEnum, NIEden)
    ## total effect
    TE <- NDE * NIE
    ## proportion mediated
    PM <- (NDE * (NIE - 1)) / (NDE * NIE - 1)

    # calculate GAMMA for SE
    if (boot_rep == 0) {
      ## controlled direct effect
      gCDE <- c(0, 0,
                if (length(betameans)) rep(0, length(betameans)) else NA,
                0,
                (a - a_star),
                0,
                m * (a - a_star),
                if (length(betameans)) rep(0, length(betameans)) else NA)
      gCDE <- gCDE[!is.na(gCDE)]
      ## natural direct effect
      A <- exp_sum(theta2, theta3 * a, beta0, beta1 * a_star, betasum) /
        exp_sum_1(theta2, theta3 * a, beta0, beta1 * a_star, betasum)
      B <- exp_sum(theta2, theta3 * a_star, beta0, beta1 * a_star, betasum) /
        exp_sum_1(theta2, theta3 * a_star, beta0, beta1 * a_star, betasum)
      # A_total <- (exp(out.model$coefficients[mediator] + out.model$coefficients[paste0(treat,":",mediator)]*a_star + med.model$coefficients["(Intercept)"] + med.model$coefficients[treat]*a + betasum))/
      #   (1 + exp(out.model$coefficients[mediator] + out.model$coefficients[paste0(treat,":",mediator)]*a_star + med.model$coefficients["(Intercept)"] + med.model$coefficients[treat]*a + betasum))
      # B_total <- (exp(out.model$coefficients[mediator] + out.model$coefficients[paste0(treat,":",mediator)]*a + med.model$coefficients["(Intercept)"] + med.model$coefficients[treat]*a + betasum))/
      #   (1 + exp(out.model$coefficients[mediator] + out.model$coefficients[paste0(treat,":",mediator)]*a + med.model$coefficients["(Intercept)"] + med.model$coefficients[treat]*a + betasum))
      gNDE <- c(A - B,
                a_star * (A - B),
                if (length(betameans)) t(betameans) * (A - B) else NA, # problem child
                0,
                a - a_star,
                A - B,
                a * A - a_star * B,
                if (length(betameans)) rep(0, length(betameans)) else NA)
      gNDE <- gNDE[!is.na(gNDE)]
      # ## pure NDE needed for calculating CI for TE
      # gNDE_pure <- c(A_pure-B_pure,
      #                a_star*(A_pure-B_pure),
      #                ifelse(length(betameans)>0,t(betameans)*(A_pure-B_pure),NA), # problem child
      #                0,
      #                a-a_star,
      #                A_pure-B_pure,
      #                a*A_pure-a_star*B_pure,
      #                ifelse(length(betameans)>0,rep(0,length(betameans)),NA))
      # gNDE_pure <- gNDE_pure[!is.na(gNDE_pure)]
      rm(A, B)
      ## natural indirect effect
      A <- exp_sum(theta2, theta3 * a, beta0, beta1 * a, betasum) /
        exp_sum_1(theta2, theta3 * a, beta0, beta1 * a, betasum)

      B <- exp_sum(theta2, theta3 * a, beta0, beta1 * a_star, betasum) /
        exp_sum_1(theta2, theta3 * a, beta0, beta1 * a_star, betasum)

      K <- exp_sum(beta0, beta1 * a, betasum) /
        exp_sum_1(beta0, beta1 * a, betasum)

      D <- exp_sum(beta0, beta1 * a_star, betasum) /
        exp_sum_1(beta0, beta1 * a_star, betasum)

      gNIE <- c((D + A) - (K + B),
                a_star * (D - B) + a * (A - K),
                if (length(betameans)) t(betameans) * ((D + A) - (K + B)) else NA, # problem child
                0, 0,
                A - B,
                a * (A - B),
                if (length(betameans)) rep(0, length(betameans)) else NA)

      gNIE <- gNIE[!is.na(gNIE)]

      rm(A, B, K, D)

      ## total effect
      gTE <- gNDE + gNIE # problem child

      # delta method of calculating confidence intervals

      ## controlled direct effect
      varCDE <- t(gCDE) %*% Sigma %*% gCDE
      CI_CDE <- exp(log(CDE) + c(-1, 1) * stats::qnorm(.975) *
                      as.vector(sqrt(varCDE)))

      ## natural direct effect
      varNDE <- t(gNDE) %*% Sigma %*% gNDE
      CI_NDE <- exp(as.vector(log(NDE)) + c(-1, 1) * stats::qnorm(.975) *
                      as.vector(sqrt(varNDE)))

      ## natural indirect effect
      varNIE <- t(gNIE) %*% Sigma %*% gNIE
      CI_NIE <- exp(log(NIE) + c(-1, 1) * stats::qnorm(.975) *
                      as.vector(sqrt(varNIE)))

      ## total effect
      varTE <- t(gTE) %*% Sigma %*% gTE
      CI_TE <- exp(as.vector(log(TE)) + c(-1, 1) * stats::qnorm(.975) *
                     as.vector(sqrt(varTE)))

    }

  } else if (out.reg == "linear" & med.reg == "logistic") {

    # calculate effect estimates

    ## controlled direct effect
    CDE <-  as.numeric(theta1 * (a - a_star) + theta3 * m * (a - a_star))

    ## natural direct effect
    NDE <- theta1 * (a - a_star) + (theta3 * (a - a_star)) *
      (exp_sum(beta0, beta1 * a_star, betasum) /
         exp_sum_1(beta0, beta1 * a_star, betasum))

    ## natural indirect effect
    NIE <- (theta2 + theta3 * a) * ((exp_sum(beta0, beta1 * a, betasum) /
                                       exp_sum_1(beta0, beta1 * a, betasum)) -
                                      (exp_sum(beta0, beta1 * a_star, betasum) /
                                         exp_sum_1(beta0, beta1 * a_star, betasum)))

    ## total effect
    TE <- NDE + NIE

    ## proportion mediated
    PM <- NIE / TE

    # calculate GAMMA for SE
    if (boot_rep == 0) {

      ## controlled direct effect
      gCDE <- c(0, 0,
                if (length(betameans)) rep(0, length(betameans)) else NA,
                0, 1, 0,
                m,
                if (length(betameans)) rep(0, length(betameans)) else NA)
      gCDE <- gCDE[!is.na(gCDE)]

      ## natural direct effect
      dex <- exp_sum(beta0, beta1 * a_star, betasum)
      d1 <- (theta3 * dex * (1 + dex) - theta3 * (dex^2)) / ((1 + dex)^2)
      d2 <- (theta3 * a_star * dex * (1 + dex) - dex^2) / ((1 + dex)^2)
      d3 <- if (length(betameans)) (theta3 * t(betameans) * dex * (1 + dex) - dex^2) /
                     ((1 + dex)^2) else NA
      d7 <- dex / (1 + dex)

      gNDE <- c(d1,
                d2,
                if (length(betameans)) d3 else NA,
                0, 1, 0,
                d7,
                if (length(betameans)) rep(0, length(betameans)) else NA)
      gNDE <- gNDE[!is.na(gNDE)]

      rm(d1, d2, d3, d7)

      ## natural indirect effect
      Q <- ((exp_sum(beta0, beta1 * a, betasum) *
               exp_sum_1(beta0, beta1 * a, betasum)) -
              ((exp_sum(beta0, beta1 * a, betasum))^2)) /
        (exp_sum_1(beta0, beta1 * a, betasum)^2)
      B <- ((exp_sum(beta0, beta1 * a_star, betasum) *
               exp_sum_1(beta0, beta1 * a_star, betasum)) -
              ((exp_sum(beta0, beta1 * a_star, betasum))^2)) /
        (exp_sum_1(beta0, beta1 * a_star, betasum)^2)
      K <- exp_sum(beta0, beta1 * a, betasum) /
        exp_sum_1(beta0, beta1 * a, betasum)
      D <- exp_sum(beta0, beta1 * a_star, betasum) /
        exp_sum_1(beta0, beta1 * a_star, betasum)

      gNIE <- c((theta2 + theta3 * a) * (Q - B),
                (theta2 + theta3 * a) * (a * Q - a_star * B),
                if (length(betameans)) (theta2 + theta3 * a) * t(betameans) * (Q - B) else NA,
                0, 0,
                K - D,
                a * (K - D),
                if (length(betameans)) rep(0, length(betameans)) else NA)
      gNIE <- gNIE[!is.na(gNIE)]

      ## total effect - exchange A for Q 07/01/2019
      gTE <- c(theta3 * (a - a_star) * B + (theta2 + theta3 * a) * (Q - B),
               a_star * theta3 * (a - a_star) * B + (theta2 + theta3 * a) *
                 (a * Q - a_star * B),
               if (length(betameans)) t(betameans) * theta3 * (a - a_star) *
                 B + (theta2 + theta3 * a) * (Q - B) else NA, # problem child
               0,
               a - a_star,
               K - D,
               (a - a_star) * D + a * (K - D),
               if (length(betameans)) rep(0, length(betameans)) else NA)
      gTE <- gTE[!is.na(gTE)]

      # delta method of calculating confidence intervals

      ## controlled direct effect
      varCDE <- t(gCDE) %*% Sigma %*% gCDE
      CI_CDE <- CDE + c(-1, 1) * stats::qnorm(.975) * as.vector(sqrt(varCDE)) *
        abs(a - a_star)

      ## natural direct effect
      varNDE <- t(gNDE) %*% Sigma %*% gNDE
      CI_NDE <- as.vector(NDE) + c(-1, 1) * stats::qnorm(.975) *
        as.vector(sqrt(varNDE)) * abs(a - a_star)

      ## natural indirect effect
      varNIE <- t(gNIE) %*% Sigma %*% gNIE
      CI_NIE <- as.vector(NIE) + c(-1, 1) * stats::qnorm(.975) *
        as.vector(sqrt(varNIE))

      ## total effect
      varTE <- t(gTE) %*% Sigma %*% gTE
      CI_TE <- as.vector(TE) + c(-1, 1) * stats::qnorm(.975) *
        as.vector(sqrt(varTE))
    }

  } else if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {

    # calculate effect estimates

    ## controlled direct effect
    CDE <-  as.numeric(theta1 * (a - a_star) + theta3 * m * (a - a_star))
    CDE <- exp(CDE)

    ## natural direct effect
    NDE <- (theta1 + theta3 * (beta0 + beta1 * a_star + betasum +
                                 (theta2 * sigmaV))) * (a - a_star) +
      (0.5 * (theta3^2) * sigmaV) * (a^2 - a_star^2)
    NDE <- exp(NDE)

    ## natural indirect effect
    NIE <- (theta2 * beta1 + theta3 * beta1 * a) * (a - a_star)
    NIE <- exp(NIE)

    ## total effect
    TE <- NDE * NIE
    ## proportion mediated
    PM <- (NDE * (NIE - 1)) / (NDE * NIE - 1)

    # calculate GAMMA for SE
    if (boot_rep == 0) {
      ## controlled direct effect
      gCDE <- c(0, 0,
                if (length(betameans)) rep(0, length(betameans)) else NA,
                0, 1, 0, m,
                if (length(betameans)) rep(0, length(betameans)) else NA,
                0)
      gCDE <- gCDE[!is.na(gCDE)]
      ## natural direct effect
      gNDE <- c(theta3,
                theta3 * a_star,
                if (length(betameans)) theta3 * t(betameans) else NA, # problem child
                0, 1,
                theta3 * sigmaV,
                (beta0 + beta1 * a_star + betasum + theta2 * sigmaV +
                   theta3 * sigmaV * (a + a_star)), # C = c?
                if (length(betameans)) rep(0, length(betameans)) else NA,
                theta3 * theta2 + 0.5 * (theta3^2) * (a + a_star))
      gNDE <- gNDE[!is.na(gNDE)]

      ## natural indirect effect
      gNIE <- c(0,
                theta2 + theta3 * a,
                if (length(betameans)) rep(0, length(betameans)) else NA,
                0, 0,
                beta1,
                beta1 * a,
                if (length(betameans)) rep(0, length(betameans)) else NA,
                0)
      gNIE <- gNIE[!is.na(gNIE)]

      ## total effect
      gTE <- c(theta3,
               theta3 * (a + a_star) + theta2,
               if (length(betameans)) theta3 * t(betameans) else NA, # problem child
               0, 1,
               theta3 * sigmaV + beta1,
               beta0 + beta1 * (a + a_star) + betasum + theta2 * sigmaV +
                 theta3 * sigmaV * (a^2 - a_star^2),
               if (length(betameans)) rep(0, length(betameans)) else NA,
               0.5 * (theta3^2) * (a^2 - a_star^2))
      gTE <- gTE[!is.na(gTE)]

      # calculate CI using delta method

      ## controlled direct effect
      varCDE <- t(gCDE) %*% Sigma %*% gCDE
      CI_CDE <- log(CDE) + c(-1, 1) * stats::qnorm(.975) *
        as.vector(sqrt(varCDE)) * (a - a_star)
      CI_CDE <- exp(CI_CDE)

      ## natural direct effect
      varNDE <- t(gNDE) %*% Sigma %*% gNDE
      CI_NDE <- as.vector(log(NDE)) + c(-1, 1) * stats::qnorm(.975) *
        as.vector(sqrt(varNDE)) * (a - a_star)
      CI_NDE <- exp(CI_NDE)

      ## natural indirect effect
      varNIE <- t(gNIE) %*% Sigma %*% gNIE
      CI_NIE <- exp(as.vector(log(NIE)) + c(-1, 1) * stats::qnorm(.975) *
                      as.vector(sqrt(varNIE)) * (a - a_star))

      ## total effect
      varTE <- t(gTE) %*% Sigma %*% gTE
      CI_TE <- exp(as.vector(log(TE)) + c(-1, 1) * stats::qnorm(.975) *
                     as.vector(sqrt(varTE)) * (a - a_star))

    }

  } else {

    # calculate effect estimates

    ## controlled direct effect
    CDE <- as.numeric(theta1 * (a - a_star) + theta3 * m * (a - a_star))

    ## natural direct effect
    NDE <- (theta1 + theta3 * beta0 + theta3 * beta1 * a_star + theta3 *
              betasum) * (a - a_star)

    ## natural indirect effect
    NIE <- (theta2 * beta1 + theta3 * beta1 * a) * (a - a_star)

    ## total effect
    TE <- NDE + NIE

    ## proportion mediated
    PM <- NIE / TE

    # calculate GAMMA for SE
    if (boot_rep == 0) {

      ## controlled direct effect
      gCDE <- c(0, 0,
                if (length(betameans)) rep(0, length(betameans)) else NA,
                0, 1, 0,
                m,
                if (length(betameans)) rep(0, length(betameans)) else NA)
      gCDE <- gCDE[!is.na(gCDE)]

      ## natural direct effect
      gNDE <- c(theta3,
                theta3 * a_star,
                if (length(betameans)) t(theta3 * t(betameans)) else NA, # problem child
                0, 1, 0,
                beta0 + beta1 * a_star + betasum,
                if (length(betameans)) rep(0, length(betameans)) else NA)
      gNDE <- gNDE[!is.na(gNDE)]

      ## pure natural indirect effect
      ### for total NIE - substitute a and a*
      gNIE <- c(0,
                theta2 + theta3 * a_star,
                if (length(betameans)) rep(0, length(betameans)) else NA,
                0, 0,
                beta1,
                beta1 * a_star,
                if (length(betameans)) rep(0, length(betameans)) else NA)
      gNIE <- gNIE[!is.na(gNIE)]

      ## total effects
      gTE <- c(theta3,
               theta3 * (a + a_star) + theta2,
               if (length(betameans)) theta3 %*% t(betameans) else NA, # problem child : c' and C' = ?
               0, 1,
               beta1,
               beta0 + beta1 * (a + a_star) + betasum,
               if (length(betameans)) rep(0, length(betameans)) else NA)
      gTE <- gTE[!is.na(gTE)]

      # delta method of calculating confidence intervals

      ## controlled direct effect
      varCDE <- t(gCDE) %*% Sigma %*% gCDE
      CI_CDE <- CDE + c(-1, 1) * stats::qnorm(.975) *
        as.vector(sqrt(varCDE)) * abs(a - a_star)

      ## natural direct effect - still being problematic even without covars
      varNDE <- t(gNDE) %*% Sigma %*% gNDE
      CI_NDE <- as.vector(NDE) + c(-1, 1) * stats::qnorm(.975) *
        as.vector(sqrt(varNDE)) * abs(a - a_star)

      ## natural indirect effect
      varNIE <- t(gNIE) %*% Sigma %*% gNIE
      CI_NIE <- NIE + c(-1, 1) * stats::qnorm(.975) *
        as.vector(sqrt(varNIE)) * abs(a - a_star)

      ## total effect
      varTE <- t(gTE) %*% Sigma %*% gTE
      CI_TE <- as.vector(TE) + c(-1, 1) * stats::qnorm(.975) *
        as.vector(sqrt(varTE)) * abs(a - a_star)

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

        cmeans <- d %>%
          dplyr::select_if(is.numeric) %>%
          purrr::map_dbl(~mean(.x, na.rm = TRUE)) # mean value for all numeric values
        cmodes <- d %>%
          dplyr::select_if(purrr::negate(is.numeric)) %>%
          purrr::map_chr(~{
            ux <- unique(.x)
            ux[which.max(tabulate(match(.x, ux)))]
          })

        betas <- stats::coef(med) # coefficients from mediation model
        beta_info <- cov_pred(cmeans, cmodes, treat, mediator, med, d)
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
          CDE <- as.numeric(theta1 * (a - a_star) + theta3 * m * (a - a_star))
          CDE <- exp(CDE)

          ## natural direct effect
          NDEnum <- exp(theta1 * a) * exp_sum_1(theta2, theta3 * a, beta0,
                                                beta1 * a_star, betasum)
          NDEden <- exp(theta1 * a_star) * exp_sum_1(theta2, theta3 * a_star,
                                                     beta0, beta1 * a_star,
                                                     betasum)
          NDE <- NDEnum / NDEden

          rm(NDEnum, NDEden)

          ## natural indirect effect
          NIEnum <- exp_sum_1(beta0, beta1 * a_star, betasum) *
            exp_sum_1(theta2, theta3 * a, beta0, beta1 * a, betasum)
          NIEden <- exp_sum_1(beta0, beta1 * a, betasum) *
            exp_sum_1(theta2, theta3 * a, beta0, beta1 * a_star, betasum)

          NIE <- as.vector(NIEnum / NIEden)

          rm(NIEnum, NIEden)

          ## total effect
          TE <- NDE * NIE

        } else if (out.reg == "linear" & med.reg == "logistic") {

          # calculate effect estimates

          ## controlled direct effect
          CDE <-  as.numeric(theta1 * (a - a_star) + theta3 * m * (a - a_star))

          ## natural direct effect
          NDE <- theta1 * (a - a_star) + (theta3 * (a - a_star)) *
            ((exp_sum(beta0, beta1 * a_star, betasum)) /
               exp_sum_1(beta0, beta1 * a_star, betasum))

          ## natural indirect effect
          NIE <- (theta2 + theta3 * a) * ((exp_sum(beta0, beta1 * a, betasum) /
                                             exp_sum_1(beta0, beta1 * a, betasum)) -
                                            (exp_sum(beta0, beta1 * a_star, betasum) /
                                               exp_sum_1(beta0, beta1 * a_star, betasum)))

          ## total effect
          TE <- NDE + NIE

        } else if (out.reg == "logistic" & med.reg == "linear") {

          # calculate effect estimates

          ## controlled direct effect
          CDE <-  as.numeric(theta1 * (a - a_star) + theta3 * m * (a - a_star))
          CDE <- exp(CDE)

          ## natural direct effect
          NDE <- (theta1 + theta3 * (beta0 + beta1 * a_star + betasum +
                                       (theta2 * sigmaV))) * (a - a_star) +
            (0.5 * (theta3^2) * sigmaV) * (a^2 - a_star^2)
          NDE <- exp(NDE)

          ## natural indirect effect
          NIE <- (theta2 * beta1 + theta3 * beta1 * a) * (a - a_star)
          NIE <- exp(NIE)

          ## total effect
          TE <- NDE * NIE

        } else {

          # calculate effect estimates

          ## controlled direct effect
          CDE <- as.numeric(theta1 * (a - a_star) + theta3 * m * (a - a_star))

          ## natural direct effect
          NDE <- (theta1 + theta3 * beta0 + theta3 * beta1 * a_star + theta3 *
                    betasum) * (a - a_star)

          ## natural indirect effect
          NIE <- (theta2 * beta1 + theta3 * beta1 * a) * (a - a_star)

          ## total effect
          TE <- NDE + NIE

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
  # return(print(CI_CDE))
}
