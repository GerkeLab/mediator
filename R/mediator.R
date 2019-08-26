#' Causal Mediation Analysis
#'
#' The `mediator` R function conducts mediation analysis under
#' the counterfactual framework assuming interation between the exposure
#' and mediator. Currently the function works for binary and continuous
#' outcomes and mediators.
#'
#' @param data Data set to use for analysis
#' @param out.model A fitted model object for the outcome.
#'   Can be of class 'glm','lm'.
#' @param med.model A fitted model object for the mediator.
#'   Can be of class 'glm','lm'.
#' @param treat A character string indicating the name of the
#'   treatment/exposure variable used.
#' @param mediator A character string indicating the name of the
#'   mediator variable used.
#' @param out.reg A character string indicating the type of
#'   regression used in the outcome model. Can be either "logisitic" or
#'   "linear".
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
mediator <- function(data = dat,
                     out.model = glm(cens ~ x + m + c + x*m,
                                     data = dat,
                                     family = "binomial"),
                     med.model = glm(m ~ x + c,
                                     data = dat,
                                     family = "binomial"),
                     treat = "x",
                     mediator = "m",
                     out.reg = "logistic",
                     med.reg = "logistic",
                     a = 1,
                     a_star = 0,
                     m = 0,
                     boot_rep = 0){

  # calculating covariate values to use later on
  betas <- coef(med.model) # coefficients from mediation model
  cmeans <- apply(data, 2, function(x) mean(as.numeric(x), na.rm = TRUE)) # mean value for all values
  betameans <- cmeans[which(names(cmeans) %in%
                              names(betas)[!(names(betas) %in%
                                               c("(Intercept)", treat))])] # subset to only covariates
  betameans <- betameans[match(names(betas)[!(names(betas) %in% c("(Intercept)", treat))],
                               names(betameans))] # put in order to match coefficients
  betasum <- betameans %*% betas[names(betas)[!(names(betas) %in% c("(Intercept)", treat))]] # mean * coefficient from model

  # get covariate names
  covmeans <- cmeans[which(names(cmeans) %in% names(betas)[!(names(betas) %in% c("(Intercept)", treat))])]
  cnames <- names(covmeans)

  # Covariance matrix for standar errors
  SigmaB <- stats::vcov(med.model)
  SigmaT <- stats::vcov(out.model)
  # including 0 for variance for interaction if missing
  if(is.na(out.model$coefficients[paste0(treat, ":", mediator)])){
    SigmaT <- rbind(cbind(SigmaT,rep(0,nrow(SigmaT))),rep(0,nrow(SigmaT)))
    dimnames(SigmaT)[[1]][nrow(SigmaT)] <- paste0(treat, ":", mediator)
    dimnames(SigmaT)[[2]][nrow(SigmaT)] <- paste0(treat, ":", mediator)
  } else {
    SigmaT <- SigmaT
  }
  SigmaT <- SigmaT[c("(Intercept)", treat, mediator, paste0(treat, ":", mediator), cnames),
                   c("(Intercept)", treat, mediator, paste0(treat, ":", mediator), cnames)]
  Sigma <- rbind(cbind(SigmaB, matrix(0, ncol = ncol(SigmaT), nrow = nrow(SigmaB))),
                 cbind(matrix(0, ncol = ncol(SigmaB), nrow = nrow(SigmaT)), SigmaT))
  # Sigma includes standard error only for logistic/linear and no others
  if (out.reg == "logistic" & med.reg == "linear") {
     sigmaV <- sigma(med.model)^2
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

  ##### ----------------------------------------------------------------- #####
  # Calculate effect estimates and confidence intervals (delta method) --------
  ##### ----------------------------------------------------------------- #####
  if (out.reg == "logistic" & med.reg == "logistic") {

    # calculate effect estimates
    ## controlled direct effect
    CDE <- as.numeric(out.model$coefficients[treat] * (a - a_star) +
                        out.model$coefficients[paste0(treat, ":", mediator)] * m *
                        (a - a_star))
    CDE <- exp(CDE)
    ## natural direct effect
    NDEnum <- exp(out.model$coefficients[treat] * a) *
      (1 + exp(out.model$coefficients[mediator] +
                 out.model$coefficients[paste0(treat, ":", mediator)] *
                 a + med.model$coefficients["(Intercept)"] +
                 med.model$coefficients[treat] * a_star + betasum))
    NDEden <- exp(out.model$coefficients[treat] * a_star) *
      (1 + exp(out.model$coefficients[mediator] +
                 out.model$coefficients[paste0(treat, ":", mediator)] *
                 a_star + med.model$coefficients["(Intercept)"] +
                 med.model$coefficients[treat] * a_star + betasum))
    NDE <- NDEnum / NDEden
    rm(NDEnum, NDEden)
    ## natural indirect effect
    NIEnum <- (1 + exp(med.model$coefficients["(Intercept)"] +
                         med.model$coefficients[treat] * a_star + betasum)) *
      (1 + exp(out.model$coefficients[mediator] +
                 out.model$coefficients[paste0(treat, ":", mediator)] *
                 a + med.model$coefficients["(Intercept)"] +
                 med.model$coefficients[treat] * a + betasum))
    NIEden <- (1 + exp(med.model$coefficients["(Intercept)"] +
                         med.model$coefficients[treat] * a + betasum)) *
      (1 + exp(out.model$coefficients[mediator] +
                 out.model$coefficients[paste0(treat, ":", mediator)] *
                 a + med.model$coefficients["(Intercept)"] +
                 med.model$coefficients[treat] * a_star + betasum))
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
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA),
                0,
                (a - a_star),
                0,
                m * (a - a_star),
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA))
      gCDE <- gCDE[!is.na(gCDE)]
      ## natural direct effect
      A <- (exp(out.model$coefficients[mediator] +
                  out.model$coefficients[paste0(treat, ":", mediator)] *
                  a + med.model$coefficients["(Intercept)"] +
                  med.model$coefficients[treat] * a_star + betasum)) /
        (1 + exp(out.model$coefficients[mediator] +
                   out.model$coefficients[paste0(treat, ":", mediator)] *
                   a + med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] * a_star + betasum))
      B <- (exp(out.model$coefficients[mediator] +
                  out.model$coefficients[paste0(treat, ":", mediator)] *
                  a_star + med.model$coefficients["(Intercept)"] +
                  med.model$coefficients[treat] * a_star + betasum)) /
        (1 + exp(out.model$coefficients[mediator] +
                   out.model$coefficients[paste0(treat, ":", mediator)] *
                   a_star + med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] * a_star + betasum))
      # A_total <- (exp(out.model$coefficients[mediator] + out.model$coefficients[paste0(treat,":",mediator)]*a_star + med.model$coefficients["(Intercept)"] + med.model$coefficients[treat]*a + betasum))/
      #   (1 + exp(out.model$coefficients[mediator] + out.model$coefficients[paste0(treat,":",mediator)]*a_star + med.model$coefficients["(Intercept)"] + med.model$coefficients[treat]*a + betasum))
      # B_total <- (exp(out.model$coefficients[mediator] + out.model$coefficients[paste0(treat,":",mediator)]*a + med.model$coefficients["(Intercept)"] + med.model$coefficients[treat]*a + betasum))/
      #   (1 + exp(out.model$coefficients[mediator] + out.model$coefficients[paste0(treat,":",mediator)]*a + med.model$coefficients["(Intercept)"] + med.model$coefficients[treat]*a + betasum))
      gNDE <- c(A - B,
                a_star * (A - B),
                ifelse(length(betameans) > 0, t(betameans) * (A - B), NA), # problem child
                0,
                a - a_star,
                A - B,
                a * A - a_star * B,
                ifelse(length(betameans) > 0, rep(0,length(betameans)), NA))
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
      A <- (exp(out.model$coefficients[mediator] +
                  out.model$coefficients[paste0(treat, ":", mediator)] *
                  a + med.model$coefficients["(Intercept)"] +
                  med.model$coefficients[treat] * a + betasum)) /
        (1 + exp(out.model$coefficients[mediator] +
                   out.model$coefficients[paste0(treat, ":", mediator)] *
                   a + med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] * a + betasum))

      B <- (exp(out.model$coefficients[mediator] +
                  out.model$coefficients[paste0(treat, ":", mediator)] *
                  a + med.model$coefficients["(Intercept)"] +
                  med.model$coefficients[treat] * a_star + betasum)) /
        (1 + exp(out.model$coefficients[mediator] +
                   out.model$coefficients[paste0(treat, ":", mediator)] *
                   a + med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] * a_star + betasum))

      K <- (exp(med.model$coefficients["(Intercept)"] +
                  med.model$coefficients[treat] * a + betasum)) /
        (1 + exp(med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] * a + betasum))

      D <- (exp(med.model$coefficients["(Intercept)"] +
                  med.model$coefficients[treat] * a_star + betasum)) /
        (1 + exp(med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] * a_star + betasum))

      gNIE <- c((D + A) - (K + B),
                a_star * (D - B) + a * (A - K),
                ifelse(length(betameans) > 0, t(betameans) * ((D + A) - (K + B)), NA), # problem child
                0, 0,
                A - B,
                a * (A - B),
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA))

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
    CDE <-  as.numeric(out.model$coefficients[treat] * (a - a_star) +
                         out.model$coefficients[paste0(treat, ":", mediator)] *
                         m * (a - a_star))

    ## natural direct effect
    NDE <- out.model$coefficients[treat] * (a - a_star) +
             (out.model$coefficients[paste0(treat, ":", mediator)] *
                (a - a_star)) * ((exp(med.model$coefficients["(Intercept)"] +
                                        med.model$coefficients[treat] *
                                        a_star + betasum)) /
                                   (1 + exp(med.model$coefficients["(Intercept)"] +
                                              med.model$coefficients[treat] *
                                              a_star + betasum)))

    ## natural indirect effect
    NIE <- (out.model$coefficients[mediator] +
              out.model$coefficients[paste0(treat, ":", mediator)] * a) *
      ((exp(med.model$coefficients["(Intercept)"] +
              med.model$coefficients[treat] * a + betasum) /
          (1 + exp(med.model$coefficients["(Intercept)"] +
                     med.model$coefficients[treat] * a + betasum))) -
         (exp(med.model$coefficients["(Intercept)"] +
                med.model$coefficients[treat] * a_star + betasum) /
            (1 + exp(med.model$coefficients["(Intercept)"] +
                       med.model$coefficients[treat] * a_star + betasum))))

    ## total effect
    TE <- NDE + NIE

    ## proportion mediated
    PM <- NIE / TE

    # calculate GAMMA for SE
    if (boot_rep == 0) {

      ## controlled direct effect
      gCDE <- c(0, 0,
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA),
                0, 1, 0,
                m,
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA))
      gCDE <- gCDE[!is.na(gCDE)]

      ## natural direct effect
      dex <- exp(med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] * a_star + betasum)
      d1 <- (out.model$coefficients[paste0(treat, ":", mediator)] * dex *
               (1 + dex) - out.model$coefficients[paste0(treat, ":", mediator)] *
               (dex^2)) / ((1 + dex)^2)
      d2 <- (out.model$coefficients[paste0(treat, ":", mediator)] *
               a_star * dex * (1 + dex) - dex^2) / ((1 + dex)^2)
      d3 <- ifelse(length(betameans) > 0,
                   (out.model$coefficients[paste0(treat, ":", mediator)] *
                      t(betameans) * dex * (1 + dex) - dex^2) /
                     ((1 + dex)^2),
                   NA)
      d7 <- dex / (1 + dex)

      gNDE <- c(d1,
                d2,
                ifelse(length(betameans) > 0, d3, NA),
                0, 1, 0,
                d7,
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA))
      gNDE <- gNDE[!is.na(gNDE)]

      rm(d1, d2, d3, d7)

      ## natural indirect effect
      Q <- ((exp(med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] * a + betasum) *
               (1 + exp(med.model$coefficients["(Intercept)"] +
                          med.model$coefficients[treat] * a + betasum))) -
              ((exp(med.model$coefficients["(Intercept)"] +
                      med.model$coefficients[treat] * a + betasum))^2)) /
        ((1 + exp(med.model$coefficients["(Intercept)"] +
                    med.model$coefficients[treat] * a + betasum))^2)
      B <- ((exp(med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] * a_star + betasum) *
               (1 + exp(med.model$coefficients["(Intercept)"] +
                          med.model$coefficients[treat] * a_star + betasum))) -
              ((exp(med.model$coefficients["(Intercept)"] +
                      med.model$coefficients[treat] * a_star + betasum))^2)) /
        ((1 + exp(med.model$coefficients["(Intercept)"] +
                    med.model$coefficients[treat] * a_star + betasum))^2)
      K <- exp(med.model$coefficients["(Intercept)"] +
                 med.model$coefficients[treat] * a + betasum) /
        (1 + exp(med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] * a + betasum))
      D <- exp(med.model$coefficients["(Intercept)"] +
                 med.model$coefficients[treat] * a_star + betasum) /
        (1 + exp(med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] * a_star + betasum))

      gNIE <- c((out.model$coefficients[mediator] +
                   out.model$coefficients[paste0(treat, ":", mediator)] * a) *
                  (Q - B),
                (out.model$coefficients[mediator] +
                   out.model$coefficients[paste0(treat, ":", mediator)] * a) *
                  (a * Q - a_star * B),
                ifelse(length(betameans) > 0,
                       (out.model$coefficients[mediator] +
                          out.model$coefficients[paste0(treat, ":", mediator)] *
                          a) * t(betameans) * (Q - B),
                       NA),
                0, 0,
                K - D,
                a * (K - D),
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA))
      gNIE <- gNIE[!is.na(gNIE)]

      ## total effect - exchange A for Q 07/01/2019
      gTE <- c(out.model$coefficients[paste0(treat, ":", mediator)] * (a - a_star) * B + (out.model$coefficients[mediator] + out.model$coefficients[paste0(treat, ":", mediator)] * a) * (Q - B),
               a_star * out.model$coefficients[paste0(treat, ":", mediator)] *
                 (a - a_star) * B + (out.model$coefficients[mediator] +
                                       out.model$coefficients[paste0(treat, ":", mediator)] *
                                       a) * (a * Q - a_star * B),
               ifelse(length(betameans) > 0,
                      t(betameans) * out.model$coefficients[paste0(treat, ":", mediator)] *
                        (a - a_star) * B + (out.model$coefficients[mediator] +
                                              out.model$coefficients[paste0(treat, ":", mediator)] *
                                              a) * (Q - B),
                      NA), # problem child
               0,
               a - a_star,
               K - D,
               (a - a_star) * D + a * (K - D),
               ifelse(length(betameans) > 0, rep(0, length(betameans)), NA))
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

  } else if (out.reg == "logistic" & med.reg == "linear") {

    # calculate effect estimates

    ## controlled direct effect
    CDE <-  as.numeric(out.model$coefficients[treat] * (a - a_star) +
                         out.model$coefficients[paste0(treat, ":", mediator)] *
                         m * (a - a_star))
    CDE <- exp(CDE)

    ## natural direct effect
    NDE <- (out.model$coefficients[treat] +
              out.model$coefficients[paste0(treat, ":", mediator)] *
              (med.model$coefficients["(Intercept)"] +
                 med.model$coefficients[treat] * a_star + betasum +
                 (out.model$coefficients[mediator] * sigmaV))) *
      (a - a_star) + (0.5 * (out.model$coefficients[paste0(treat, ":", mediator)]^2) *
                        sigmaV) * (a^2 - a_star^2)
    NDE <- exp(NDE)

    ## natural indirect effect
    NIE <- (out.model$coefficients[mediator] *
              med.model$coefficients[treat] +
              out.model$coefficients[paste0(treat, ":", mediator)] *
              med.model$coefficients[treat] * a) * (a - a_star)
    NIE <- exp(NIE)

    ## total effect
    TE <- NDE * NIE
    ## proportion mediated
    PM <- (NDE * (NIE - 1)) / (NDE * NIE - 1)

    # calculate GAMMA for SE
    if (boot_rep == 0) {
      ## controlled direct effect
      gCDE <- c(0, 0,
                ifelse(length(betameans) > 0,
                       rep(0, length(betameans)),
                       NA),
                0, 1, 0, m,
                ifelse(length(betameans) > 0,
                       rep(0, length(betameans)),
                       NA),
                0)
      gCDE <- gCDE[!is.na(gCDE)]
      ## natural direct effect
      gNDE <- c(out.model$coefficients[paste0(treat, ":", mediator)],
                out.model$coefficients[paste0(treat, ":", mediator)] *
                  a_star,
                ifelse(length(betameans) > 0,
                       out.model$coefficients[paste0(treat, ":", mediator)] *
                         t(betameans),
                       NA), # problem child
                0, 1,
                out.model$coefficients[paste0(treat, ":", mediator)] * sigmaV,
                (med.model$coefficients["(Intercept)"] +
                   med.model$coefficients[treat] *
                   a_star + betasum + out.model$coefficients[mediator] *
                   sigmaV + out.model$coefficients[paste0(treat, ":", mediator)] *
                   sigmaV * (a + a_star)), # C = c?
                ifelse(length(betameans) > 0,
                       rep(0, length(betameans)),
                       NA),
                out.model$coefficients[paste0(treat, ":", mediator)] *
                  out.model$coefficients[mediator] + 0.5 *
                  (out.model$coefficients[paste0(treat, ":", mediator)]^2) *
                  (a + a_star))
      gNDE <- gNDE[!is.na(gNDE)]

      ## natural indirect effect
      gNIE <- c(0,
                out.model$coefficients[mediator] +
                  out.model$coefficients[paste0(treat, ":", mediator)] * a,
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA),
                0, 0,
                med.model$coefficients[treat],
                med.model$coefficients[treat] * a,
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA), 0)
      gNIE <- gNIE[!is.na(gNIE)]

      ## total effect
      gTE <- c(out.model$coefficients[paste0(treat, ":", mediator)],
               out.model$coefficients[paste0(treat, ":", mediator)] *
                 (a + a_star) + out.model$coefficients[mediator],
               ifelse(length(betameans) > 0,
                      out.model$coefficients[paste0(treat, ":", mediator)] *
                        t(betameans),
                      NA), # problem child
               0, 1,
               out.model$coefficients[paste0(treat, ":", mediator)] *
                 sigmaV + med.model$coefficients[treat],
               med.model$coefficients["(Intercept)"] +
                 med.model$coefficients[treat] * (a + a_star) + betasum +
                 out.model$coefficients[mediator] * sigmaV +
                 out.model$coefficients[paste0(treat, ":", mediator)] *
                 sigmaV * (a^2 - a_star^2),
               ifelse(length(betameans) > 0, rep(0, length(betameans)), NA),
               0.5 * (out.model$coefficients[paste0(treat, ":", mediator)]^2) *
                 (a^2 - a_star^2))
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
    CDE <- as.numeric(out.model$coefficients[treat] * (a - a_star) +
                        out.model$coefficients[paste0(treat, ":", mediator)] *
                        m * (a - a_star))

    ## natural direct effect
    NDE <- (out.model$coefficients[treat] +
              out.model$coefficients[paste0(treat, ":", mediator)] *
              med.model$coefficients["(Intercept)"] +
              out.model$coefficients[paste0(treat, ":", mediator)] *
              med.model$coefficients[treat] * a_star +
              out.model$coefficients[paste0(treat, ":", mediator)] * betasum) *
      (a - a_star)

    ## natural indirect effect
    NIE <- (out.model$coefficients[mediator] * med.model$coefficients[treat] +
              out.model$coefficients[paste0(treat, ":", mediator)] *
              med.model$coefficients[treat] * a) * (a - a_star)

    ## total effect
    TE <- NDE + NIE

    ## proportion mediated
    PM <- NIE / TE

    # calculate GAMMA for SE
    if (boot_rep == 0) {

      ## controlled direct effect
      gCDE <- c(0, 0,
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA),
                0, 1, 0,
                m,
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA))
      gCDE <- gCDE[!is.na(gCDE)]

      ## natural direct effect
      gNDE <- c(out.model$coefficients[paste0(treat, ":", mediator)],
                out.model$coefficients[paste0(treat, ":", mediator)] * a_star,
                ifelse(length(betameans) > 0,
                       t(out.model$coefficients[paste0(treat, ":", mediator)] *
                           t(betameans)),
                       NA), # problem child
                0, 1, 0,
                med.model$coefficients["(Intercept)"] +
                  med.model$coefficients[treat] * a_star + betasum,
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA))
      gNDE <- gNDE[!is.na(gNDE)]

      ## pure natural indirect effect
      ### for total NIE - substitute a and a*
      gNIE <- c(0,
                out.model$coefficients[mediator] +
                  out.model$coefficients[paste0(treat, ":", mediator)] * a_star,
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA),
                0, 0,
                med.model$coefficients[treat],
                med.model$coefficients[treat] * a_star,
                ifelse(length(betameans) > 0, rep(0, length(betameans)), NA))
      gNIE <- gNIE[!is.na(gNIE)]

      ## total effects
      gTE <- c(out.model$coefficients[paste0(treat, ":", mediator)],
               out.model$coefficients[paste0(treat, ":", mediator)] *
                 (a + a_star) + out.model$coefficients[mediator],
               ifelse(length(betameans) > 0,
                      out.model$coefficients[paste0(treat, ":", mediator)] %*%
                        t(betameans),
                      NA), # problem child : c' and C' = ?
               0, 1,
               med.model$coefficients[treat],
               med.model$coefficients["(Intercept)"] +
                 med.model$coefficients[treat] * (a + a_star) + betasum,
               ifelse(length(betameans) > 0, rep(0, length(betameans)), NA))
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

        # calculating covariate values to use later on
        betas <- coef(med) # coefficients from mediation model
        cmeans <- apply(d, 2, function(x) mean(as.numeric(x), na.rm = TRUE)) # mean value for all values
        betameans <- cmeans[which(names(cmeans) %in%
                                    names(betas)[!(names(betas) %in%
                                                     c("(Intercept)", treat))])] # subset to only covariates
        betameans <- betameans[match(names(betas)[!(names(betas) %in%
                                                      c("(Intercept)", treat))],
                                     names(betameans))] # put in order to match coefficients
        betasum <- betameans %*% betas[names(betas)[!(names(betas) %in%
                                                        c("(Intercept)", treat))]] # mean * coefficient from model

        # get covariate names
        covmeans <- cmeans[which(names(cmeans) %in%
                                   names(betas)[!(names(betas) %in%
                                                    c("(Intercept)", treat))])]
        cnames <- names(covmeans)

        # Covariance matrix for standar errors
        SigmaB <- stats::vcov(med)
        SigmaT <- stats::vcov(out)
        SigmaT <- SigmaT[c("(Intercept)", treat, mediator,
                           paste0(treat, ":", mediator), cnames),
                         c("(Intercept)", treat, mediator,
                           paste0(treat, ":", mediator), cnames)]

        Sigma <- rbind(cbind(SigmaB, matrix(0, ncol = ncol(SigmaT),
                                            nrow = nrow(SigmaB))),
                       cbind(matrix(0, ncol = ncol(SigmaB),
                                    nrow = nrow(SigmaT)), SigmaT))

        # Sigma includes standard error only for logistic/linear and no others
        if (out.reg == "logistic" & med.reg == "linear") {

          sigmaV <- stats::sigma(med)^2

          Sigma <- rbind(cbind(Sigma, rep(0, nrow(Sigma))),
                         c(rep(0, ncol(Sigma)), sigmaV))


        } else {
          Sigma <- Sigma
        }

        rm(SigmaB, SigmaT)

        if (out.reg == "logistic" & med.reg == "logistic") {

          # calculate effect estimates

          ## controlled direct effect
          CDE <- as.numeric(out$coefficients[treat] *
                              (a - a_star) +
                              out$coefficients[paste0(treat, ":", mediator)] *
                              m * (a - a_star))
          CDE <- exp(CDE)

          ## natural direct effect
          NDEnum <- exp(out$coefficients[treat] * a) *
            (1 + exp(out$coefficients[mediator] +
                       out$coefficients[paste0(treat, ":", mediator)] *
                       a + med$coefficients["(Intercept)"] +
                       med$coefficients[treat] * a_star + betasum))
          NDEden <- exp(out$coefficients[treat] * a_star) *
            (1 + exp(out$coefficients[mediator] +
                       out$coefficients[paste0(treat, ":", mediator)] *
                       a_star + med$coefficients["(Intercept)"] +
                       med$coefficients[treat] * a_star + betasum))
          NDE <- NDEnum / NDEden

          rm(NDEnum, NDEden)

          ## natural indirect effect
          NIEnum <- (1 + exp(med$coefficients["(Intercept)"] +
                               med$coefficients[treat] * a_star + betasum)) *
            (1 + exp(out$coefficients[mediator] +
                       out$coefficients[paste0(treat, ":", mediator)] *
                       a + med$coefficients["(Intercept)"] +
                       med$coefficients[treat] * a + betasum))
          NIEden <- (1 + exp(med$coefficients["(Intercept)"] +
                               med$coefficients[treat] * a + betasum)) *
            (1 + exp(out$coefficients[mediator] +
                       out$coefficients[paste0(treat, ":", mediator)] *
                       a + med$coefficients["(Intercept)"] +
                       med$coefficients[treat] * a_star + betasum))

          NIE <- as.vector(NIEnum / NIEden)

          rm(NIEnum, NIEden)

          ## total effect
          TE <- NDE * NIE

        } else if (out.reg == "linear" & med.reg == "logistic") {

          # calculate effect estimates

          ## controlled direct effect
          CDE <-  as.numeric(out$coefficients[treat] * (a - a_star) +
                               out$coefficients[paste0(treat, ":", mediator)] *
                               m * (a - a_star))

          ## natural direct effect
          NDE <- out$coefficients[treat] * (a - a_star) +
            (out$coefficients[paste0(treat, ":", mediator)] * (a - a_star)) *
            ((exp(med$coefficients["(Intercept)"] + med$coefficients[treat] *
                    a_star + betasum)) /
               (1 + exp(med$coefficients["(Intercept)"] +
                          med$coefficients[treat] * a_star + betasum)))

          ## natural indirect effect
          NIE <- (out$coefficients[mediator] +
                    out$coefficients[paste0(treat, ":", mediator)] * a) *
            ((exp(med$coefficients["(Intercept)"] + med$coefficients[treat] *
                    a + betasum) /
                (1 + exp(med$coefficients["(Intercept)"] +
                           med$coefficients[treat] * a + betasum))) -
               (exp(med$coefficients["(Intercept)"] + med$coefficients[treat] *
                      a_star + betasum) /
                  (1 + exp(med$coefficients["(Intercept)"] +
                             med$coefficients[treat] * a_star + betasum))))

          ## total effect
          TE <- NDE + NIE

        } else if (out.reg == "logistic" & med.reg == "linear") {

          # calculate effect estimates

          ## controlled direct effect
          CDE <-  as.numeric(out$coefficients[treat] * (a - a_star) +
                               out$coefficients[paste0(treat, ":", mediator)] *
                               m * (a - a_star))
          CDE <- exp(CDE)

          ## natural direct effect
          NDE <- (out$coefficients[treat] +
                    out$coefficients[paste0(treat, ":", mediator)] *
                    (med$coefficients["(Intercept)"] + med$coefficients[treat] *
                       a_star + betasum + (out$coefficients[mediator] *
                                             sigmaV))) * (a - a_star) +
            (0.5 * (out$coefficients[paste0(treat, ":", mediator)]^2) * sigmaV) *
            (a^2 - a_star^2)
          NDE <- exp(NDE)

          ## natural indirect effect
          NIE <- (out$coefficients[mediator] * med$coefficients[treat] +
                    out$coefficients[paste0(treat, ":", mediator)] *
                    med$coefficients[treat] * a) * (a - a_star)
          NIE <- exp(NIE)

          ## total effect
          TE <- NDE * NIE

        } else {

          # calculate effect estimates

          ## controlled direct effect
          CDE <- as.numeric(out$coefficients[treat] * (a - a_star) +
                              out$coefficients[paste0(treat, ":", mediator)] *
                              m * (a - a_star))

          ## natural direct effect
          NDE <- (out$coefficients[treat] +
                    out$coefficients[paste0(treat, ":", mediator)] *
                    med$coefficients["(Intercept)"] +
                    out$coefficients[paste0(treat, ":", mediator)] *
                    med$coefficients[treat] * a_star +
                    out$coefficients[paste0(treat, ":", mediator)] * betasum) *
            (a - a_star)

          ## natural indirect effect
          NIE <- (out$coefficients[mediator] * med$coefficients[treat] +
                    out$coefficients[paste0(treat, ":", mediator)] *
                    med$coefficients[treat] * a) * (a - a_star)

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
