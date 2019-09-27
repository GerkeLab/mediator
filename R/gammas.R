gamma_cde <- function(betameans, a, a_star, m, out.reg, med.reg, ...){
  if (out.reg %in% c("logistic","coxph") & med.reg == "logistic"){
    gCDE <- c(0, 0,
              if (length(betameans)) rep(0, length(betameans)) else NA,
              0,
              (a - a_star),
              0,
              m * (a - a_star),
              if (length(betameans)) rep(0, length(betameans)) else NA)
  } else if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {
    gCDE <- c(0, 0,
              if (length(betameans)) rep(0, length(betameans)) else NA,
              0, 1, 0,
              m,
              if (length(betameans)) rep(0, length(betameans)) else NA,
              0)
  } else {
    gCDE <- c(0, 0,
              if (length(betameans)) rep(0, length(betameans)) else NA,
              0, 1, 0,
              m,
              if (length(betameans)) rep(0, length(betameans)) else NA)
  }

  gCDE <- gCDE[!is.na(gCDE)]

}

gamma_nde <- function(out.reg, med.reg, theta2, theta3, a, beta0,
                      beta1, a_star, betasum, betameans, sigmaV, ...){
  if (out.reg %in% c("logistic","coxph") & med.reg == "logistic") {
    A <- exp_sum(theta2, theta3 * a, beta0, beta1 * a_star, betasum) /
      exp_sum_1(theta2, theta3 * a, beta0, beta1 * a_star, betasum)
    B <- exp_sum(theta2, theta3 * a_star, beta0, beta1 * a_star, betasum) /
      exp_sum_1(theta2, theta3 * a_star, beta0, beta1 * a_star, betasum)
    gNDE <- c(A - B,
              a_star * (A - B),
              if (length(betameans)) t(betameans) * (A - B) else NA,
              0,
              a - a_star,
              A - B,
              a * A - a_star * B,
              if (length(betameans)) rep(0, length(betameans)) else NA)
  } else if (out.reg == "linear" & med.reg == "logistic") {
    dex <- exp_sum(beta0, beta1 * a_star, betasum)
    d1 <- (theta3 * dex * (1 + dex) - theta3 * (dex^2)) / ((1 + dex)^2)
    d2 <- (theta3 * a_star * dex * (1 + dex) - dex^2) / ((1 + dex)^2)
    d3 <- if (length(betameans)) (theta3 * t(betameans) * dex * (1 + dex) -
                                    dex^2) / ((1 + dex)^2) else NA
    d7 <- dex / (1 + dex)

    gNDE <- c(d1,
              d2,
              if (length(betameans)) d3 else NA,
              0, 1, 0,
              d7,
              if (length(betameans)) rep(0, length(betameans)) else NA)
  } else if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {
    gNDE <- c(theta3,
              theta3 * a_star,
              if (length(betameans)) theta3 * t(betameans) else NA,
              0, 1,
              theta3 * sigmaV,
              (beta0 + beta1 * a_star + betasum + theta2 * sigmaV +
                 theta3 * sigmaV * (a + a_star)), # C = c?
              if (length(betameans)) rep(0, length(betameans)) else NA,
              theta3 * theta2 + 0.5 * (theta3^2) * (a + a_star))
  } else {
    gNDE <- c(theta3,
              theta3 * a_star,
              if (length(betameans)) t(theta3 * t(betameans)) else NA,
              0, 1, 0,
              beta0 + beta1 * a_star + betasum,
              if (length(betameans)) rep(0, length(betameans)) else NA)
  }

  gNDE <- gNDE[!is.na(gNDE)]

}

gamma_nie <- function(beta0, beta1, a, betasum, a_star, out.reg,
                      med.reg, theta2, theta3, betameans, ...){

  K <- exp_sum(beta0, beta1 * a, betasum) /
    exp_sum_1(beta0, beta1 * a, betasum)

  D <- exp_sum(beta0, beta1 * a_star, betasum) /
    exp_sum_1(beta0, beta1 * a_star, betasum)

  if (out.reg %in% c("logistic","coxph") & med.reg == "logistic") {
    A <- exp_sum(theta2, theta3 * a, beta0, beta1 * a, betasum) /
      exp_sum_1(theta2, theta3 * a, beta0, beta1 * a, betasum)

    B <- exp_sum(theta2, theta3 * a, beta0, beta1 * a_star, betasum) /
      exp_sum_1(theta2, theta3 * a, beta0, beta1 * a_star, betasum)

    gNIE <- c((D + A) - (K + B),
              a_star * (D - B) + a * (A - K),
              if (length(betameans)) t(betameans) * ((D + A) - (K + B)) else NA,
              0, 0,
              A - B,
              a * (A - B),
              if (length(betameans)) rep(0, length(betameans)) else NA)

  } else if (out.reg == "linear" & med.reg == "logistic") {
    Q <- ((exp_sum(beta0, beta1 * a, betasum) *
             exp_sum_1(beta0, beta1 * a, betasum)) -
            ((exp_sum(beta0, beta1 * a, betasum))^2)) /
      (exp_sum_1(beta0, beta1 * a, betasum)^2)

    B <- ((exp_sum(beta0, beta1 * a_star, betasum) *
             exp_sum_1(beta0, beta1 * a_star, betasum)) -
            ((exp_sum(beta0, beta1 * a_star, betasum))^2)) /
      (exp_sum_1(beta0, beta1 * a_star, betasum)^2)

    gNIE <- c((theta2 + theta3 * a) * (Q - B),
              (theta2 + theta3 * a) * (a * Q - a_star * B),
              if (length(betameans)) (theta2 + theta3 * a) * t(betameans) * (Q - B) else NA,
              0, 0,
              K - D,
              a * (K - D),
              if (length(betameans)) rep(0, length(betameans)) else NA)

  } else if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {
    gNIE <- c(0,
              theta2 + theta3 * a,
              if (length(betameans)) rep(0, length(betameans)) else NA,
              0, 0,
              beta1,
              beta1 * a,
              if (length(betameans)) rep(0, length(betameans)) else NA,
              0)
  } else {
    gNIE <- c(0,
              theta2 + theta3 * a_star,
              if (length(betameans)) rep(0, length(betameans)) else NA,
              0, 0,
              beta1,
              beta1 * a_star,
              if (length(betameans)) rep(0, length(betameans)) else NA)
  }

  gNIE <- gNIE[!is.na(gNIE)]

}

gamma_te <- function(gNDE, gNIE, beta0, beta1, a, betasum, a_star,
                     theta3, theta2, betameans, out.reg, med.reg, sigmaV, ...){
  if (out.reg %in% c("logistic","coxph") & med.reg == "logistic") {

    gTE <- gNDE + gNIE

  } else if (out.reg == "linear" & med.reg == "logistic") {

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

    gTE <- c(theta3 * (a - a_star) * B + (theta2 + theta3 * a) * (Q - B),
             a_star * theta3 * (a - a_star) * B + (theta2 + theta3 * a) *
               (a * Q - a_star * B),
             if (length(betameans)) t(betameans) * theta3 * (a - a_star) *
               B + (theta2 + theta3 * a) * (Q - B) else NA,
             0,
             a - a_star,
             K - D,
             (a - a_star) * D + a * (K - D),
             if (length(betameans)) rep(0, length(betameans)) else NA)

  } else if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {

    gTE <- c(theta3,
             theta3 * (a + a_star) + theta2,
             if (length(betameans)) theta3 * t(betameans) else NA,
             0, 1,
             theta3 * sigmaV + beta1,
             beta0 + beta1 * (a + a_star) + betasum + theta2 * sigmaV +
               theta3 * sigmaV * (a^2 - a_star^2),
             if (length(betameans)) rep(0, length(betameans)) else NA,
             0.5 * (theta3^2) * (a^2 - a_star^2))

  } else {

    gTE <- c(theta3,
             theta3 * (a + a_star) + theta2,
             if (length(betameans)) theta3 %*% t(betameans) else NA,
             0, 1,
             beta1,
             beta0 + beta1 * (a + a_star) + betasum,
             if (length(betameans)) rep(0, length(betameans)) else NA)

  }

  gTE <- gTE[!is.na(gTE)]
}
