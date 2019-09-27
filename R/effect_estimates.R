controlled_direct_effect <- function(theta1, a, a_star, theta3, m, out.reg, ...) {

  CDE <- as.numeric(theta1 * (a - a_star) + theta3 * m * (a - a_star))

  if (out.reg %in% c("logistic","coxph")){

    CDE <- exp(CDE)

  } else {

    CDE

  }
}

total_effect <- function(NDE, NIE, out.reg) {

  if (out.reg %in% c("logistic","coxph")){

    TE <- NDE * NIE

  } else {

    TE <- NDE + NIE

  }
}

prop_mediated <- function(NDE, NIE, out.reg, TE) {

  if (out.reg %in% c("logistic","coxph")) {

    PM <- (NDE * (NIE - 1)) / (NDE * NIE - 1)

  } else {

    PM <- NIE / TE

  }
}

natural_direct_effect <- function(theta1, a, theta2, theta3, beta0,
                                  beta1, a_star, betasum, sigmaV,
                                  out.reg, med.reg, ...) {

  if (out.reg %in% c("logistic","coxph") & med.reg == "logistic") {

    NDEnum <- exp(theta1 * a) * exp_sum_1(theta2, theta3 * a, beta0,
                                          beta1 * a_star, betasum)

    NDEden <- exp(theta1 * a_star) * exp_sum_1(theta2, theta3 * a_star,
                                               beta0, beta1 * a_star,
                                               betasum)

    NDE <- NDEnum / NDEden

  } else if (out.reg == "linear" & med.reg == "logistic") {

    NDE <- theta1 * (a - a_star) + (theta3 * (a - a_star)) *
      (exp_sum(beta0, beta1 * a_star, betasum) /
         exp_sum_1(beta0, beta1 * a_star, betasum))

  } else if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {

    NDE <- (theta1 + theta3 * (beta0 + beta1 * a_star + betasum +
                                 (theta2 * sigmaV))) * (a - a_star) +
      (0.5 * (theta3^2) * sigmaV) * (a^2 - a_star^2)
    NDE <- exp(NDE)

  } else {

    NDE <- (theta1 + theta3 * beta0 + theta3 * beta1 * a_star + theta3 *
              betasum) * (a - a_star)

  }
}

natural_indirect_effect <- function(out.reg, med.reg, beta0, beta1, a_star,
                                    betasum, theta2, theta3, a, ...){
  if (out.reg %in% c("logistic","coxph") & med.reg == "logistic") {

    NIEnum <- exp_sum_1(beta0, beta1 * a_star, betasum) *
      exp_sum_1(theta2, theta3 * a, beta0, beta1 * a, betasum)

    NIEden <- exp_sum_1(beta0, beta1 * a, betasum) *
      exp_sum_1(theta2, theta3 * a, beta0, beta1 * a_star, betasum)

    NIE <- as.vector(NIEnum / NIEden)

  } else if (out.reg == "linear" & med.reg == "logistic") {

    NIE <- (theta2 + theta3 * a) * ((exp_sum(beta0, beta1 * a, betasum) /
                                       exp_sum_1(beta0, beta1 * a, betasum)) -
                                      (exp_sum(beta0, beta1 * a_star, betasum) /
                                         exp_sum_1(beta0, beta1 * a_star, betasum)))

  } else if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {

    NIE <- (theta2 * beta1 + theta3 * beta1 * a) * (a - a_star)
    NIE <- exp(NIE)

  } else {

    NIE <- (theta2 * beta1 + theta3 * beta1 * a) * (a - a_star)

  }
}
