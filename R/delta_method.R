delta_cde_nde <- function(gamma, Sigma, effect, a, a_star, out.reg, med.reg) {
  var <- t(gamma) %*% Sigma %*% gamma

  if (out.reg %in% c("logistic","coxph") & med.reg == "logistic") {
    CI <- exp(log(effect) + c(-1, 1) * stats::qnorm(.975) *
                as.vector(sqrt(var)))
  } else if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {
    CI <- exp(log(effect) + c(-1, 1) * stats::qnorm(.975) *
                as.vector(sqrt(var)) * (a - a_star))
  } else {
    CI <- effect + c(-1, 1) * stats::qnorm(.975) *
      as.vector(sqrt(var)) * abs(a - a_star)
  }

}

delta_nie <- function(gamma, Sigma, effect, a, a_star, out.reg, med.reg) {
  var <- t(gamma) %*% Sigma %*% gamma

  if (out.reg %in% c("logistic","coxph") & med.reg == "logistic") {
    CI <- exp(log(effect) + c(-1, 1) * stats::qnorm(.975) *
                as.vector(sqrt(var)))
  } else if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {
    CI <- exp(log(effect) + c(-1, 1) * stats::qnorm(.975) *
                as.vector(sqrt(var)) * (a - a_star))
  } else if (out.reg == "linear" & med.reg == "logistic") {
    CI <- effect + c(-1, 1) * stats::qnorm(.975) * as.vector(sqrt(var))
  } else {
    CI <- effect + c(-1, 1) * stats::qnorm(.975) * as.vector(sqrt(var)) *
      abs(a - a_star)
  }

}

delta_te <- function(gamma, Sigma, effect, a, a_star, out.reg, med.reg) {
  var <- t(gamma) %*% Sigma %*% gamma

  if (out.reg %in% c("logistic","coxph") & med.reg == "logistic") {
    CI <- exp(log(effect) + c(-1, 1) * stats::qnorm(.975) *
                as.vector(sqrt(var)))
  } else if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {
    CI <- exp(log(effect) + c(-1, 1) * stats::qnorm(.975) *
                as.vector(sqrt(var)) * (a - a_star))
  } else if (out.reg == "linear" & med.reg == "logistic") {
    CI <- effect + c(-1, 1) * stats::qnorm(.975) * as.vector(sqrt(var))
  } else {
    CI <- effect + c(-1, 1) * stats::qnorm(.975) * as.vector(sqrt(var)) *
      abs(a - a_star)
  }

}
