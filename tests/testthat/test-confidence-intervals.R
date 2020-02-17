context("confidence intervals")

load("test_dat.RData")

# set uniform values to test across all types ---------------------------------
data = test_dat
a = 1
a_star = 0
m = 0

test_that("continuous outcome and mediator match SAS macro",{

  # setting options (specificed by user) --------------------------------------
  out.model = glm(y ~ x + c + m + cens + x*c, data = test_dat)
  med.model = glm(c ~ x + m + cens, data = test_dat)
  treat = "x"
  mediator = "c"
  out.reg = "linear"
  med.reg = "linear"

  # calculating values to use later on ----------------------------------------
  out_vars <- if (out.reg=="coxph") names(attr(out.model$terms,"dataClasses"))[-1] else
    names(attr(out.model$terms,"dataClasses"))

  var_set <- unique(c(out_vars,
                      names(attr(med.model$terms,"dataClasses"))))
  data <- data %>% dplyr::select(var_set)

  betas <- stats::coef(med.model) # coefficients from mediation model
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

  # calculate effect estimates ------------------------------------------------

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

  # calculate gammas ----------------------------------------------------------
  gCDE <- do.call(gamma_cde, arg_list)

  ## natural direct effect
  gNDE <- do.call(gamma_nde, arg_list)

  ## natural indirect effect
  gNIE <- do.call(gamma_nie, arg_list)

  ## total effect
  gTE <- do.call(gamma_te,
                 c(arg_list, list("gNDE" = gNDE, "gNIE" = gNIE)))

  # delta method of calculating confidence intervals ---------------------------

  ## controlled direct effect
  CI_CDE <- delta_cde_nde(gCDE, Sigma, CDE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_CDE[[1]],5), -0.18535, tolerance = 0.00001)
  expect_equal(round(CI_CDE[[2]],5), 0.36890, tolerance = 0.00001)

  ## natural direct effect
  CI_NDE <- delta_cde_nde(gNDE, Sigma, NDE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_NDE[[1]], 5), -0.19302, tolerance = 0.00001)
  expect_equal(round(CI_NDE[[2]], 5), 0.36116, tolerance = 0.00001)

  ## natural indirect effect
  CI_NIE <- delta_nie(gNIE, Sigma, NIE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_NIE[[1]], 5), -0.02038, tolerance = 0.00001)
  expect_equal(round(CI_NIE[[2]], 5), 0.00083, tolerance = 0.00001)

  ## total effect
  CI_TE <- delta_te(gTE, Sigma, TE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_TE[[1]], 5), -0.20349, tolerance = 0.00001)
  expect_equal(round(CI_TE[[2]], 5), 0.35208, tolerance = 0.00001)
})

test_that("binary outcome and continuous mediator match SAS macro",{

  # setting options (specificed by user) --------------------------------------
  out.model = glm(cens ~ x + y + m + c + x*y, data = test_dat,
                  family="binomial")
  med.model = glm(y ~ x + m + c, data = test_dat)
  treat = "x"
  mediator = "y"
  out.reg = "logistic"
  med.reg = "linear"

  # calculating values to use later on ----------------------------------------
  out_vars <- if (out.reg=="coxph") names(attr(out.model$terms,"dataClasses"))[-1] else
    names(attr(out.model$terms,"dataClasses"))

  var_set <- unique(c(out_vars,
                      names(attr(med.model$terms,"dataClasses"))))
  data <- data %>% dplyr::select(var_set)

  betas <- stats::coef(med.model) # coefficients from mediation model
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
  # calculate effect estimates ------------------------------------------------

  ## controlled direct effect
  CDE <-  do.call(controlled_direct_effect, arg_list)

  ## natural direct effect
  NDE <- do.call(natural_direct_effect, arg_list)

  ## natural indirect effect
  NIE <- do.call(natural_indirect_effect, arg_list)

  ## total effect
  TE <- total_effect(NDE, NIE, out.reg)

  ## proportion mediated
  PM <- prop_mediated(NDE, NIE, out.reg, TE)

  # calculate gammas ----------------------------------------------------------
  gCDE <- do.call(gamma_cde, arg_list)

  ## natural direct effect
  gNDE <- do.call(gamma_nde, arg_list)

  ## natural indirect effect
  gNIE <- do.call(gamma_nie, arg_list)

  ## total effect
  gTE <- do.call(gamma_te,
                 c(arg_list, list("gNDE" = gNDE, "gNIE" = gNIE)))

  # delta method of calculating confidence intervals ---------------------------

  ## controlled direct effect
  CI_CDE <- delta_cde_nde(gCDE, Sigma, CDE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_CDE[[1]],5), 0.14373, tolerance = 0.00001)
  expect_equal(round(CI_CDE[[2]],5), 1.28613, tolerance = 0.00001)

  ## natural direct effect
  CI_NDE <- delta_cde_nde(gNDE, Sigma, NDE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_NDE[[1]], 5), 0.20471, tolerance = 0.00001)
  expect_equal(round(CI_NDE[[2]], 5), 2.55513, tolerance = 0.00001)

  ## natural indirect effect
  CI_NIE <- delta_nie(gNIE, Sigma, NIE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_NIE[[1]], 5), 0.81617, tolerance = 0.00001)
  expect_equal(round(CI_NIE[[2]], 5), 1.37127, tolerance = 0.00001)

  ## total effect
  CI_TE <- delta_te(gTE, Sigma, TE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_TE[[1]], 5), 0.25162, tolerance = 0.00001)
  expect_equal(round(CI_TE[[2]], 5), 2.32655, tolerance = 0.00001)

})

test_that("continuous outcome and binary mediator match SAS macro",{

  # setting options (specificed by user) --------------------------------------
  out.model = glm(y ~ x + m + c + cens + x*m, data = test_dat)
  med.model = glm(m ~ x + c + cens, data = test_dat, family = "binomial")
  treat = "x"
  mediator = "m"
  out.reg = "linear"
  med.reg = "logistic"

  # calculating values to use later on ----------------------------------------
  out_vars <- if (out.reg=="coxph") names(attr(out.model$terms,"dataClasses"))[-1] else
    names(attr(out.model$terms,"dataClasses"))

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

  # calculate effect estimates ------------------------------------------------

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

  # calculate gammas ----------------------------------------------------------
  gCDE <- do.call(gamma_cde, arg_list)

  ## natural direct effect
  gNDE <- do.call(gamma_nde, arg_list)

  ## natural indirect effect
  gNIE <- do.call(gamma_nie, arg_list)

  ## total effect
  gTE <- do.call(gamma_te,
                 c(arg_list, list("gNDE" = gNDE, "gNIE" = gNIE)))

  # delta method of calculating confidence intervals ---------------------------

  ## controlled direct effect
  CI_CDE <- delta_cde_nde(gCDE, Sigma, CDE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_CDE[[1]],5), -0.12956, tolerance = 0.00001)
  expect_equal(round(CI_CDE[[2]],5), 0.65040, tolerance = 0.00001)

  ## natural direct effect
  CI_NDE <- delta_cde_nde(gNDE, Sigma, NDE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_NDE[[1]], 5), -0.25349, tolerance = 0.00001)
  expect_equal(round(CI_NDE[[2]], 5), 0.44499, tolerance = 0.00001)

  ## natural indirect effect
  CI_NIE <- delta_nie(gNIE, Sigma, NIE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_NIE[[1]], 5), -0.03495, tolerance = 0.00001)
  expect_equal(round(CI_NIE[[2]], 5), 0.02756, tolerance = 0.00001)

  ## total effect
  CI_TE <- delta_te(gTE, Sigma, TE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_TE[[1]], 5), -0.18416, tolerance = 0.00001)
  expect_equal(round(CI_TE[[2]], 5), 0.36827, tolerance = 0.00001)

})

test_that("binary outcome and mediator match SAS macro",{

  # setting options (specificed by user) --------------------------------------
  out.model = glm(cens ~ x + m + c + y + x*m, data = test_dat, family = "binomial")
  med.model = glm(m ~ x + c + y,data = test_dat, family = "binomial")
  treat = "x"
  mediator = "m"
  out.reg = "logistic"
  med.reg = "logistic"

  # calculating values to use later on ----------------------------------------
  out_vars <- if (out.reg=="coxph") names(attr(out.model$terms,"dataClasses"))[-1] else
    names(attr(out.model$terms,"dataClasses"))

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

  # calculate effect estimates ------------------------------------------------

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

  # calculate gammas ----------------------------------------------------------
  gCDE <- do.call(gamma_cde, arg_list)

  ## natural direct effect
  gNDE <- do.call(gamma_nde, arg_list)

  ## natural indirect effect
  gNIE <- do.call(gamma_nie, arg_list)

  ## total effect
  gTE <- do.call(gamma_te,
                 c(arg_list, list("gNDE" = gNDE, "gNIE" = gNIE)))

  # delta method of calculating confidence intervals ---------------------------

  ## controlled direct effect
  CI_CDE <- delta_cde_nde(gCDE, Sigma, CDE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_CDE[[1]],5), 0.16286, tolerance = 0.00001)
  expect_equal(round(CI_CDE[[2]],5), 1.80266, tolerance = 0.00001)

  ## natural direct effect
  CI_NDE <- delta_cde_nde(gNDE, Sigma, NDE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_NDE[[1]], 5), 0.25922, tolerance = 0.00001)
  expect_equal(round(CI_NDE[[2]], 5), 1.47183, tolerance = 0.00001)

  ## natural indirect effect
  CI_NIE <- delta_nie(gNIE, Sigma, NIE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_NIE[[1]], 5), 0.90767, tolerance = 0.00001)
  expect_equal(round(CI_NIE[[2]], 5), 1.10615, tolerance = 0.00001)

  ## total effect
  CI_TE <- delta_te(gTE, Sigma, TE, a, a_star, out.reg, med.reg)
  expect_equal(round(CI_TE[[1]], 5), 0.25987, tolerance = 0.00001)
  expect_equal(round(CI_TE[[2]], 5), 1.47404, tolerance = 0.00001)

})
