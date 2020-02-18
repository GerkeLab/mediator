
context("survival")

# import data for use in tests - roughly based off data from SAS macro --------
load("test_dat.RData")
set.seed(8675309)
test_dat$time <- round(abs(rnorm(n=nrow(test_dat)))*100)

# library(survival)

# test_dat2 <- test_dat
# test_dat2$cens <- ifelse(test_dat2$cens==0,1,0)
# write.table(test_dat2,
#             file="/Volumes/Lab_Gerke/gformulaPackage/MediationPsychMethods/test_data_survival.txt",
#             row.names = FALSE, quote=FALSE, sep="\t")

# set uniform values to test across all types ---------------------------------
data = test_dat
a = 1
a_star = 0
m = 0

test_that("coxph outcome and continuous mediator match SAS",{

  # setting options (specificed by user) --------------------------------------
  out.model = survival::coxph(survival::Surv(time,cens) ~ x + m + c + x*c,
                              data = test_dat, ties = "breslow")
  med.model = lm(c ~ x + m, data = test_dat)
  treat = "x"
  mediator = "c"
  out.reg = "coxph"
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
  expect_equal(round(CDE, 5), 0.75787)

  ## natural direct effect
  NDE <- do.call(natural_direct_effect, arg_list)
  expect_equal(round(as.numeric(NDE), 5), 0.69060)

  ## natural indirect effect
  NIE <- do.call(natural_indirect_effect, arg_list)
  expect_equal(round(as.numeric(NIE), 5), 0.99912, tolerance = 0.00001)

  ## total effect
  TE <- total_effect(NDE, NIE, out.reg)
  expect_equal(round(as.numeric(TE), 5), 0.68999)

  ## proportion mediated
  PM <- prop_mediated(NDE, NIE, out.reg, TE)
  expect_equal(round(as.numeric(PM), 9), 0.001961653)

})

test_that("coxph outcome and  binary mediator match SAS",{

  # setting options (specificed by user) --------------------------------------
  out.model = survival::coxph(survival::Surv(time,cens) ~ x + m + c + x*m,
                  data = test_dat)
  med.model = glm(m ~ x + c ,data = test_dat, family = "binomial")
  treat = "x"
  mediator = "m"
  out.reg = "coxph"
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
  expect_equal(round(CDE, 5),0.63850, tolerance = 0.001)

  ## natural direct effect
  NDE <- do.call(natural_direct_effect, arg_list)
  expect_equal(round(as.numeric(NDE), 5), 0.66819, tolerance = 0.001)

  ## natural indirect effect
  NIE <- do.call(natural_indirect_effect, arg_list)
  expect_equal(round(as.numeric(NIE), 5), 1.00014, tolerance = 0.001)

  ## total effect
  TE <- total_effect(NDE, NIE, out.reg)
  expect_equal(round(as.numeric(TE), 5), 0.66828, tolerance = 0.001)

  ## proportion mediated
  PM <- prop_mediated(NDE, NIE, out.reg, TE)
  expect_equal(round(as.numeric(PM), 9), -0.000274389, tolerance = 0.0001)
})

