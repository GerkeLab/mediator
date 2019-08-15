context("multiple covariates")

# import data for use in tests - roughly based off data from SAS macro --------
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
  betas <- coef(med.model) # coefficients from mediation model
  cmeans <- apply(data, 2, function(x) mean(as.numeric(x), na.rm = TRUE)) # mean value for all values
  betameans <- cmeans[which(names(cmeans) %in%
                              names(betas)[!(names(betas) %in%
                                               c("(Intercept)", treat))])] # subset to only covariates
  betameans <- betameans[match(names(betas)[!(names(betas) %in% c("(Intercept)", treat))],
                               names(betameans))] # put in order to match coefficients
  betasum <- betameans %*% betas[names(betas)[!(names(betas) %in% c("(Intercept)", treat))]] # mean * coefficient from model

  # get covariate names -------------------------------------------------------
  covmeans <- cmeans[which(names(cmeans) %in% names(betas)[!(names(betas) %in% c("(Intercept)", treat))])]
  cnames <- names(covmeans)

  # Covariance matrix for standar errors --------------------------------------
  # set mediator covariance as 0
  SigmaB <- vcov(med.model)
  SigmaT <- vcov(out.model)
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
  rm(SigmaB, SigmaT)

  # calculate effect estimates ------------------------------------------------

  ## controlled direct effect
  CDE <- as.numeric(out.model$coefficients[treat] * (a - a_star) +
                      out.model$coefficients[paste0(treat, ":", mediator)] *
                      m * (a - a_star))
  expect_equal(round(CDE, 6), 0.091774)

  ## natural direct effect
  NDE <- (out.model$coefficients[treat] +
            out.model$coefficients[paste0(treat, ":", mediator)] *
            med.model$coefficients["(Intercept)"] +
            out.model$coefficients[paste0(treat, ":", mediator)] *
            med.model$coefficients[treat] * a_star +
            out.model$coefficients[paste0(treat, ":", mediator)] * betasum) *
    (a - a_star)
  expect_equal(round(as.numeric(NDE), 6), 0.084069)

  ## natural indirect effect
  NIE <- (out.model$coefficients[mediator] * med.model$coefficients[treat] +
            out.model$coefficients[paste0(treat, ":", mediator)] *
            med.model$coefficients[treat] * a) * (a - a_star)
  expect_equal(round(as.numeric(NIE), 6), -0.009772)

  ## total effect
  TE <- NDE + NIE
  expect_equal(round(as.numeric(TE), 6), 0.074297)

  ## proportion mediated
  PM <- NIE / TE
  expect_equal(round(as.numeric(PM), 5), -0.13153

  )

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
  betas <- coef(med.model) # coefficients from mediation model
  cmeans <- apply(data, 2, function(x) mean(as.numeric(x), na.rm = TRUE)) # mean value for all values
  betameans <- cmeans[which(names(cmeans) %in%
                              names(betas)[!(names(betas) %in%
                                               c("(Intercept)", treat))])] # subset to only covariates
  betameans <- betameans[match(names(betas)[!(names(betas) %in% c("(Intercept)", treat))],
                               names(betameans))] # put in order to match coefficients
  betasum <- betameans %*% betas[names(betas)[!(names(betas) %in% c("(Intercept)", treat))]] # mean * coefficient from model

  # get covariate names -------------------------------------------------------
  covmeans <- cmeans[which(names(cmeans) %in% names(betas)[!(names(betas) %in% c("(Intercept)", treat))])]
  cnames <- names(covmeans)

  # Covariance matrix for standar errors --------------------------------------
  # set mediator covariance as 0
  SigmaB <- vcov(med.model)
  SigmaT <- vcov(out.model)
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
  rm(SigmaB, SigmaT)
  # Sigma includes standard error only for logistic/linear and no others
  if (out.reg == "logistic" & med.reg == "linear") {
    sigmaV <- sigma(med.model)^2
    Sigma <- rbind(cbind(Sigma, rep(0, nrow(Sigma))),
                   c(rep(0, ncol(Sigma)), sigmaV))

  } else {
    Sigma <- Sigma
  }

  # calculate effect estimates ------------------------------------------------

  ## controlled direct effect
  CDE <-  as.numeric(out.model$coefficients[treat] * (a - a_star) +
                       out.model$coefficients[paste0(treat, ":", mediator)] *
                       m * (a - a_star))
  CDE <- exp(CDE)
  expect_equal(round(CDE, 5), 0.42995)

  ## natural direct effect
  NDE <- (out.model$coefficients[treat] +
            out.model$coefficients[paste0(treat, ":", mediator)] *
            (med.model$coefficients["(Intercept)"] +
               med.model$coefficients[treat] * a_star + betasum +
               (out.model$coefficients[mediator] * sigmaV))) *
    (a - a_star) + (0.5 * (out.model$coefficients[paste0(treat, ":", mediator)]^2) *
                      sigmaV) * (a^2 - a_star^2)
  NDE <- exp(NDE)
  expect_equal(round(as.numeric(NDE), 5), 0.72323)

  ## natural indirect effect
  NIE <- (out.model$coefficients[mediator] *
            med.model$coefficients[treat] +
            out.model$coefficients[paste0(treat, ":", mediator)] *
            med.model$coefficients[treat] * a) * (a - a_star)
  NIE <- exp(NIE)
  expect_equal(round(as.numeric(NIE), 5), 1.05791, tolerance = 0.00001)

  ## total effect
  TE <- NDE * NIE
  expect_equal(round(as.numeric(TE), 5), 0.76512)

  ## proportion mediated
  PM <- (NDE * (NIE - 1)) / (NDE * NIE - 1)
  expect_equal(round(as.numeric(PM), 5), -0.17833)

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
  betas <- coef(med.model) # coefficients from mediation model
  cmeans <- apply(data, 2, function(x) mean(as.numeric(x), na.rm = TRUE)) # mean value for all values
  betameans <- cmeans[which(names(cmeans) %in%
                              names(betas)[!(names(betas) %in%
                                               c("(Intercept)", treat))])] # subset to only covariates
  betameans <- betameans[match(names(betas)[!(names(betas) %in% c("(Intercept)", treat))],
                               names(betameans))] # put in order to match coefficients
  betasum <- betameans %*% betas[names(betas)[!(names(betas) %in% c("(Intercept)", treat))]] # mean * coefficient from model

  # get covariate names -------------------------------------------------------
  covmeans <- cmeans[which(names(cmeans) %in% names(betas)[!(names(betas) %in% c("(Intercept)", treat))])]
  cnames <- names(covmeans)

  # Covariance matrix for standar errors --------------------------------------
  # set mediator covariance as 0
  SigmaB <- vcov(med.model)
  SigmaT <- vcov(out.model)
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
  rm(SigmaB, SigmaT)
  # Sigma includes standard error only for logistic/linear and no others
  if (out.reg == "logistic" & med.reg == "linear") {
    sigmaV <- sigma(med.model)^2
    Sigma <- rbind(cbind(Sigma, rep(0, nrow(Sigma))),
                   c(rep(0, ncol(Sigma)), sigmaV))

  } else {
    Sigma <- Sigma
  }

  # calculate effect estimates ------------------------------------------------

  ## controlled direct effect
  CDE <-  as.numeric(out.model$coefficients[treat] * (a - a_star) +
                       out.model$coefficients[paste0(treat, ":", mediator)] *
                       m * (a - a_star))
  expect_equal(round(CDE, 6), 0.26042, tolerance = 0.00001)

  ## natural direct effect
  NDE <- out.model$coefficients[treat] * (a - a_star) +
    (out.model$coefficients[paste0(treat, ":", mediator)] *
       (a - a_star)) * ((exp(med.model$coefficients["(Intercept)"] +
                               med.model$coefficients[treat] *
                               a_star + betasum)) /
                          (1 + exp(med.model$coefficients["(Intercept)"] +
                                     med.model$coefficients[treat] *
                                     a_star + betasum)))
  expect_equal(round(as.numeric(NDE), 6), 0.09575, tolerance = 0.00001)

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
  expect_equal(round(as.numeric(NIE), 6), -0.00369, tolerance = 0.00001)

  ## total effect
  TE <- NDE + NIE
  expect_equal(round(as.numeric(TE), 6), 0.09206, tolerance = 0.00001)

  ## proportion mediated
  PM <- NIE / TE
  expect_equal(round(as.numeric(PM), 5), -0.040117, tolerance = 0.00001)
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
  betas <- coef(med.model) # coefficients from mediation model
  cmeans <- apply(data, 2, function(x) mean(as.numeric(x), na.rm = TRUE)) # mean value for all values
  betameans <- cmeans[which(names(cmeans) %in%
                              names(betas)[!(names(betas) %in%
                                               c("(Intercept)", treat))])] # subset to only covariates
  betameans <- betameans[match(names(betas)[!(names(betas) %in% c("(Intercept)", treat))],
                               names(betameans))] # put in order to match coefficients
  betasum <- betameans %*% betas[names(betas)[!(names(betas) %in% c("(Intercept)", treat))]] # mean * coefficient from model

  # get covariate names -------------------------------------------------------
  covmeans <- cmeans[which(names(cmeans) %in% names(betas)[!(names(betas) %in% c("(Intercept)", treat))])]
  cnames <- names(covmeans)

  # Covariance matrix for standar errors --------------------------------------
  # set mediator covariance as 0
  SigmaB <- vcov(med.model)
  SigmaT <- vcov(out.model)
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
  rm(SigmaB, SigmaT)
  # Sigma includes standard error only for logistic/linear and no others
  if (out.reg == "logistic" & med.reg == "linear") {
    sigmaV <- sigma(med.model)^2
    Sigma <- rbind(cbind(Sigma, rep(0, nrow(Sigma))),
                   c(rep(0, ncol(Sigma)), sigmaV))

  } else {
    Sigma <- Sigma
  }

  # calculate effect estimates ------------------------------------------------

  ## controlled direct effect
  CDE <- as.numeric(out.model$coefficients[treat] * (a - a_star) +
                      out.model$coefficients[paste0(treat, ":", mediator)] * m *
                      (a - a_star))
  CDE <- exp(CDE)
  expect_equal(round(CDE, 5), 0.54184, tolerance = 0.00001)

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
  expect_equal(round(as.numeric(NDE), 5), 0.61768, tolerance = 0.00001)

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
  expect_equal(round(as.numeric(NIE), 5), 1.00201, tolerance = 0.00001)

  ## total effect
  TE <- NDE * NIE
  expect_equal(round(as.numeric(TE), 5), 0.61892, tolerance = 0.00001)

  ## proportion mediated
  PM <- (NDE * (NIE - 1)) / (NDE * NIE - 1)
  expect_equal(round(as.numeric(PM), 6), -0.003255, tolerance = 0.00001)
})
