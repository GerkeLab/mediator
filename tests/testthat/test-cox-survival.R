context("survival")

# import data for use in tests - roughly based off data from SAS macro --------
load("test_dat.RData")
set.seed(8675309)
test_dat$time <- round(abs(rnorm(n=nrow(test_dat)))*100)

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
  } else if (out.reg=="coxph") {
    SigmaT <- suppressWarnings(rbind(cbind(rep(0,nrow(SigmaT)),SigmaT),rep(0,nrow(SigmaT))))
    dimnames(SigmaT)[[1]][nrow(SigmaT)] <- "(Intercept)"
    dimnames(SigmaT)[[2]][1] <- "(Intercept)"
  } else {
    SigmaT <- SigmaT
  }
  SigmaT <- SigmaT[c("(Intercept)", treat, mediator, paste0(treat, ":", mediator), cnames),
                   c("(Intercept)", treat, mediator, paste0(treat, ":", mediator), cnames)]
  Sigma <- rbind(cbind(SigmaB, matrix(0, ncol = ncol(SigmaT), nrow = nrow(SigmaB))),
                 cbind(matrix(0, ncol = ncol(SigmaB), nrow = nrow(SigmaT)), SigmaT))
  rm(SigmaB, SigmaT)
  # Sigma includes standard error only for logistic/coxph outcome and linear mediator
  if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {
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
  expect_equal(round(CDE, 5), 0.75787)

  ## natural direct effect
  NDE <- (out.model$coefficients[treat] +
            out.model$coefficients[paste0(treat, ":", mediator)] *
            (med.model$coefficients["(Intercept)"] +
               med.model$coefficients[treat] * a_star + betasum +
               (out.model$coefficients[mediator] * sigmaV))) *
    (a - a_star) + (0.5 * (out.model$coefficients[paste0(treat, ":", mediator)]^2) *
                      sigmaV) * (a^2 - a_star^2)
  NDE <- exp(NDE)
  expect_equal(round(as.numeric(NDE), 5), 0.69060)

  ## natural indirect effect
  NIE <- (out.model$coefficients[mediator] *
            med.model$coefficients[treat] +
            out.model$coefficients[paste0(treat, ":", mediator)] *
            med.model$coefficients[treat] * a) * (a - a_star)
  NIE <- exp(NIE)
  expect_equal(round(as.numeric(NIE), 5), 0.99912, tolerance = 0.00001)

  ## total effect
  TE <- NDE * NIE
  expect_equal(round(as.numeric(TE), 5), 0.68999)

  ## proportion mediated
  PM <- (NDE * (NIE - 1)) / (NDE * NIE - 1)
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
  } else if (out.reg=="coxph") {
    SigmaT <- suppressWarnings(rbind(cbind(rep(0,nrow(SigmaT)),SigmaT),rep(0,nrow(SigmaT))))
    dimnames(SigmaT)[[1]][nrow(SigmaT)] <- "(Intercept)"
    dimnames(SigmaT)[[2]][1] <- "(Intercept)"
  } else {
    SigmaT <- SigmaT
  }
  SigmaT <- SigmaT[c("(Intercept)", treat, mediator, paste0(treat, ":", mediator), cnames),
                   c("(Intercept)", treat, mediator, paste0(treat, ":", mediator), cnames)]
  Sigma <- rbind(cbind(SigmaB, matrix(0, ncol = ncol(SigmaT), nrow = nrow(SigmaB))),
                 cbind(matrix(0, ncol = ncol(SigmaB), nrow = nrow(SigmaT)), SigmaT))
  rm(SigmaB, SigmaT)
  # Sigma includes standard error only for logistic/linear and no others
  if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {
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
  expect_equal(round(CDE, 5),0.63850, tolerance = 0.001)

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
  expect_equal(round(as.numeric(NDE), 5), 0.66819, tolerance = 0.001)

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
  expect_equal(round(as.numeric(NIE), 5), 1.00014, tolerance = 0.001)

  ## total effect
  TE <- NDE * NIE
  expect_equal(round(as.numeric(TE), 5), 0.66828, tolerance = 0.001)

  ## proportion mediated
  PM <- (NDE * (NIE - 1)) / (NDE * NIE - 1)
  expect_equal(round(as.numeric(PM), 9), -0.000274389, tolerance = 0.0001)
})

