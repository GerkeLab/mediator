# exp_sum simply sums all numeric arguments and exponentiates
exp_sum <- function(...) {
  exp(sum(c(...)))
}
# exp_sum(1, 3, 4) returns result of exp(1 + 3 + 4)

# exp_sum_1 sums all numeric arguments, exponentiates, then adds 1
exp_sum_1 <- function(...) {
  1 + exp(sum(c(...)))
}
# exp_sum_1(1, 3, 4) returns result of 1 + exp(1 + 3 + 4)

# cov_pred gets betas*mean or betas*mode for covariates
## function still needs to be cleaned up a bit
cov_pred <- function(cmeans, cmodes, treat, mediator, med.model, data){

  # combine means and modes into data frame -----------------------------------
  pred_data <- data.frame(t(cmeans))
  pred_data <- cbind(pred_data,t(cmodes))

  # set treatment and mediator to 0 -------------------------------------------
  pred_data[treat] <- 0
  pred_data[mediator] <- 0

  # set factors to levels from data -------------------------------------------
  if(length(med.model$xlevels) >= 1){
    for(i in names(med.model$xlevels)){
      pred_data[,i] <- factor(pred_data[,i],
                              levels = med.model$xlevels[[i]])
    }
  }

  # add noise from original data ----------------------------------------------
  noise <- data[,names(pred_data)]
  pred_data <- rbind(pred_data,noise)

  # get model matrix and multiply by coef from model --------------------------
  ## subset to only those needed - not intercept or treatment
  tmp <- stats::model.matrix(med.model, data = pred_data)
  cov_vals <- data.frame(beta = stats::coef(med.model),
                         val = t(tmp)[,1])
  drop <- c("(Intercept)",treat)
  cov_vals <- cov_vals[!(rownames(cov_vals) %in% drop),]
  betasum <- cov_vals$beta * cov_vals$val
  names(betasum) <- row.names(cov_vals)

  betamean <- cov_vals$val
  names(betamean) <- row.names(cov_vals)
  # betasum <- betasum[!(names(betasum) %in% c("(Intercept)",treat))]

  return(list(betasum = betasum, betamean = betamean))

}
