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

# cov pred gets betas*mean or mode for covariates
## function still needs to be cleaned up a bit
cov_pred <- function(cmeans, cmodes, treat, mediator, med.model){

  pred_data <- data.frame(t(cmeans))
  pred_data <- cbind(pred_data,cmodes)

  colnames(pred_data) <- c(names(cmeans),
                           names(cmodes))

  pred_data[treat] <- 0
  pred_data[mediator] <- 0

  # okay so TG recommended setting levels
  # if(length(med.model$xlevels) >= 1){
  #   for(i in names(med.model$xlevels)){
  #     pred_data[,i] <- factor(pred_data[,i],
  #                            levels = med.model$xlevels[[i]])
  #   }
  # }

  noise <- data[,names(pred_data)]

  pred_data <- rbind(pred_data,noise)

  tmp <- model.matrix(med.model, pred_data)

  cov_vals <- data.frame(beta = coef(med.model),
                         val = t(tmp)[,1])

  cov_vals$pred <- cov_vals$beta * cov_vals$val

  betasum <- cov_vals$pred
  names(betasum) <- row.names(cov_vals)
  betasum <- betasum[!(names(betasum) %in% c("(Intercept)",treat))]

}
