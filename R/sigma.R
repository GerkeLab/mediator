comb_sigma <- function(med.model, out.model, treat, mediator_name,
                       out.reg, cnames, med.reg){

  SigmaB <- stats::vcov(med.model)
  SigmaT <- stats::vcov(out.model)
  # including 0 for variance for interaction if missing
  if(is.na(out.model$coefficients[paste0(treat, ":", mediator_name)])){
    SigmaT <- rbind(cbind(SigmaT,rep(0,nrow(SigmaT))),rep(0,nrow(SigmaT)+1))
    dimnames(SigmaT)[[1]][nrow(SigmaT)] <- paste0(treat, ":", mediator_name)
    dimnames(SigmaT)[[2]][nrow(SigmaT)] <- paste0(treat, ":", mediator_name)
  } else if (out.reg == "coxph") {
    SigmaT <- rbind(cbind(rep(0,nrow(SigmaT)),SigmaT),rep(0,nrow(SigmaT)+1))
    dimnames(SigmaT)[[1]][nrow(SigmaT)] <- "(Intercept)"
    dimnames(SigmaT)[[2]][1] <- "(Intercept)"
  } else {
    SigmaT <- SigmaT
  }
  SigmaT <- SigmaT[c("(Intercept)", treat, mediator_name, paste0(treat, ":", mediator_name), cnames),
                   c("(Intercept)", treat, mediator_name, paste0(treat, ":", mediator_name), cnames)]
  Sigma <- rbind(cbind(SigmaB, matrix(0, ncol = ncol(SigmaT), nrow = nrow(SigmaB))),
                 cbind(matrix(0, ncol = ncol(SigmaB), nrow = nrow(SigmaT)), SigmaT))
  # Sigma includes standard error only for logistic/linear and no others
  if (out.reg %in% c("logistic","coxph") & med.reg == "linear") {
    sigmaV <- stats::sigma(med.model)^2
    Sigma <- rbind(cbind(Sigma, rep(0, nrow(Sigma))),
                   c(rep(0, ncol(Sigma)), sigmaV))
  } else {
    Sigma <- Sigma
  }

}
