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



