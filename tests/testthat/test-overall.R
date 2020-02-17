context("overall")

load("test_dat.RData")

test_that("R function matches SAS macro",{

  # first test ------------------------------------------------------
  first_test <- mediator(data = test_dat,
                         out.model = glm(y ~ x + c + m + cens + x*c,
                                         data = test_dat),
                         med.model = glm(c ~ x + m + cens, data = test_dat),
                         treat = "x",
                         a = 1, a_star = 0, m = 0, boot_rep = 0)

  expect_equal(first_test, tibble(Effect = c("CDE", "NDE", "NIE",
                                             "Total Effect",
                                             "Proportion Mediated"),
                                  Estimate = c(0.09177, 0.08407, -0.00977,
                                               0.07430, -0.13153),
                                  `Lower 95% CI` = c(-0.18535, -0.19302, -0.02038,
                                                     -0.20349, NA),
                                  `Upper 95% CI` = c(0.36890, 0.36116, 0.00083,
                                                     0.35208, NA)),
               tolerance = 0.00001)


   # second test ----------------------------------------------------
   second_test <- mediator(data = test_dat,
                           out.model = glm(cens ~ x + y + m + c + x*y,
                                           data = test_dat, family="binomial"),
                           med.model = glm(y ~ x + m + c, data = test_dat),
                           treat = "x",
                           a = 1, a_star = 0, m = 0, boot_rep = 0)

   expect_equal(second_test, tibble(Effect = c("CDE", "NDE", "NIE",
                                               "Total Effect",
                                               "Proportion Mediated"),
                                    Estimate = c(0.42995, 0.72323, 1.05792,
                                                 0.76512, -0.17833),
                                    `Lower 95% CI` = c(0.14373, 0.20471, 0.81617,
                                                       0.25162, NA),
                                    `Upper 95% CI` = c(1.28613, 2.55513, 1.37127,
                                                       2.32655, NA)),
                tolerance = 0.00001)

   # third test -----------------------------------------------------
   third_test <- mediator(data = test_dat,
                          out.model = glm(y ~ x + m + c + cens + x*m,
                                          data = test_dat),
                          med.model = glm(m ~ x + c + cens,
                                          data = test_dat, family = "binomial"),
                          treat = "x",
                          a = 1, a_star = 0, m = 0, boot_rep = 0)

   expect_equal(third_test, tibble(Effect = c("CDE", "NDE", "NIE",
                                               "Total Effect",
                                               "Proportion Mediated"),
                                   Estimate = c(0.26042, 0.09575, -0.00369,
                                                0.09206, -0.04012),
                                   `Lower 95% CI` = c(-0.12956, -0.25349, -0.03495,
                                                      -0.18416, NA),
                                   `Upper 95% CI` = c(0.65040, 0.44499, 0.02756,
                                                      0.36827, NA)),
                tolerance = 0.00001)

   # fourth test ----------------------------------------------------
   fourth_test <- mediator(data = test_dat,
                           out.model = glm(cens ~ x + m + c + y + x*m,
                                           data = test_dat, family = "binomial"),
                           med.model = glm(m ~ x + c + y,
                                           data = test_dat, family = "binomial"),
                           treat = "x",
                           a = 1, a_star = 0,
                           m = 0,
                           boot_rep = 0)

   expect_equal(fourth_test, tibble(Effect = c("CDE", "NDE", "NIE",
                                              "Total Effect",
                                              "Proportion Mediated"),
                                    Estimate = c(0.54184, 0.61768, 1.00201,
                                                 0.61892, -0.00325),
                                    `Lower 95% CI` = c(0.16286, 0.25922, 0.90767,
                                                       0.25987, NA),
                                    `Upper 95% CI` = c(1.80266, 1.47183, 1.10615,
                                                       1.47404, NA)),
                tolerance = 0.00001)
})
