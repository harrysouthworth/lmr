library(testthat)
library(lmr)

test_that("Pena-Yohai initial estimates are the same every time", {
  liver <- texmex::liver

  for (run in 1:10){
    i <- sample(1:nrow(liver), size = nrow(liver), replace = TRUE)
    d <- liver[i, ]

    m1 <- py(log(ALT.M) ~ log(ALT.B) + log(ALP.B) + log(AST.B) + as.numeric(dose),
             data = d)
    m2 <- py(log(ALT.M) ~ log(ALT.B) + log(ALP.B) + log(AST.B) + as.numeric(dose),
             data = d)

    expect_identical(m1, m2, label = "py produces the same output every time")

    m1 <- py(log(ALT.M) ~ log(ALT.B) + log(ALP.B) + log(AST.B) + as.numeric(dose),
             data = d, resid_keep_method = "proportion",
             psc_keep = .75, resid_keep_prop = .75)

    expect_false(all(m1$coefficients == m2$coefficients),
                 label = "py produces different output with different settings")

    m2 <- py(log(ALT.M) ~ log(ALT.B) + log(ALP.B) + log(AST.B) + as.numeric(dose),
             data = d, resid_keep_method = "proportion",
             psc_keep = .75, resid_keep_prop = .75)

    expect_identical(m1, m2, label = "py produces the same output every time with user-supplied settings")
  }
})

test_that("robustbase and MASS functions give the same result", {
  liver <- texmex::liver

  set.seed(20210615)
  for (run in 1:10){
    i <- sample(1:nrow(liver), size = nrow(liver), replace = TRUE)
    d <- liver[i, ]

    mrb <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data = d,
               engine = "lmrob")
    mr <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data = d,
              engine = "rlm")
    expect_equal(coef(mr), coef(mrb), tol = .002,
                 label = "lmrob and rlm give same coefficients")
    expect_equal(mr$init$coefficients, mrb$init$coefficients, tol = 1e-7,
                 label = "lmrob and rlm use same initial coefficients")
    expect_equal(mr$init$scale, mrb$init$scale, tol = 1e-7,
                 label = "lmrob and rlm use same initial scale")
    expect_gt(cor(resid(mrb), resid(mr)), .999,
              label = "lmrob and rlm produce the same residuals")

    mrr <- robust::lmRob(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data = d,
                         control = robust::lmRob.control(efficiency = .85,
                                                         weight = c("bisquare", "bisquare"),
                                                         mxr = 200))
    expect_equal(coef(mrr), coef(mrb), tol = .002,
                 label = "lmRob gives same coefficients")
    expect_gt(cor(resid(mrb), resid(mrr)), .999,
              label = "lmRob and lmr produce same residuals")
  }
})


test_that("rfpe works correctly", {
  liver <- texmex::liver

  set.seed(20210615)
  for (i in 1:10){
    i <- sample(1:nrow(liver), size = nrow(liver), replace = TRUE)
    d <- liver[i, ]


    mrb0 <- lmr(log(ALT.M) ~ 1, data = d,
                engine = "lmrob")
    mrb1 <- lmr(log(ALT.M) ~ log(ALT.B), data = d,
               engine = "lmrob")
    mrb2 <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data = d,
                engine = "lmrob")

    s <- mrb2$scale

    expect_error(rfpe(mrb2),
                   label = "rfpe: failure if scale not provided")

    expect_lt(rfpe(mrb1, scale = s), rfpe(mrb0, scale = s),
              label = "rfpe: ok with lmrob, models 0 and 1")
    expect_lt(rfpe(mrb2, scale = s), rfpe(mrb1, scale = s),
              label = "rfpe: ok with lmrob, models 1 and 2")

    mrb0 <- lmr(log(ALT.M) ~ 1, data = d,
                engine = "rlm")
    mrb1 <- lmr(log(ALT.M) ~ log(ALT.B), data = d,
                engine = "rlm")
    mrb2 <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data = d,
                engine = "rlm")

    s <- mrb2$scale

    expect_lt(rfpe(mrb1, scale = s), rfpe(mrb0, scale = s),
              label = "rfpe: ok with rlm, models 0 and 1")
    expect_lt(rfpe(mrb2, scale = s), rfpe(mrb1, scale = s),
              label = "rfpe: ok with rlm, models 1 and 2")
  }
})



