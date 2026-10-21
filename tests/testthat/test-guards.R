test_that("residual df is n minus the rank of the design matrix", {
  d <- data.frame(a = factor(c("x", "x", "y", "y")), b = factor(c("p", "q", "p", "q")))
  expect_equal(residual_df(~ a + b, d), 1L)
  expect_equal(residual_df(~ a * b, d), 0L)
  # A confounded design: two factors that coincide add no rank, so df is not 4 - 3.
  d$c <- d$a
  expect_equal(residual_df(~ a + c, d), 2L)
})

test_that("assert_df blocks sample-level tests at 0 df and flags 1 df as exploratory", {
  d <- data.frame(g = factor(c("a", "b")))
  expect_error(assert_df(~g, d, "one vs one"), "0 residual degrees of freedom")
  d4 <- data.frame(a = factor(c("x", "x", "y", "y")), b = factor(c("p", "q", "p", "q")))
  expect_warning(res <- assert_df(~ a + b, d4, "additive"), "exploratory")
  expect_true(res$exploratory)
  d6 <- data.frame(g = factor(rep(c("a", "b"), 3)))
  expect_silent(res <- assert_df(~g, d6, "six samples"))
  expect_equal(res$df, 4L)
})
