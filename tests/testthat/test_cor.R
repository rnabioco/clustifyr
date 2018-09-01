context("")

test_that("pearson works", {
  res <- lapply(
    c(
      "pearson",
      "cosine",
      "spearman"
    ),
    function(x) corr_coef(1:10, 1:10, method = x)
  )
  expect_equal(res[[1]], 1)
  expect_equal(res[[2]], 1)
  expect_equal(res[[3]], 1)
})
