library(testthat)
library(PathwayVote)  # 替换为你包的名称

test_that("eQTM object can be created with empty metadata", {
  test_data <- data.frame(
    cpg = c("cg000001", "cg000002"),
    statistics = c(0.5, -0.7),
    p_value = c(0.01, 0.05),
    distance = c(500, 1200),
    entrez = c("1", "2"),
    stringsAsFactors = FALSE
  )

  # metadata 为空
  test_obj <- create_eQTM(data = test_data, metadata = list())

  expect_s4_class(test_obj, "eQTM")
  expect_true(is.data.frame(getData(test_obj)))
  expect_true(is.list(getMetadata(test_obj)))
  expect_equal(length(getMetadata(test_obj)), 0)
})

test_that("eQTM object fails when required columns are missing", {
  bad_data <- data.frame(
    statistics = c(0.1),
    p_value = c(0.01),
    distance = c(100),
    entrez = c("1"),
    stringsAsFactors = FALSE
  )

  expect_error(create_eQTM(bad_data), "Missing required columns")
})

test_that("eQTM object requires at least entrez or ensembl", {
  bad_data <- data.frame(
    cpg = "cg000003",
    statistics = 0.2,
    p_value = 0.05,
    distance = 1000,
    stringsAsFactors = FALSE
  )

  expect_error(create_eQTM(bad_data), "At least one of entrez or ensembl must be provided")
})
