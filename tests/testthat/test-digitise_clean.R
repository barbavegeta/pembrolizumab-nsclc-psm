test_that("cleaning enforces monotone survival, increasing time and (0, 1) start", {
  raw <- data.frame(
    time_months = c(2, 1, 3, 3, 5, 70),
    surv        = c(0.80, 0.95, 0.82, 0.70, 1.02, 0.30)
  )
  res <- clean_digitised_curve(raw, x_max = 66, label = "toy")
  cv <- res$curve

  expect_equal(cv$time_months[1], 0)
  expect_equal(cv$surv[1], 1)
  expect_true(all(diff(cv$time_months) >= 0))
  expect_true(all(diff(cv$surv) <= 0))
  expect_true(all(cv$surv >= 0 & cv$surv <= 1))
  expect_false(any(cv$time_months > 66))
  expect_setequal(
    res$log$rule,
    c("beyond_axis", "surv_clipped_to_unit_interval", "reordered_by_time",
      "running_minimum", "origin_added")
  )
})

test_that("points after the last numbers-at-risk time but within the axis are kept", {
  raw <- data.frame(time_months = c(0, 10, 61, 64), surv = c(1, 0.5, 0.3, 0.3))
  res <- clean_digitised_curve(raw, x_max = 66, label = "toy")
  expect_equal(nrow(res$curve), 4)
  expect_equal(nrow(res$log), 0)
})

test_that("readers reject malformed input", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("x,y", "0,1"), f)
  expect_error(read_digitised_curve(f), "missing column")
  writeLines(c("time_months,surv", "0,100", "1,90"), f)
  expect_error(read_digitised_curve(f), "percentage")
  writeLines(c("time_months,pembro,chemo", "0,154,151", "6,160,100"), f)
  expect_error(read_nar(f), "increase")
})
