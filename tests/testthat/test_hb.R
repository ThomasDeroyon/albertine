
context('hb')

test_that("hb works", {

 phat <-  suppressWarnings(
   glm(data = ict_sample,
              resp ~ as.factor(division) + turnover + employees,
              family = 'binomial')
   )

 expect_error(hb(phat = phat$fitted.values, parameters = 0.9),
              regexp = NA)

}
)

