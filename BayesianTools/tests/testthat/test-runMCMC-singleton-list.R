# Suggested regression test for BayesianTools/tests/testthat/.
for (sampler in c("DEzs", "DREAMzs")) {
  testthat::test_that(paste(sampler, "resumes singleton and longer sampler lists"), {
    set.seed(1)
    setup <- createBayesianSetup(
      likelihood = function(x) sum(dnorm(x, log = TRUE)),
      lower = rep(-10, 2),
      upper = rep(10, 2)
    )
    fit <- runMCMC(
      setup, sampler = sampler,
      settings = list(iterations = 30, burnin = 0, nrChains = 1,
                      startValue = 3, message = FALSE)
    )
    testthat::expect_s3_class(fit, "mcmcSampler")

    for (settings in list(NULL, list(iterations = 15, burnin = 0,
                                    message = FALSE))) {
      # A singleton list should resume just like the first entry of a longer list.
      set.seed(2)
      longer <- runMCMC(createMcmcSamplerList(list(fit, fit)), settings = settings)
      set.seed(2)
      singleton <- runMCMC(createMcmcSamplerList(list(fit)), settings = settings)

      testthat::expect_s3_class(singleton, "mcmcSamplerList")
      testthat::expect_length(singleton, 1L)
      testthat::expect_s3_class(longer, "mcmcSamplerList")
      testthat::expect_length(longer, 2L)
      testthat::expect_identical(singleton[[1]]$chain, longer[[1]]$chain)
      testthat::expect_identical(singleton[[1]]$X, longer[[1]]$X)
      testthat::expect_identical(singleton[[1]]$Z, longer[[1]]$Z)

      for (output in list(singleton, longer)) {
        for (i in seq_along(output)) {
          testthat::expect_s3_class(output[[i]], "mcmcSampler")
          for (j in seq_along(fit$chain)) {
            before <- as.matrix(fit$chain[[j]])
            after <- as.matrix(output[[i]]$chain[[j]])
            testthat::expect_gt(nrow(after), nrow(before))
            testthat::expect_identical(
              unname(after[seq_len(nrow(before)), , drop = FALSE]), unname(before)
            )
          }
        }
      }
    }
  })
}
