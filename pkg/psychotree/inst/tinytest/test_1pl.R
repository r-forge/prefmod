### Check pltree with 1pl model on the SPISA dataset ###

## load data and fit model ##
data("SPISA", package = "psychotree")
fit_1PL <- pltree(spisa ~ age, data = SPISA, type = "1PL", maxit = 5000)
# saveRDS(fit_1PL, "inst/tinytest/fit_1PL.RDS") #dont run, old test
# dont run, creation of smaller test:
# test_1PL <- list()
# test_1PL[[1]] <- coef(fit_1PL)
# test_1PL[[2]] <- itempar(fit_1PL)
# test_1PL[[3]] <- threshpar(fit_1PL)
# test_1PL[[4]] <- guesspar(fit_1PL)
# test_1PL[[5]] <- upperpar(fit_1PL)
# saveRDS(test_1PL, "inst/tinytest/test_1PL.RDS")

## test methods
expect_equal(coef(fit_1PL), readRDS("test_1PL.RDS")[[1]],
             tolerance = 0.0001)
expect_equal(itempar(fit_1PL), readRDS("test_1PL.RDS")[[2]],
             tolerance = 0.0001)
expect_equal(threshpar(fit_1PL), readRDS("test_1PL.RDS")[[3]],
             tolerance = 0.0001)
expect_equal(guesspar(fit_1PL), readRDS("test_1PL.RDS")[[4]],
             tolerance = 0.0001)
expect_equal(upperpar(fit_1PL), readRDS("test_1PL.RDS")[[5]],
             tolerance = 0.0001)

### Check pltree with 1pl model and an impact factor on the SPISA dataset ###

## load data and fit model ##
data("SPISA", package = "psychotree")
fit_1PL_mg <- pltree(spisa ~ gender | gender, data = SPISA, type = "1PL", maxit = 2000)
# saveRDS(fit_1PL_mg, "inst/tinytest/fit_1PL_mg.RDS") #dont run

# dont run, creation of smaller tests:
# test_1PL_mg <- list()
# test_1PL_mg[[1]] <- coef(fit_1PL_mg)
# test_1PL_mg[[2]] <- itempar(fit_1PL_mg)
# test_1PL_mg[[3]] <- threshpar(fit_1PL_mg)
# test_1PL_mg[[4]] <- guesspar(fit_1PL_mg)
# test_1PL_mg[[5]] <- upperpar(fit_1PL_mg)
# saveRDS(test_1PL_mg, "inst/tinytest/test_1PL_mg.RDS")

## test methods
expect_equal(coef(fit_1PL_mg), readRDS("test_1PL_mg.RDS")[[1]])
expect_equal(itempar(fit_1PL_mg), readRDS("test_1PL_mg.RDS")[[2]])
expect_equal(threshpar(fit_1PL_mg), readRDS("test_1PL_mg.RDS")[[3]])
expect_equal(guesspar(fit_1PL_mg), readRDS("test_1PL_mg.RDS")[[4]])
expect_equal(upperpar(fit_1PL_mg), readRDS("test_1PL_mg.RDS")[[5]])

