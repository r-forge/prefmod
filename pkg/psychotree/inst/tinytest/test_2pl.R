### Check pltree with 2pl model on the SPISA dataset ###

## load data and fit model ##
data("SPISA", package = "psychotree")
fit_2PL <- pltree(spisa[,1:9] ~ age + gender, data = SPISA, type = "2PL", maxit = 5000)
# saveRDS(fit_2PL, "inst/tinytest/fit_2PL.RDS") #dont run - old test

# dont run, creation of smaller tests:
#test_2PL <- list()
#test_2PL[[1]] <- coef(fit_2PL)
#test_2PL[[2]] <- itempar(fit_2PL)
#test_2PL[[3]] <- threshpar(fit_2PL)
#test_2PL[[4]] <- guesspar(fit_2PL)
#test_2PL[[5]] <- upperpar(fit_2PL)
#saveRDS(test_2PL, "test_2PL.RDS")

## test methods
expect_equal(coef(fit_2PL), readRDS("test_2PL.RDS")[[1]],
             tolerance = 0.0001)
expect_equal(itempar(fit_2PL), readRDS("test_2PL.RDS")[[2]],
             tolerance = 0.0001)
expect_equal(threshpar(fit_2PL), readRDS("test_2PL.RDS")[[3]],
             tolerance = 0.0001)
expect_equal(guesspar(fit_2PL), readRDS("test_2PL.RDS")[[4]],
             tolerance = 0.0001)
expect_equal(upperpar(fit_2PL), readRDS("test_2PL.RDS")[[5]],
             tolerance = 0.0001)


### Check pltree with 2pl model and an impact factor on the SPISA dataset ###

## load data and fit model ##
data("SPISA", package = "psychotree")
fit_2PL_mg <- pltree(spisa[,1:9] ~ gender | age + gender, data = SPISA, type = "2PL", maxit = 5000)
# saveRDS(fit_2PL_mg, "inst/tinytest/fit_2PL_mg.RDS") #dont run

# dont run, creation of smaller tests:
# test_2PL_mg <- list()
# test_2PL_mg[[1]] <- coef(fit_2PL_mg)
# test_2PL_mg[[2]] <- itempar(fit_2PL_mg)
# test_2PL_mg[[3]] <- threshpar(fit_2PL_mg)
# test_2PL_mg[[4]] <- guesspar(fit_2PL_mg)
# test_2PL_mg[[5]] <- upperpar(fit_2PL_mg)
# saveRDS(test_2PL_mg, "test_2PL_mg.RDS")

## test methods
expect_equal(coef(fit_2PL_mg), readRDS("test_2PL_mg.RDS")[[1]],
             tolerance = 0.0001)
expect_equal(itempar(fit_2PL_mg), readRDS("test_2PL_mg.RDS")[[2]],
             tolerance = 0.0001)
expect_equal(threshpar(fit_2PL_mg), readRDS("test_2PL_mg.RDS")[[3]],
             tolerance = 0.0001)
expect_equal(guesspar(fit_2PL_mg), readRDS("test_2PL_mg.RDS")[[4]],
             tolerance = 0.0001)
expect_equal(upperpar(fit_2PL_mg), readRDS("test_2PL_mg.RDS")[[5]],
             tolerance = 0.0001)
