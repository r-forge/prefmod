### Check pltree with rasch model on the SPISA dataset ###

## load data and fit model ##
data("SPISA", package = "psychotree")
fit_rasch <- pltree(spisa ~ age + gender + semester + elite + spon, data = SPISA, type = "Rasch")
# saveRDS(fit_rasch, "inst/tinytest/fit_rasch.RDS") #dont run
# dont run, creation of smaller test:
# test_rasch <- list()
# test_rasch[[1]] <- coef(fit_rasch)
# test_rasch[[2]] <- itempar(fit_rasch)
# test_rasch[[3]] <- threshpar(fit_rasch)
# test_rasch[[4]] <- guesspar(fit_rasch)
# test_rasch[[5]] <- upperpar(fit_rasch)
# saveRDS(test_rasch, "inst/tinytest/test_rasch.RDS")

## compare to previously fitted model ##
# Removed - takes too much spache
# expect_equal_to_reference(fit_rasch, "fit_rasch.RDS")

## check that a warning is emitted if I try to supply an impact factor, but I am doing CML
expect_warning(pltree(spisa ~ age + gender + semester + elite + spon | gender, data = SPISA, type = "Rasch"))

## test methods
expect_equal(coef(fit_rasch), readRDS("test_rasch.RDS")[[1]],
             tolerance = 0.0001)
expect_equal(itempar(fit_rasch), readRDS("test_rasch.RDS")[[2]],
             tolerance = 0.0001)
expect_equal(threshpar(fit_rasch), readRDS("test_rasch.RDS")[[3]],
             tolerance = 0.0001)
expect_equal(guesspar(fit_rasch), readRDS("test_rasch.RDS")[[4]],
             tolerance = 0.0001)
expect_equal(upperpar(fit_rasch), readRDS("test_rasch.RDS")[[5]],
             tolerance = 0.0001)

