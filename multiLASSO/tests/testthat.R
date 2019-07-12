library(testthat)
library(multiLASSO)

test_check("multiLASSO")

#-------------------------------------------------------# 
context("Check Output Structure")
#-------------------------------------------------------# 
output = multiLASSO(data = data("HRFsim"), header = FALSE, plot = FALSE)

expect_equal(length(output), 50)