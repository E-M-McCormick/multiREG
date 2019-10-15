library(multiREG)

test_check("multiREG")

#-------------------------------------------------------# 
context("Check Output Structure")
#-------------------------------------------------------# 
output = multiREG(data = data("HRFsim"), header = FALSE, plot = FALSE)

expect_equal(length(output), 27)