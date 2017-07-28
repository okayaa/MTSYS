context("generates_transformation_functions_Tb works correctly")

test_that("generates_transformation_functions_Tb", {
  
  tmp <- generates_transformation_functions_Tb(stackloss[-c(2, 12, 19), ])
  center_subtraction_function <- tmp[[1]]
  subtracts_M_0 <- tmp[[2]]
  adds_M_0 <- tmp[[3]]
  
  transformed_data <- center_subtraction_function(stackloss[-c(2, 12, 19), -4])
  derived_center <- attr(transformed_data, "scaled:center")
  calculated_M <- subtracts_M_0(stackloss[-c(2, 12, 19), 4])
  calculated_y <- adds_M_0(calculated_M)
  
  correct_center <- c(50, 18, 72)
  correct_M <- c(34.531, 29.531, 20.531, 10.531, 10.531, 11.531, 12.531, 7.531, 
                 6.531, 6.531, 3.531, 4.531, 0.531, -0.469, 0.531, 0.531, 
                 7.531, 7.531)
  
  expect_equal(derived_center, correct_center)  
  expect_equal(round(calculated_M, 3), correct_M)
  expect_equal(calculated_y, stackloss[-c(2, 12, 19), 4])
  
})

