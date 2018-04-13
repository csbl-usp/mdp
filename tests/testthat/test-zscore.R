context("correct output for z-score function")


data(example_data)
data(example_pheno)

control_samples <- example_pheno[example_pheno$Class == "baseline","Sample"] 



test_that("same z score values are returned when data is made negative",{

example_data_neg <- example_data*(-1)

zscore_pos <- compute_zscore(data = example_data,
                             control_samples = control_samples,
                             measure = "mean",
                             std = 2)

zscore_neg <- compute_zscore(data = example_data_neg,
                             control_samples = control_samples,
                             measure = "mean"
                             ,std = 2)

expect_equal(zscore_pos,zscore_neg)

zscore_pos <- compute_zscore(data = example_data,
                             control_samples = control_samples,
                             measure = "median",
                             std = 2)

zscore_neg <- compute_zscore(data = example_data_neg,
                             control_samples = control_samples,
                             measure = "median",
                             std = 2)

expect_equal(zscore_pos,zscore_neg)


})