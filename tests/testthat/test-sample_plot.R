context("correct input for sample_plot")


data(example_data)
data(example_pheno)

mdp.results <- mdp(data=example_data, pdata=example_pheno, control_lab = "baseline")
sample_scores_list <- mdp.results$sample_scores
sample_data <- sample_scores_list[["perturbedgenes"]]


test_that("function does not run if sample data is missing correct classes",{
    
    names(sample_data) <- c("wrong","Score","Class")
    expect_error(sample_plot(sample_data,control_lab = "baseline"),
                 "Sample data")
    
    names(sample_data) <- c("Sample","wrong","Class")
    expect_error(sample_plot(sample_data,control_lab = "baseline"),
                 "Sample data")
    
    names(sample_data) <- c("Sample","Score","wrong")
    expect_error(sample_plot(sample_data,control_lab = "baseline")
                 ,"Sample data")
})
