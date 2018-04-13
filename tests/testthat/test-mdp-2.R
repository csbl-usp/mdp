context("correct output for mdp function")



data(example_data)
data(example_pheno)

file_address <- system.file("extdata", "ReactomePathways.gmt", package = "mdp")
pathways <- fgsea::gmtPathways(file_address)

test_that("list of correct size is returned",{
    
    output <- mdp(data=example_data,pdata=example_pheno,control_lab="baseline")
    expect_that(length(output),equals(5))  
    
    output <- mdp(data=example_data,pdata=example_pheno,control_lab="baseline",pathways=pathways)
    expect_that(length(output),equals(6))  

})    
    
   