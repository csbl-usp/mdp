context("correct input files for mdp function")

data(example_data)
data(example_pheno)


test_that("errors are returned when there are missing inputs",{
    
    expect_error(mdp(),
                 "Please*")
    expect_error(mdp(data=example_data),
                 "Please*")
    expect_error(mdp(pdata=example_pheno),
                 "Please*")
    expect_error(mdp(data=example_data,pdata=example_pheno),
                 "Please*")  
    
})



test_that("errors are returned when pheno data is not in correct format",{
    
    pheno_wrong <- example_pheno
    names(pheno_wrong) <- c("Label1","Label2")
    
    expect_error(mdp(data=example_data,pdata=pheno_wrong,control_lab="baseline"),
                 "*phenodata*")  
    
    expect_error(mdp(data=example_data,pdata=example_pheno,control_lab="healthy"),
                 "*control*")  
    expect_error(mdp(data=example_data,pdata=example_pheno[20:40,],control_lab="baseline"),
                 "*two control*")  
    
})


test_that("errors are returned when data and input parameters are not numeric",{
    
    data_wrong <- example_data
    data_wrong[,1] <- as.character(example_data[,1])
    
    expect_error(mdp(data=data_wrong,pdata=example_pheno,control_lab="baseline"),
                 "*numeric*")  
    
    expect_error(mdp(data=example_data,pdata=example_pheno,control_lab="baseline",fraction_genes = "n"),
                 "*numeric*")  
    
    expect_error(mdp(data=example_data,pdata=example_pheno,control_lab="baseline",std = "n"),
                 "*numeric*")  
    
    
})

test_that("errors are returned when pathways are not a list",{
    
    expect_error(mdp(data=example_data,pdata=example_pheno,control_lab="baseline",pathways="pathways"),
                 "*list*")  

})


