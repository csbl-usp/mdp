context("correct output for compute_gene function")


data(example_data)
data(example_pheno)

control_lab <- "baseline"
control_samples <- example_pheno[example_pheno$Class == control_lab,"Sample"] 

class_lab <- "Asymptomatic"
class_samples <- example_pheno[example_pheno$Class == class_lab,"Sample"] 

class2_lab <- "Symptomatic"
class2_samples <- example_pheno[example_pheno$Class == class2_lab,"Sample"] 

zscore <- compute_zscore(example_data, control_samples)
zscore_binary <- zscore[,2:ncol(zscore)]
zscore_binary[zscore_binary >= 2] <- 1
zscore_binary <- cbind(zscore[,1],zscore_binary)

gmdp_freq_results <- compute_gene_score(zscore, example_pheno, control_lab,"gene_freq")
gmdp_results <- compute_gene_score(zscore, example_pheno, control_lab)


test_that("check that gene score is the mean z-score value for that class",{
    
 expect_equal(as.numeric(gmdp_results[,"baseline"]),
              as.numeric(rowMeans(zscore[,as.character(control_samples)])))
    
 expect_equal(gmdp_results[,class_lab],
              as.numeric(rowMeans(zscore[,as.character(class_samples)])))
 
 expect_equal(gmdp_results[,class2_lab],
              as.numeric(rowMeans(zscore[,as.character(class2_samples)])))
 
})

test_that("check that gene frequency score is the frequency that gene is perturbed in each class",{
    
    expect_equal(as.numeric(gmdp_freq_results[,"baseline"]),
                 as.numeric(rowMeans(zscore_binary[,as.character(control_samples)])))
    
    expect_equal(gmdp_freq_results[,class_lab],
                 as.numeric(rowMeans(zscore_binary[,as.character(class_samples)])))
    
    expect_equal(gmdp_freq_results[,class2_lab],
                 as.numeric(rowMeans(zscore_binary[,as.character(class2_samples)])))
    
    
})

