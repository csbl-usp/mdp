#' Molecular Degree of Perturbation
#'
#' Based on the Molecular Distance to Health, this function calculates scores 
#' to each sample based on their perturbation from healthy
#'
#' @export
#' @param data \code{data frame} of gene expression data with the gene symbols 
#' in the row names
#' @param pdata \code{data frame} of phenodata with a column headed Class and the 
#' other headed Sample.
#' @param control_lab character \code{vector} specifying the control class
#' @param print set as default to TRUE for pdfs of the sample scores 
#' to be saved
#' @param directory (optional) character string of output directory
#' @param pathways (optional) \code{list} whose names are pathways and elements are 
#' genes in the pathway. see details section for more information
#' @param measure 'medan' as default, can change to 'median'.
#' \code{mean} will select for z-score and \code{median} will select for modified
#' z-score. (see details)
#' @param std \code{numeric} set as default to 2, this governs the thresholding of 
#' expression data. z-scored expression values with absolute value less than 'std' 
#' will be set to 0.
#' @param fraction_genes \code{numeric} fraction of genes that will contribute to the top perturbed genes. Set as default to 0.25
#' @param save_tables Set as default to TRUE. Tables of zscore and gene and 
#' sample scores will be saved.
#' @param file_name (optional) character string that will be added to the saved file names
#' @return A list: zscore, gene_scores, gene_freq, sample_scores, perturbed_genes
#' \itemize{
#' \item Z-score - z-score is calculated using the control samples to compute
#' the average and the standard deviation. The absolute value of this
#' matrix is taken and values less than the std are set to zero.
#' This z-score data frame is used to compute the gene and sample scores.
#' \item Gene scores - mean z-score value for each gene in each class
#' \item Gene frequency - frequency with which a gene has a non zero z-score 
#' value in each class
#' \item Sample scores - list containing sample scores for different genesets.
#' Sample scores are the sum of the z-scored gene values for each sample,
#' averaged for the number of genes in that geneset.
#' \item Perturbed genes - vector of the top fraction of genes that have
#' higher gene scores in the test classes compared to the control.
#' \item Pathways - if genesets are provided, they are ranked according to the 
#' signal-to-noise #' ratio of test sample scores versus control sample scores 
#' calculated using that geneset.
#'}
#' @examples
#' # basic run
#' mdp(example_data,example_pheno,'baseline')
#' # run with pathways
#' pathway_file <- system.file('extdata', 'ReactomePathways.gmt', 
#' package = 'mdp')
#' mypathway <- fgsea::gmtPathways(pathway_file) # load a gmt file
#' mdp(data=example_data,pdata=example_pheno,control_lab='baseline',
#' pathways=mypathway)
#' @importFrom utils write.table
#' @importFrom stats mad median sd
#' @section Loading pathways:
#' a \code{list} of pathways can be loaded from a .gmt file using the
#' \code{fgsea} function using \code{fgsea::gmtPathways('gmt.file.location')}
#' @section Selecting mean or median:
#' if \code{median} is selected, the z-score will be calculated using the
#' median, and the standard deviation will be estimated using the median
#'  absolute deviation, utilising the \code{mad} function.
mdp <- function(data, pdata, control_lab, directory = "", pathways, 
                print = TRUE, 
                measure = c("mean","median"), std = 2, fraction_genes = 0.25, 
                save_tables = TRUE, file_name = "") {
    
    
    # ---------- create directory ---------------------------------#######
    
    if (directory != "") {
        path = directory
        if (dir.exists(directory) == FALSE) {
            dir.create(directory)
        }
    } else {
        path = "."
    }
    
    
    # --------------------- ORGANISE DATA -------------------------#######
    
    if (missing(pdata)) {
        stop("Please include phenodata")
    }    
    if (missing(data)) {
        stop("Please include expression data")
    }
    if (missing(control_lab)) {
        stop("Please include control label")
    }   
        
    if (!("Sample" %in% names(pdata) && "Class" %in% names(pdata))) {
        stop("Please label phenodata columns as 'Sample' and 'Class'")
    } else if (sum(pdata$Sample %in% names(data)) == 0) {
        stop("Please provide phenodata sample names that match expression data columns")
    } else if (sum(control_lab %in% pdata$Class) == 0) {
        stop("Please provide a control label that matches a class in the phenodata")
    } else if (sum(pdata$Class %in% control_lab) < 2) {
        stop("Please provide at least two control samples (around 10 is an ideal minimum)")
    }
    
    if (!all(apply(data[, ], 2, is.numeric))) {
        stop("Please provide numeric values in expression data columns")
    }
    
    if (!is.numeric(fraction_genes)){
        stop("Select a numeric gene fraction of value less than 1")
    } else if (fraction_genes > 1) {
        stop("Select a gene fraction of less than 1")
    }
    
    
    if (!is.numeric(std)){
      stop("Select a numeric standard deviation threshold")
    }
    
    
    if (!missing(pathways)) {
        if (!is.list(pathways)) {
            stop("Please provide pathways in a list format (see help for more details")
        }
    }
    
    
    
    measure <- match.arg(measure)

    # Only keep the samples that have both pdata and data
    pdata <- pdata[as.character(pdata$Sample) %in% colnames(data), ]
    rownames(pdata) <- pdata$Sample
    # Expression data has only samples that are in pdata
    data <- data[, as.character(pdata$Sample)]
    
    control_samples <- as.character(pdata$Sample[pdata$Class == control_lab])
    test_samples <- as.character(pdata$Sample[pdata$Class != control_lab])
    
    message("Calculating Z score")
    zscore <- compute_zscore(data, control_samples, measure, std)
    
    message("Calculating gene scores")
    
    # find gene scores and gene frequency for each class
    
    gmdp_results <- compute_gene_score(zscore, 
                                        pdata, 
                                        control_lab, 
                                        "gene_score")
    
    gmdp_freq_results <- compute_gene_score(zscore, 
                                            pdata, 
                                            control_lab, 
                                            "gene_freq")
    
    # find perturbed genes
    perturbed_genes <- compute_perturbed_genes(gmdp_results, 
                                                control_lab, 
                                                fraction_genes)
    
    message("Calculating sample scores")
    # calculate sample scores
    sample_results <- compute_sample_scores(zscore, perturbed_genes, 
                                            control_samples, test_samples, 
                                            pathways, pdata)
    
    
    if (print == TRUE) {
        
        message("printing")
        
        sample_plot(sample_results[["allgenes"]], 
                    filename = paste0(file_name, "allgenes"), 
                    directory = path, 
                    print = TRUE, 
                    display = TRUE, 
                    title = "allgenes", 
                    control_lab = control_lab)
        
        sample_plot(sample_results[["perturbedgenes"]], 
                    filename = paste0(file_name, "perturbed"), 
                    directory = path, 
                    print = TRUE, 
                    display = FALSE, 
                    title = "perturbedgenes", 
                    control_lab = control_lab)
        
        if (!missing(pathways)) {
            pathway_results <- pathway_summary(sample_results, 
                                                path, file_name, 
                                                control_samples, 
                                                control_lab)
            
        }
        
    }
    
    if (save_tables == TRUE) {
        
        write.table(x = zscore, 
                    file = file.path(path, paste0(file_name, "zscore.tsv")), 
                    row.names = FALSE)
        write.table(x = gmdp_results, 
                    file = file.path(path, paste0(file_name, "gene_scores.tsv")), 
                    row.names = FALSE)
        write.table(x = sample_results, 
                    file = file.path(path, paste0(file_name, "sample_scores.tsv")), 
                    row.names = FALSE)
        
    }
    
    # ---------------- OUTPUT ---------------- ####
    if (missing(pathways)) {
        
        output <- list(zscore, 
                        gmdp_results, 
                        gmdp_freq_results, 
                        sample_results, 
                        perturbed_genes)
        
        names(output) <- c("zscore", 
                            "gene_scores", 
                            "gene_freq", 
                            "sample_scores", 
                            "perturbed_genes")
        
    } else {
        output <- list(zscore, 
                        gmdp_results, 
                        gmdp_freq_results, 
                        sample_results, 
                        perturbed_genes, 
                        pathway_results)
        
        names(output) <- c("zscore", 
                            "gene_scores", 
                            "gene_freq", 
                            "sample_scores", 
                            "perturbed_genes", 
                            "pathways")
        
    }
    
    
    
    return(output)
}


#' Computes the thresholded Z score
#' Plots the Z score using control samples to compute the average and standard 
#' deviation
#' @export
#' @param data Gene expression data with gene symbols in rows, sample names 
#' in columns
#' @param control_samples Character vector specifying the control sample names
#' @param measure Either 'mean' or 'median'. 'mean' uses mean and standard 
#' deviation. 'median' uses the median and the median absolute deviation to 
#' estimate the standard devation (modified z-score).
#' @param std Set as default to 2. This controls the standard deviation threshold 
#' for the Z-score calculation. #' Normalised expression values less than 'std' 
#' will be set to 0.
#' @examples
#' control_samples <- example_pheno$Sample[example_pheno$Class == 'baseline']
#' compute_zscore(example_data, control_samples,'median',2)
#' @return zscore data frame
compute_zscore <- function(data, control_samples, measure = c("mean","median"), std = 2) {
    
    measure <- match.arg(measure)
    
    
    if (measure == "mean") {
        stats = data.frame(mean = rowMeans(data[, control_samples]), 
                            std = apply(data[, control_samples], 1, sd))
        
    } else if (measure == "median") {
        # use the median absolute deviation to estimate the standard deviation
        
        stats = data.frame(mean = apply(data[, control_samples], 1, median),
                            std = apply(data[, control_samples], 1, mad))
        
    } else {
        stop("Please specify either 'mean' or 'median' as a measure input")
    }
    
    # compute Z score
    zscore <- (data - stats[, 1])/stats[, 2]
    rownames(zscore) <- rownames(data)
    zscore <- zscore[!is.infinite(rowSums(zscore)), ]
    zscore <- zscore[!is.na(rowSums(zscore)), ]
    zscore <- abs(zscore)
    zscore[zscore < std] <- 0
    zscore <- data.frame(Symbol = rownames(zscore), zscore)
    
    return(zscore)
}

#' Compute gene score
#' Computes gene scores for each gene within each class and perturbation freq
#' @param zscore zscore data frame
#' @param pdata phenotypic data with Class and Sample columns
#' @param control_lab character specifying control class
#' @param score_type set to 'gene_score' or 'gene_freq' to compute gene scores 
#' or frequencies
#' @return data frame of gene scores or gene frequencies
compute_gene_score <- function(zscore, pdata, control_lab, score_type = c("gene_score","gene_freq")) {
    
    score_type <- match.arg(score_type)
    
    all_groups <- unique(pdata$Class)  # find all groups
    # find all_groups, without the control label
    groups_no_control <- all_groups[-grep(control_lab, all_groups)]  
    
    gene_results <- zscore$Symbol  # table to save the gene scores
    
    for (group in all_groups) {
        if (score_type == "gene_score") {
            # find average expression of gene in each class
            gene_average <- rowMeans(zscore[, as.character(pdata$Sample[pdata$Class == group])])
            
        } else if (score_type == "gene_freq") {
            
            # find frequency that gene is perturbed in each class
            gene_average <- rowMeans(zscore[, as.character(pdata$Sample[pdata$Class == group])] > 0)
        }
        
        gene_results <- data.frame(gene_results, gene_average)
        
    }
    
    names(gene_results) <- c("Symbol", as.character(all_groups))
    
    return(gene_results)
}



#' Compute perturbed genes
#' Find the top fraction of genes that are more perturbed 
#' in test versus controls
#' @param gmdp_results results table of gene scores
#' @param control_lab label specificying control class
#' @param fraction_genes fraction of top perturbed genes that 
#' will make the set of perturbed genes
#' @return vector of perturbed genes
compute_perturbed_genes <- function(gmdp_results, 
                                    control_lab, fraction_genes) {
    
    control_idx <- grep(control_lab, names(gmdp_results))
    
    num_test_groups <- dim(gmdp_results)[2] - 2
    
    # score is proportional to: av. test group gene score - control
    score <- rowSums(gmdp_results[, -c(1, control_idx), drop = FALSE]) - 
                    (num_test_groups * gmdp_results[, control_idx])
    
    # order gene scores according to the difference in test vs control scores
    gmdp_results <- gmdp_results[order(-score), ]
    
    # perturbed genes are top fraction of ordered genes
    perturbed_genes <- gmdp_results$Symbol[1:round(fraction_genes * 
                                            dim(gmdp_results)[1])]
    
    
    return(perturbed_genes)
}





#' Compute sample scores for each pathway
#' @param zscore zscore data frame
#' @param perturbed_genes list of pertured genes
#' @param control_samples vector of control sample names
#' @param test_samples vector of test sample names
#' @param pathways list of pathways
#' @param pdata phenotypic data with Sample and Class columns
#' @return data frame of sample scores
compute_sample_scores <- function(zscore, perturbed_genes, control_samples, 
                                    test_samples, pathways, pdata) {
    
    genesets <- list()  # list of all gene sets
    genesets[[1]] <- zscore$Symbol
    genesets[[2]] <- perturbed_genes
    names(genesets)[1:2] <- c("allgenes", "perturbedgenes")
    
    if (!missing(pathways)) {
        genesets <- c(genesets, pathways)
    }
    
    # calculate sample scores for all genesets
    sample_results <- lapply(seq_along(genesets), function(idx) {
        sample_scores <- colMeans(zscore[zscore$Symbol %in% genesets[[idx]],
                                         2:ncol(zscore)])
        data.frame(Sample = names(sample_scores),
                   Score = sample_scores,
                   Class = pdata[names(sample_scores),
                                 "Class"])
    })
    
    names(sample_results) <- names(genesets)
    
    return(sample_results)
}


#' Plot sample scores
#' Plots the sample scores data.frame for a given geneset. 
#' Data frame must have Score, Sample and Class columns
#' @export
#' @param sample_data \code{data frame} of sample score information for a geneset. 
#' Must have columns 'Sample', 'Score' and 'Class'
#' @param filename (optional) character string that will be added to the saved pdf filename
#' @param directory (optional) character string of directory to save file
#' @param title (optional) character string of title name for graph
#' @param print (default TRUE) Save as a pdf file
#' @param display (default TRUE) Display plot
#' @param control_lab (optional) character string Specifying control_lab will set the control 
#' class as light blue as a default
#' @examples
#' sample_plot(sample_data = sample_data, control_lab = 'baseline')
#' @return generates a plot of the sample scores
sample_plot <- function(sample_data, filename = "", 
                        directory = "", title = "", print = TRUE, 
                        display = TRUE, control_lab) {
    

    if (!("Sample" %in% names(sample_data)) | 
        !("Class" %in% names(sample_data)) |
        !("Score" %in% names(sample_data))){
        stop("Sample data must be data frame with colnames 'Class' 'Sample' 'Score'")
    }
    
    
    if (missing(control_lab)) {
        stop("Please include control label")
    }   
    

    
    if (directory != "") {
        path = directory
        if (dir.exists(directory) == FALSE) {
            dir.create(directory)
        }
        
    } else {
        path = "."
    }
    
    if (length(unique(sample_data$Geneset)) > 1) {
        stop("Please subset sample scores for one geneset only")
    }
    
    
    
    # make color for each class, with control class as light blue
    groups <- unique(sample_data$Class)
    palette <- c("#86cce0", "#4da566", "#d67048", "#d67048", "#b59519", 
                "#b59519", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
   
    if (length(groups) > length(palette)) {
        palette <- rep(palette, ceiling(length(groups)/length(palette)))
    }
    groups_coloured <- palette[1:length(groups)]
    
    if (!missing(control_lab)) {
        if (!(control_lab %in% sample_data$Class)) {
            
            stop("Please provide control label that features in the sample data")
        } else {
            
            groups_reordered <- c(control_lab, 
                                    as.character(groups[-grep(control_lab, 
                                                        groups)]))
            
            names(groups_coloured) <- groups_reordered
        }
    }
    
    
    sample_data <- sample_data[order(sample_data$Score), ]
    sample_data$Sample <- factor(sample_data$Sample, 
                                levels = sample_data$Sample[order(sample_data$Score)])
    
    # Plot scores as bar graphs
    plot1 <- ggplot2::ggplot(data = sample_data,
                            ggplot2::aes_string(y = "Score", 
                                                x = "Sample", 
                                                fill = "Class")) +
                ggplot2::geom_bar(stat = "identity", 
                                    width = 0.8, 
                                    alpha = 0.7) + 
                ggplot2::labs(title = title, 
                                    x = "Samples", 
                                    y = "Sample score") + 
                ggplot2::theme(axis.line = ggplot2::element_line(size = 0.5, linetype = "solid"),
                                panel.grid.major = 
                                ggplot2::element_line(colour = "black", 
                                                    linetype = "blank"),
                                panel.grid.minor = 
                                        ggplot2::element_line(linetype = "blank"),
                                axis.title = ggplot2::element_text(size = 10),
                                axis.text.x = ggplot2::element_text(size = 8,
                                                                    angle = 60,
                                                                    hjust = 1),
                                plot.title = ggplot2::element_text(size = 12),
                                panel.background = 
                                        ggplot2::element_rect(fill = "grey100"),
                                legend.text = ggplot2::element_text(size = 8),
                                legend.key.size = ggplot2::unit(0.9, "line"),
                                legend.position = c(0.23, 0.85),
                                legend.direction = "vertical") + 
                ggplot2::theme(panel.background = 
                                ggplot2::element_rect(fill = NA,
                                                    linetype = "solid")) +
                ggplot2::scale_fill_manual(values = groups_coloured)
    
    
    
    
    # find means of sample scores
    class_means <- vector()
    for (j in unique(sample_data$Class)) {
        class_means = c(class_means, 
                        mean(sample_data[sample_data$Class == j, 
                                        "Score"]))
    }
    names(class_means) <- unique(sample_data$Class)
    class_means <- class_means[order(class_means)]
    
    
    ## Boxplot of scores
    plot2 <- ggplot2::ggplot(data = sample_data, 
                            ggplot2::aes_string(y = "Score", 
                                                x = "Class", 
                                                fill = "Class")) + 
            ggplot2::geom_boxplot(outlier.shape = NA) + 
            ggplot2::stat_summary(fun.y = mean, 
                                    geom = "point", 
                                    shape = 23, 
                                    size = 6) + 
            ggplot2::labs(title = title,
                            x = "Class", 
                            y = "Sample score") + 
            ggplot2::theme(legend.position = "null") + 
            ggplot2::scale_x_discrete(limits = names(class_means)) + 
            ggplot2::geom_jitter(shape = 16, 
                                position = ggplot2::position_jitter(0.2), 
                                size = 2, color = "grey10", alpha = 0.7) + 
        ggplot2::theme(axis.line = ggplot2::element_line(size = 0.5, 
                                                        linetype = "solid"), 
                        axis.title = 
                                    ggplot2::element_text(size = 10), 
                        panel.grid.major = 
                                    ggplot2::element_line(linetype = "blank"), 
                        panel.grid.minor = 
                                    ggplot2::element_line(linetype = "blank"), 
                        panel.background = 
                                    ggplot2::element_rect(fill = "white"), 
                        axis.text.x = 
                                    ggplot2::element_text(angle = 60, hjust = 1), 
                        legend.text = 
                                    ggplot2::element_text(size = 8), 
                        axis.text = ggplot2::element_text(size = 8)) + 
        ggplot2::scale_fill_manual(values = groups_coloured)
    
    
    if (print == TRUE) {
        sample_name <- paste(filename, "samples.pdf", sep = "")
        grDevices::pdf(file.path(path, sample_name))
        multiplot(plot1, plot2, cols = 2)
        grDevices::dev.off()
    }
    
    if (display == TRUE) {
        multiplot(plot1, plot2, cols = 2)
    }
    
}





# -------- Multiple plot function (R cookbook) ---- #### 
# ggplot objects can be passed in ..., or to plotlist 
# (as a list of ggplot2::ggplot
# objects) - cols: Number of columns in layout - layout: 
# A matrix specifying the layout. 
# If present, 'cols' is ignored.  If the layout is
# something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE), then 
# plot 1 will go in the upper left, 
# 2 will go in the upper right, and 3 will go all
# the way across the bottom.  
#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist = NULL, 
                    file, cols = 1, layout = NULL) {
    
    # Make a list from the ... arguments and plotlist
    
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel ncol: Number of columns of plots nrow: 
        # Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), 
                        ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots == 1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = 
                                    grid::grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        
        for (i in 1:numPlots) {
            
            # Get the i,j matrix positions of the regions 
            # that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], 
                  vp = grid::viewport(layout.pos.row = matchidx$row, 
                                    layout.pos.col = matchidx$col))
        }
    }
}


#' print pathways
#' generates a summary plot for pathways and sample score plot of best gene set
#' @param sample_results list of sample scores for each geneset
#' @param path directory to save images
#' @param file_name name of saved imaged
#' @param control_samples list of control sample names
#' @param control_lab label that specifies control class
#' @return data frame of signal to noise ratio of control vc test sample scores 
#' for each pathway
pathway_summary <- function(sample_results, path, file_name, 
                            control_samples, control_lab) {
    signal_noise <- vector()
    
    signal_noise <- vapply(seq_along(sample_results), function(idx) {
        sample_sub <- sample_results[[idx]]
        control_scores <- sample_sub[sample_sub$Sample %in% control_samples,
                                     "Score"]
        test_scores <- sample_sub[!(sample_sub$Sample %in% control_samples),
                                  "Score"]
        signal_noise[idx] <- (mean(test_scores) - mean(control_scores))/
            (sd(test_scores) + sd(control_scores))
    }, numeric(1)
    )
    
    pathway_scores <- data.frame(Geneset = names(sample_results), 
                                Sig2noise = signal_noise)
    
    pathway_scores <- pathway_scores[order(-pathway_scores$Sig2noise), ]
    
    # find best pathway
    top_pathway <- pathway_scores[1:3, "Geneset"]
    top_pathway <- top_pathway[top_pathway != "allgenes" & 
                                top_pathway != "perturbedgenes"]
    top_pathway <- top_pathway[1]
    
    
    sample_plot(sample_results[[top_pathway]], 
                filename = paste0(file_name, "bestPathway"), 
                directory = path, title = top_pathway, 
                print = TRUE, 
                display = TRUE, control_lab)
    
    return(pathway_scores)
    
}
