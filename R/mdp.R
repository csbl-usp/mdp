#' Molecular Degree of Perturbation
#'
#' Based on the Molecular Distance to Health, this function allows the inspection of gene and sample heterogeneity in respect to a control class. When a gmt file is
#' submitted, the MDP is run on gene subsets. The MDP returns perturbation scores for each gene and each sample.
#'
#' @export
#' @param data A data.frame of gene expression data with the gene symbols in the row names
#' @param pdata A data.frame of phenodata with a column headed Class and the other headed Sample.
#' @param control_lab A character vector of the control class
#' @param print Set as default to TRUE if you wish graph pdfs of the geneMDP and sampleMDP values to be printed
#' @param directory The output directory (optional)
#' @param pathways A loaded gmt file in a list format (optional), use read_gmt("gmt.file.location") to load
#' @param measure Set as default to "mean", can be changed to "median". This measure is used in all Z-score calculations.
#' @param std Set as default to 2. This controls the standard deviation threshold for the Z-score calculation. Normalised expression values less than "std" will be set to 0.
#' @param save_tables Set as default to TRUE. Will save tables of zscore and gene and sample scores.
#' @return A list: [[1]] $zscore [[2]] $gene_scores  [[3]] $gene_freq [[4]] $sample_scores [[5]] perturbed_genes
#' @examples basic usage:
#' mdp(exp,pheno,"healthy_control")
#' mypathway <- read_gmt("gmt.file.location")
#' mdp(exp,pheno,"healthy_control",print=TRUE,directory="myexp",pathways=mypathway)
mdp <- function(data,pdata,control_lab,directory="",pathways,print=TRUE,measure="mean",std=2,fraction_genes = 0.25,save_tables = TRUE,file_name=""){


# ---------- create directory ---------------------------------#######
if (directory != ""){
  path = directory
  if (dir.exists(directory) == F){
    dir.create(directory)
  }
} else {
  path = "."
}

if (print == T){
progress <- progress_bar(4)
} else {
progress <- progress_bar(3)
}



# --------------------- ORGANISE DATA -------------------------------------#######
if (missing(pdata)){
  stop("Please include phenodata")
}  else if (!("Sample" %in% names(pdata) && "Class" %in% names(pdata))){
  stop("Please label phenodata columns as 'Sample' and 'Class'")
} else if (sum(pdata$Sample %in% names(data)) == 0){
  stop("Please provide phenodata sample names that match expression data columns")
} else if  (sum(control_lab %in% pdata$Class) == 0){
  stop("Please provide a control label that matches a class in the phenodata")
} else if (sum(pdata$Class %in% control_lab) < 2){
  stop("Please provide at least two control samples (around 10 is an ideal minimum)")
}
if (missing(data)){
  stop("Please include expression data")
} else if (!all(apply(data[,],2,is.numeric))){
  stop("Please provide numeric values in expression data columns")
}
if (fraction_genes > 1){
  stop("Select a gene fraction of less than 1")
}
if (!missing(pathways)){
  if (!is.list(pathways)){
    stop("Please provide pathways in a list format")
   }
}

pdata <- pdata[as.character(pdata$Sample) %in% colnames(data),]  # Only keep the samples that have both pdata and data
rownames(pdata) <- pdata$Sample
data <- data[,as.character(pdata$Sample)]  # Expression data has only samples that are in pdata

control_samples <- as.character(pdata$Sample[pdata$Class == control_lab])
test_samples <- as.character(pdata$Sample[pdata$Class != control_lab])


progress("Calculating Z score")
zscore <- compute_zscore(data,control_samples,measure,std) # make z score calculation

progress("Calculating gene scores")
gmdp_results <- compute_gene_score(zscore,pdata,control_lab,"gene_score")# find gene scores and gene frequency for each class
gmdp_freq_results <- compute_gene_score(zscore,pdata,control_lab,"gene_freq")

perturbed_genes <- compute_perturbed_genes(gmdp_results,control_lab,fraction_genes) # find perturbed genes

progress("Calculating sample scores")
sample_results <- compute_sample_scores(zscore,perturbed_genes,control_samples,test_samples,pathways,pdata) # calculate sample scores


if (print == T){
progress("printing")

  smdp_plot(sample_results[["allgenes"]],filename=paste0(file_name,"allgenes"),directory=path,title="allgenes",control_lab=control_lab)
  smdp_plot(sample_results[["perturbedgenes"]],filename=paste0(file_name,"perturbed"),directory=path,title="perturbedgenes",control_lab=control_lab)

    if (!missing(pathways)){


      pathway_results <- pathway_summary(sample_results,path,file_name,control_samples,control_lab)

    }

}

if (save_tables == T){

  write.table(x=zscore,file=file.path(path,paste0(file_name,"zscore.tsv")),row.names=F)
  write.table(x=gmdp_results,file=file.path(path,paste0(file_name,"gene_scores.tsv")),row.names=F)
  write.table(x=sample_results,file=file.path(path,paste0(file_name,"sample_scores.tsv")),row.names=F)


}

# ---------------- OUTPUT ---------------- ####
if (missing(pathways)){
  output <-list(zscore,gmdp_results,gmdp_freq_results,sample_results,perturbed_genes)
  names(output) <- c("zscore","gene_scores","gene_freq","sample_scores","perturbed_genes")
} else {
  output <-list(zscore,gmdp_results,gmdp_freq_results,sample_results,perturbed_genes,pathway_results)
  names(output) <- c("zscore","gene_scores","gene_freq","sample_scores","perturbed_genes","pathways")

}



return(output)
}


#' Computes the thresholded Z score
#' Plots the Z score using control samples to compute the average and standard deviation
#' @export
#' @param data gene expression data with gene symbools in rows
#' @param control_samples a character vector of control sample names
#' @param measure either "mean" or "median". mean uses mean and standard deviation. "median" uses modified z score.
compute_zscore <- function(data,control_samples,measure,std){

  if (measure == "mean"){
    stats = data.frame("mean" = rowMeans(data[,control_samples]), "std" = apply(data[,control_samples],1,sd))
  } else if (measure == "median"){
    stats = data.frame("mean" = apply(data[,control_samples],1,median), "std" = apply(data[,control_samples],1,mad) * (1 / 0.6745) )
  } else {
    stop("Please specify either 'mean' or 'median' as a measure input")
  }

  # compute Z score
  zscore <-(data - stats[, 1]) / stats[, 2]
  rownames(zscore) <- rownames(data)
  zscore <- zscore[!is.infinite(rowSums(zscore)),]
  zscore <- zscore[!is.na(rowSums(zscore)),]
  zscore <- abs(zscore)
  zscore[zscore < std] <- 0
  zscore <- data.frame("Symbol" = rownames(zscore), zscore)

  return(zscore)
}

#' Compute gene score
#' Computes gene scores for each gene within each class and perturbation freq
compute_gene_score <- function(zscore,pdata,control_lab,score_type){

all_groups <- unique(pdata$Class) # find all groups
groups_no_control <- all_groups[-grep(control_lab,all_groups)] # all_groups, without the control label

gene_results <- zscore$Symbol # table to save the gene scores

for (group in all_groups){

  if (score_type == "gene_score"){
  gene_average <- rowMeans(zscore[,as.character(pdata$Sample[pdata$Class == group])]) # find average expression of gene in each group
  } else {
  gene_average <- rowMeans(zscore[,as.character(pdata$Sample[pdata$Class == group])] > 0)
  }

  gene_results <- data.frame(gene_results, gene_average)

}

names(gene_results) <- c("Symbol",as.character(all_groups))

return(gene_results)
}



#' Compute perturbed genes
#' Find the top fraction of genes that are more perturbed in test versus controls
compute_perturbed_genes <- function(gmdp_results,control_lab,fraction_genes){

  control_idx <- grep(control_lab,names(gmdp_results))
  num_test_groups <- dim(gmdp_results)[2] - 2
  score <- rowSums(gmdp_results[,-c(1,control_idx), drop=FALSE]) - (num_test_groups * gmdp_results[,control_idx])
  gmdp_results <- gmdp_results[order(-score),]
  perturbed_genes <- gmdp_results$Symbol[1:round(fraction_genes*dim(gmdp_results)[1])]


return(perturbed_genes)
}





#' Compute sample scores for each pathway
compute_sample_scores <- function(zscore,perturbed_genes,control_samples,test_samples,pathways,pdata){

genesets <- list() # list of all gene sets
genesets[[1]] <- zscore$Symbol
genesets[[2]] <- perturbed_genes
names(genesets)[1:2] <- c("allgenes","perturbedgenes")
if (!missing(pathways)){
  genesets <- c(genesets,pathways)
}

sample_results <- list() # calculate sample scores for all genesets
for (idx in 1:length(genesets)){
  sample_scores <- colMeans(zscore[zscore$Symbol %in% genesets[[idx]],2:ncol(zscore)]) # average gene expression for each sample
  sample_results[[idx]] <-  data.frame("Sample" = names(sample_scores),
                                       "Score" = sample_scores,
                                       "Class" = pdata[names(sample_scores),"Class"])
}
names(sample_results) <- names(genesets)

return(sample_results)
}


#' Plot sMDP scores
#'
#' Plots the sMDP scores for a given geneset. Takes a data.frame that contains the sMDP scores for one geneset, along with the Sample and Class information. Data.frame must have a sMDP, Sample and Class columns.
#' @export
#' @param sMDP.data a data.frame containing sMDP information for a geneset, with columns "Sample", "sMDP" and "Class"
#' @param filename filename
#' @param directory directory to save file
#' @param title title name for graph
smdp_plot <- function(sample_data,filename=file_name,directory="",title="",print=T,control_lab){


  if (directory != ""){
    path = directory
    if (dir.exists(directory) == F){
      dir.create(directory)
    }
  } else {
    path = "."
  }
  if (missing(filename)){
    print = F
  }
  if (length(unique(sample_data$Geneset)) > 1){
    stop("Please subset sample scores for one geneset only")
  }



    # make color for each class, with control class as light blue
  groups <- unique(sample_data$Class)
  palette <- c("#86cce0","#4da566","#d67048","#d67048","#b59519","#b59519","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A")
  if (length(groups) > length(palette)){
        palette <- rep(palette, ceiling(length(groups)/length(palette)))
  }
  groups_coloured <- palette[1:length(groups)]

  if (!missing(control_lab)){
    if (!(control_lab %in% sample_data$Class)){
        stop("Please provide control label that features in the sample data")
    } else {

   groups_reordered <- c(control_lab,as.character(groups[-grep(control_lab,groups)]))
   names(groups_coloured) <- groups_reordered
    }
  }


  sample_plot <-sample_data[order(sample_data$Score),]
  sample_plot$Sample <- factor(sample_plot$Sample, levels = sample_plot$Sample[order(sample_plot$Score)])

  #Plot sMDP as bar graphs
  plot1 <- ggplot2::ggplot(data =  sample_plot, ggplot2::aes(y = Score, x = Sample, fill = Class)) +
    ggplot2::geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
    ggplot2::labs(title = title, x = "Samples", y = "Sample score") +
    ggplot2::theme(axis.line = ggplot2::element_line(size = 0.5,
                   linetype = "solid"), panel.grid.major = ggplot2::element_line(colour = "black",
                   linetype = "blank"), panel.grid.minor = ggplot2::element_line(linetype = "blank"),
                   axis.title = ggplot2::element_text(size = 12),
                   axis.text.x = ggplot2::element_text(size = 8, angle = 60, hjust=1),
                   plot.title = ggplot2::element_text(size = 12),
                  panel.background = ggplot2::element_rect(fill = "grey100"),
                   legend.text=ggplot2::element_text(size=8),
                   legend.key.size =  ggplot2::unit(0.9,"line"),
                   legend.position = c(0.23, 0.85),
                   legend.direction = "vertical",) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = NA, linetype = "solid")) +
    ggplot2::scale_fill_manual(values = groups_coloured)




  # find means of sample scores
 class_means <- vector()
 for (j in unique(sample_plot$Class)){
   class_means = c(class_means,mean(sample_plot[sample_plot$Class == j,"Score"]))
  }
  names( class_means) <- unique(sample_plot$Class)
  class_means <-  class_means[order( class_means)]


  ##Boxplot of sMDP score for each class
  plot2 <- ggplot2::ggplot(data = sample_plot, ggplot2::aes(y = Score, x = Class, fill = Class)) +
    ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::stat_summary(fun.y = mean, geom = "point", shape = 23, size = 6) +
    ggplot2::labs(title = title, x = "Class", y = "Sample score") +
    ggplot2::theme(legend.position = "null") +
    ggplot2::scale_x_discrete(limits=names(class_means)) +
    ggplot2::geom_jitter(shape = 16, position = ggplot2::position_jitter(0.2), size=2, color = "grey10", alpha=0.7) +
    ggplot2::theme(axis.line = ggplot2::element_line(size = 0.5, linetype = "solid"),
                   panel.grid.major = ggplot2::element_line(linetype = "blank"),
                   panel.grid.minor = ggplot2::element_line(linetype = "blank"),
                   panel.background = ggplot2::element_rect(fill = "white"),
                   axis.text.x = ggplot2::element_text(angle=60, hjust=1),
                   legend.text=ggplot2::element_text(size=8),
                   axis.text = ggplot2::element_text(size = 8)) +
    ggplot2::scale_fill_manual(values = groups_coloured)


  if (print == T){
  smdp.name <- paste(filename,"samples.pdf",sep="")
  pdf(file.path(path,smdp.name))
  multiplot(plot1,plot2,cols=2)
  dev.off()
  }

  return(plot2)

}



progress_bar <- function(n_steps, nth_step=0) {
  return(function(msg){
    nth_step <<- nth_step + 1
    cat(nth_step, "/", n_steps, ": ", msg, "\n", sep="")
  })
}



# -------- Multiple plot function (R cookbook) ---- ####
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot2::ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {


  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
  }
}


#' print pathways
#' summary plot for pathways and sample score plot of best gene set
pathway_summary <- function(sample_results,path,file_name,control_samples,control_lab){


  signal_noise <- vector()
  for (idx in 1:length(sample_results)){

  sample_sub <- sample_results[[idx]]
  control_scores <- sample_sub[sample_sub$Sample %in% control_samples,"Score"]
  test_scores <- sample_sub[!(sample_sub$Sample %in% control_samples),"Score"]
  signal_noise[idx] <- ( mean(test_scores) - mean(control_scores) ) /  ( sd(test_scores) + sd(control_scores) ) # signal to noise score
  }

  pathway_scores <- data.frame("Geneset" = names(sample_results), "Sig2noise" = signal_noise)
  pathway_scores <- pathway_scores[order(-pathway_scores$Sig2noise),]
  top_pathway <- pathway_scores[1,"Geneset"]


  if ( top_pathway != "allgenes" &  top_pathway != "perturbedgenes"){
    smdp_plot(sample_results[[top_pathway]],filename=paste0(file_name,"_",top_pathway),directory=path,title=top_pathway,print=T,control_lab)
  }

  return(pathway_scores)

}


