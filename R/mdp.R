#' Molecular Degree of Perturbation
#'
#' Based on the Molecular Distance to Health, this function allows the inspection of gene and sample heterogeneity in respect to a control class. When a gmt file is
#' submitted, the MDP is run on gene subsets. The MDP returns perturbation scores for each gene and each sample.
#'
#' @export
#' @param data A data.frame of gene expression data with the gene symbols in the first column
#' @param pdata A data.frame of phenodata with a column headed Class and the other headed Sample.
#' @param control_lab A character vector of the control class
#' @param print Set as default to TRUE if you wish graph pdfs of the geneMDP and sampleMDP values to be printed
#' @param directory The output directory (optional)
#' @param pathways The location of a gmt file (optional)
#' @param measure Set as default to "median", can be changed to "mean". This measure is used in all Z-score calculations.
#' @param std Set as default to 2. This controls the standard deviation threshold for the Z-score calculation. Normalised expression values less than "std" will be set to 0.
#' @return A list where [[1]] $Zscore contains a table of Z scores, [[2]] $gMDP contains gMDP scores and [[3]] $sMDP are the sampleMDP scores
#' @examples mdp(exp,pheno,"healthy_control")
#' @examples mdp(exp,pheno,"healthy_control",print=TRUE,directory="myexp",pathways="mypathways.gmt")
mdp <- function(data,pdata,control_lab,directory="",pathways,print=TRUE,measure="median",std=2){


# --------------- FUNCTIONS - CALCULATE Z SCORE AND CALCULATE CONTROL STATS --- ###

  Z <- function(exp,health,std){
    # COMPUTE THE Z SCORE NORMALISED BY CONTROLS
    # exp <- column vector of expression values per sample, rows = genes
    # health <- mean of controls in first column, SD of controls in second column
    # n <- if Z is less than n standard deviations, set to NA
    z   <- (exp-health[,1])/health[,2]
    z[abs(z)<std] <- 0

    return(z)
  }

  GetMeanSD <- function(xx,measure) {
    controlData <- vector()

    xx <- as.numeric(xx)
    ncontrol  <- length(xx)

    if (measure == "mean"){
      meanX <- mean(xx,na.rm=TRUE)
    }
    else if (measure == "median"){
	    meanX  <- median(xx,na.rm=T)
    }


   sdX   <- sd(xx,na.rm=TRUE)

  controlData <- as.numeric(c(meanX,sdX))
  names(controlData) <- c("Mean","sd")
  return(controlData)
  }



new_progress <- function(n_steps, nth_step=0) {
	return(function(msg){
		nth_step <<- nth_step + 1
		cat(nth_step, "/", n_steps, ": ", msg, "\n", sep="")
	})
}



# ---------- create directory ---------------------------------#######


if (directory != ""){
  path = directory
  if (dir.exists(directory) == F){
    dir.create(directory)
  }
} else {
  path = "."
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
} else if (!all(apply(data,2,is.numeric))){
  stop("Please provide numeric values in expression data columns")
}

pdata <- pdata[as.character(pdata$Sample) %in% colnames(data),]  # Only keep the samples that have both pdata and data
rownames(pdata) <- pdata$Sample
data <- data[,as.character(pdata$Sample)]  # Expression data has c(1) and samples

progress <- new_progress(5)


# --------------- ORGANISE GROUPS AND CONTROLS FROM PHENO DATA -------------- ####
nGroups <- length(unique(pdata$Class))
idx <- list() # list of which indexes in pdata match each class type
for (i in 1:nGroups){
  idx[[i]] <- which(pdata$Class == unique(pdata$Class)[i])
  if (unique(pdata$Class)[i] == control_lab){
  control_idx <- i
  }
}
names(idx) <- unique(pdata$Class)


# ----------- FIND MEAN, SD FOR HEALTHY CONTROLS
progress("Computing Z-score for contrast samples")


dataHealth <- t(apply(data[,idx[[control_idx]]],1,function(x) {
  GetMeanSD(as.numeric(x),measure)}))



# --------------- COMPUTE Z SCORE ------------------------- ####

progress("Computing Z-score for reference samples")

N <- 1 # min number of standard deviations that are accepted

# Z score for all samples apart from controls
Zscore <- matrix(nrow=nrow(data), ncol=ncol(data))
rownames(Zscore) <- rownames(data)

colnames(Zscore) <- names(data)
Zscore[,-idx[[control_idx]]] <- apply(data[,-idx[[control_idx]]] ,2,function(x) Z(x,dataHealth,std))

# Z score for controls using leave-one-out
Zscore[,idx[[control_idx]]] <- sapply(1:length(idx[[control_idx]]), function(x) {
  # calculate mean using leave-one out
    dataHealthSubset.mean <- t(apply(data[, idx[[control_idx]][-x]],1,function(j) {
      GetMeanSD(as.numeric(j),measure)}))
  # calculate standard deviation using all of the control samples
    dataHealthSubset.sd <- t(apply(data[, idx[[control_idx]]],1,function(j) {
      GetMeanSD(as.numeric(j),measure)}))

    dataHealthSubset <- cbind(dataHealthSubset.mean[,1], dataHealthSubset.sd[,2])

  return(Z(data[,idx[[control_idx]][x]],dataHealthSubset,std))
})


Zscore <- Zscore[!is.infinite(rowSums(Zscore)),]
Zscore <- Zscore[!is.na(rowSums(Zscore)),]
Zscore.nonabsolute <- data.frame("Symbol" = rownames(Zscore), Zscore)
Zscore <- abs(Zscore)
Zscore.annotated <- data.frame("Symbol" = rownames(Zscore), Zscore)

# print Z score
write.table(x=Zscore.nonabsolute,file=file.path(path,"Zscore.tsv"),row.names=F)


# --------------- FIND MDP FOR GROUPS --------- ####
# Calculate gMDP score
Zgroups.list <- lapply(idx, function(x) rowSums(Zscore[,x])/length(x)) # calculates gMDP for each gene in each class
Zgroups <- abs(matrix(unlist(Zgroups.list), ncol=nGroups))   #
colnames(Zgroups) <- paste("gMDP",names(idx),sep="_")
rownames(Zgroups) <- rownames(Zscore)
Zgroups.annotated <- data.frame("Symbol" = Zscore.annotated[,"Symbol"], Zgroups)

gMDP.ref <- Zgroups.annotated[,paste("gMDP_",control_lab,sep="")]

# Calculate the frequency of perturbation
Zgroups.np.list <- lapply(idx, function(x) rowSums((Zscore[,x] > 0))/length(x)) # calculates gMDP for each gene in each class
Zgroups.np <- abs(matrix(unlist(Zgroups.np.list), ncol=nGroups))   #
colnames(Zgroups.np) <- paste("freq",names(idx),sep="_")
rownames(Zgroups.np) <- rownames(Zscore)
Zgroups.annotated <- data.frame(Zgroups.annotated,Zgroups.np)

gMDP.np.ref <- Zgroups.np[,paste("freq_",control_lab,sep="")]

# Calculate the score -> (gMDP - gMDPcontrol)*(freq - freqcontrol), if (g - gc) and (f - fc) > 0
gMDP.score <- Zgroups - gMDP.ref
np.score <- Zgroups.np - gMDP.np.ref

# make data frame of normalised gMDP and freq scores
gene.df <- data.frame()
for (i in 1:dim(gMDP.score)[2]){
  gene.df <- rbind(gene.df,data.frame("Gene" = rownames(gMDP.score), "gMDP" = gMDP.score[,i], "Freq" = np.score[,i], "Disease" = gsub("gMDP_","",colnames(gMDP.score)[i])))
}
gene.df <- gene.df[!(gene.df$Disease == control_lab),]

gMDP.score[gMDP.score < 0] <- 0
np.score[np.score < 0] <- 0
score <- gMDP.score*np.score
colnames(score) <- paste("score",names(idx),sep="_")

Zgroups.annotated <- data.frame(Zgroups.annotated,score)
Zgroups.annotated <- Zgroups.annotated[,-grep(paste("score_",control_lab,sep=""),names(Zgroups.annotated))] # remove control score
Zgroups.annotated <- Zgroups.annotated[order(-rowSums(score)),]
Zgroups.ranked <- cbind("order" = seq(1:nrow(Zgroups.annotated)),Zgroups.annotated)

# find the top perturbed genes
genes <- Zgroups.ranked$Symbol[rowSums(score) > 0]
genes <- genes[1:round(length(genes)/4)]


# --------------- FIND sMDP FOR ALL GENE SETs --------------------#####
# geneset1 = all genes, geneset2 = perturbed genes, rest = pathways
genesets <- list()
genesets[[1]] <- rownames(Zscore)
genesets[[2]] <-  genes
names(genesets)[1:2] <- c("allgenes","perturbedgenes")
if (!missing(pathways)){
  genesets <- c(genesets,gmt)
}
names(genesets)[1:2] <- c("allgenes","perturbedgenes")


# ------ find sMDP scores for gene subsets ---------------#####
progress("Calculating sMDP scores")

Zsamples.list <- data.frame()
sMDP.list <- list()
for(s in 1:length(genesets)){

  Zscore.annotated.sub <- Zscore.annotated[genesets[[s]],]
  Zsamples <- colSums(Zscore.annotated.sub[,-c(1)], na.rm=T)/nrow(Zscore.annotated.sub)

  Zsamples.df <- data.frame("Sample" = names(Zsamples), "sMDP" = Zsamples, "Class" = pdata[names(Zsamples),"Class"], "Perturbed"  = "Unperturbed",stringsAsFactors=FALSE)
  Zsamples.list <- rbind(Zsamples.list,Zsamples.df$sMDP)
  # find mean and standard deviation of healthy samples
  healthy.mean <- GetMeanSD(Zsamples[idx[[control_idx]]],measure)
  Zsamples.df$Thresh <- rep(healthy.mean[1] + healthy.mean[[2]],length(Zsamples))

  # leave one out for healthy
  leave.out.mean <- sapply(1:length(idx[[control_idx]]), function(x,y) (GetMeanSD(y[idx[[control_idx]][-x],"sMDP"],measure)), y=Zsamples.df)
  Zsamples.df[idx[[control_idx]],"Thresh"] <- leave.out.mean[1,] + healthy.mean[2]

  # mark perturbed samples
  Zsamples.df[Zsamples.df$sMDP >= Zsamples.df$Thresh,"Perturbed"] <- "Perturbed"
  sMDP.list[[s]] <- Zsamples.df

}

names(Zsamples.list) <- Zsamples.df$Sample
names(sMDP.list) <- names(genesets)


# --------------- PLOT: gMDP FOR ALL GENES and genesets--------------------####


progress("Generating gMDP scores")


  # write table gMDP

  write.table(x=Zgroups.annotated,file=file.path(path,"gMDPscore.tsv"),row.names=F)

if (print == TRUE){

    Zgroups.ranked <- cbind("order" = seq(1:nrow(Zgroups.annotated)),Zgroups.annotated)
    # print frequency versus gMDP score for all genes
    ngenes <- nrow(Zgroups.ranked)
    Zgroups.plot <- data.frame(matrix(nrow = nGroups*ngenes, ncol = 4)) # construct a matrix of dim (ngenes*classes)*3, cols show gene Order, gMDP, Class,  frac perturbed

    for (i in 1:nGroups){
      Zgroups.plot[(ngenes*(i-1)+1):(ngenes*i),1:3] <- Zgroups.ranked[,c(1,i+2, i+2+nGroups)]
      Zgroups.plot[(ngenes*(i-1)+1):(ngenes*i),4] <- rep(colnames(Zgroups.ranked)[i+2],ngenes)
    }
    names(Zgroups.plot) <- c("Order","gMDP","fracPerturbed","Class")
    head(Zgroups.plot)
    Zgroups.plot$Class <- factor(Zgroups.plot$Class)

    pdf(file.path(path,"geneMDPfreq.pdf"))
    test <-  ggplot2::ggplot(Zgroups.plot, ggplot2::aes(x=fracPerturbed, y=gMDP, alpha=1, colour=Class, group=Class)) +  ggplot2::geom_jitter(alpha=0.3, size=3) +
     ggplot2::theme_bw() + #+  #geom_line(ggplot2::aes(linetype = factor(State)))
       #geom_smooth(alpha=.2, size=1) + ggplot2::theme_bw() + #+  #geom_line(ggplot2::aes(linetype = factor(State)))
     ggplot2::labs(title="gMDP vs. fraction of samples gene perturbed in", x="Fraction of samples in which gene is perturbed", y = "gMDP score") + ggplot2::guides(alpha = FALSE)
    plot(test)
    dev.off()

    # print gMDP score
 #   pdf(file.path(folder.path,"geneMDP.pdf"))
#    test2 <-  ggplot2::ggplot(Zgroups.plot, ggplot2::aes(x=Order, y=gMDP, alpha=1, colour=Class, group=Class)) +  ggplot2::geom_jitter(alpha=0.2, size=3) +
 #    ggplot2::theme_bw() + #+  #geom_line(ggplot2::aes(linetype = factor(State)))
  #     #geom_smooth(alpha=.2, size=1) + ggplot2::theme_bw() + #+  #geom_line(ggplot2::aes(linetype = factor(State)))
  #   ggplot2::labs(title="gMDP for each class", x="Genes", y = "gMDP score") + ggplot2::guides(alpha = FALSE)
  #  plot(test2)
  #  dev.off()

}

# ------------------------- GENE SET ANALYSIS for pathways ------------ #####


progress("Ranking genesets")

 	# Calculate the number of genes in each pathway
	geneNumber <- lapply(1:length(genesets), function(x) length(genesets[[x]]))
	geneNumber <- unlist(geneNumber)

	overlapNumber <- lapply(1:length(genesets), function(x,y) sum(y %in% genesets[[x]]), y=Zscore.annotated$Symbol)
	overlapNumber <- unlist(overlapNumber)


	pathwayMDP.df <- data.frame("Pathway"= names(genesets),"Rank"=rep(NA,length(genesets)), "Test.value"=rep(NA,length(genesets)), "Pathway.size" = geneNumber, "Genes.in.data" = overlapNumber, Zsamples.list)

	smdp.h <- as.matrix(Zsamples.list[,pdata$Sample[pdata$Class == control_lab]])
	smdp.p <- as.matrix(Zsamples.list[,pdata$Sample[pdata$Class != control_lab]])



 	pathwayMDP.df$Test.value <-  sapply(1:(dim(smdp.h)[1]),function(x) (mean(smdp.p[x,])-mean(smdp.h[x,]))/(sd(smdp.h[x,], na.rm=T) + sd(smdp.p[x,], na.rm = T)))

	pathwayMDP.df$Rank <- rank(pathwayMDP.df$Test.value,ties.method="min", na.last = "keep")
  pathwayMDP.df <- pathwayMDP.df[order(-pathwayMDP.df$Rank),]

  # write table
  write.table(x=pathwayMDP.df,file=file.path(path,"sMDPscores.tsv"),row.names=F)

  # make plot for Pathway information


  pathwayMDP.df$Pathway <- factor(pathwayMDP.df$Pathway, levels = pathwayMDP.df$Pathway[order(pathwayMDP.df$Test.value)])


  pdf(file.path(path,"geneset_summary.pdf"))
  test <- ggplot2::ggplot(pathwayMDP.df, ggplot2::aes(x=Pathway, y=Test.value)) +  ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_bw() + ggplot2::coord_flip() +
    ggplot2::labs(title="Pathway summary", x="Pathways", y = "Signal to noise ratio of control versus non-control sample sMDP scores")
  plot(test)
  dev.off()




  if (print == TRUE){
    # --------------- PLOT: sMDP -------------------#####

    top.pathway <- pathwayMDP.df$Pathway[-match(c("allgenes","perturbedgenes"),pathwayMDP.df$Pathway)][1]
    pathway.names <- c("allgenes","perturbedgenes",as.character(top.pathway))
    pathway.names.idx <- match(pathway.names, names(sMDP.list))
    pathway.names <- make.names(pathway.names)

    for (j in 1:3){
      smdp.plot(sMDP.data=sMDP.list[[pathway.names.idx[j]]],filename=pathway.names[j],directory=path,title.graph=pathway.names[j])
    }

  }


# ---------------- OUTPUT ---------------- ####
output <-list(Zscore.annotated,Zgroups.annotated,sMDP.list,pathwayMDP.df)
names(output) <- c("Zscore","gMDP","sMDP","pathways")



return(output)
}




#' Read gmt
#'
#' File to read a gmt file and return as a list
#' @export
#' @param fname the path to your gene matrix transposed file format (*.gmt) file
read.gmt <- function(fname){
  res <- list(genes=list(),
              desc=list())
  gmt <- file(fname)
  gmt.lines <- readLines(gmt)
  close(gmt)
  gmt.list <- lapply(gmt.lines,
                     function(x) unlist(strsplit(x, split="\t")))
  gmt.names <- sapply(gmt.list, '[', 1)
  gmt.desc <- lapply(gmt.list, '[', 2)
  gmt.genes <- lapply(gmt.list,
                      function(x){x[3:length(x)]})
  names(gmt.desc) <- names(gmt.genes) <- gmt.names
  return(gmt.genes)
}



#' Plot sMDP scores
#'
#' Plots the sMDP scores for a given geneset. Takes a data.frame that contains the sMDP scores for one geneset, along with the Sample and Class information. Data.frame must have a sMDP, Sample and Class columns.
#' @export
#' @param sMDP.data a data.frame containing sMDP information for a geneset, with columns "Sample", "sMDP" and "Class"
#' @param filename filename
#' @param directory directory to save file
#' @param title.graph title name for graph
smdp.plot <- function(sMDP.data,filename="",directory="",title.graph=""){

  if (directory != ""){
    path = directory
    if (dir.exists(directory) == F){
      dir.create(directory)
    }
  } else {
    path = "."
  }

  # -------- Multiple plot function ---- ####
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot2::ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
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


  Zsamples.plot <-sMDP.data[order(sMDP.data$sMDP),]

  Zsamples.plot$Sample <- factor(Zsamples.plot$Sample, levels = Zsamples.plot$Sample[order(Zsamples.plot$sMDP)])

  #Plot sMDP as bar graphs
  bp1 <- ggplot2::ggplot(data = Zsamples.plot, ggplot2::aes(y = sMDP, x = Sample, fill = Class)) +
    ggplot2::geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
    ggplot2::labs(title = title.graph, x = "Samples", y = "sMDP score") +
    #  geom_hline(yintercept = MDP_cut, color = "darksalmon", linetype = "dashed") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::theme(axis.line = ggplot2::element_line(size = 0.5,
                                                     linetype = "solid"), panel.grid.major = ggplot2::element_line(colour = "black",
                                                                                                                   linetype = "blank"), panel.grid.minor = ggplot2::element_line(linetype = "blank"),
                   axis.title = ggplot2::element_text(size = 12),
                   axis.text = ggplot2::element_text(size = 12, angle = 90),
                   plot.title = ggplot2::element_text(size = 12),
                   panel.background = ggplot2::element_rect(fill = "grey100"),
                   legend.key = ggplot2::element_rect(fill = "grey85"),
                   legend.background = ggplot2::element_rect(fill = "grey94"),
                   legend.direction = "horizontal") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = NA, linetype = "solid"))


  # find order
  c=1
  meansMDP <- vector()
  for (j in unique(Zsamples.plot$Class)){
    meansMDP[c] = mean(Zsamples.plot[Zsamples.plot$Class == j,"sMDP"])
    c = c+1
  }
  names(meansMDP) <- unique(Zsamples.plot$Class)
  meansMDP <- meansMDP[order(meansMDP)]


  ##Boxplot of sMDP score for each class
  bp2 <- ggplot2::ggplot(data = Zsamples.plot, ggplot2::aes(y = sMDP, x = Class, fill = Class)) +
    ggplot2::geom_boxplot(outlier.shape=NA) + ggplot2::stat_summary(fun.y = mean, geom = "point", shape = 23, size = 6) +
    ggplot2::labs(title = title.graph, x = "Groups", y = "sMDP score") +
    ggplot2::theme(legend.position = "null") +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
    #ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), ggplot2::aes(color=Class)) +
    ggplot2::scale_x_discrete(limits=names(meansMDP)) +
    ggplot2::geom_jitter(shape = 16, position = ggplot2::position_jitter(0.2), size=2, color = "grey10", alpha=0.7) +
    ggplot2::theme(axis.line = ggplot2::element_line(size = 0.5, linetype = "solid"),
                   panel.grid.major = ggplot2::element_line(linetype = "blank"),
                   panel.grid.minor = ggplot2::element_line(linetype = "blank"),
                   panel.background = ggplot2::element_rect(fill = "white"),
                   axis.text = ggplot2::element_text(size = 10))

  smdp.name <- paste(filename,"sMDP.pdf",sep="")
  pdf(file.path(path,smdp.name))
  multiplot(bp1,bp2,cols=2)
  dev.off()


}











