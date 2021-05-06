#!/usr/bin/env Rscript

# RUN:
# transphylo_inferTTree.R -s <study_accession> -t <tree_file> -o <out_dir> -c <clusters_file> -p <processors>

# Setup ----

rm(list=ls())


library(optparse)
library(ape)
library(TransPhylo)
library(foreach)
library(doParallel)
library(parallel)



option_list = list(
  make_option(c("-s", "--study_accession"), type="character", default=NULL,
              help="enter study accession code e.g. PAKISTAN", metavar="character"),
  make_option(c("-t", "--tree_file"), type="character", default=NULL,
              help="input (dated) tree file - i.e. MCC tree from BEAST and TreeAnnotator", metavar="character"),
  make_option(c("-c", "--clusters_file"), type="character", default=NULL,
              help="input file of sample names and their clusters from cut_tree.py", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL,
              help="name of transphylo output directory - outputs files and plots", metavar="character"),
  make_option(c("-m", "--mcmc_iter"), type="character", default=1000,
              help="number of mcmc iterations of TransPhylo", metavar="character"), 
  make_option(c("-p", "--processors"), type="character", default=4,
              help="number of processors/cores to use", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Input args ----

study_accession <- opt$study_accession
tree_file <- opt$tree_file
clusters_file <- opt$clusters_file
out_dir <- opt$out_dir
mcmc_iter <- as.numeric(opt$mcmc_iter)
cores <- as.numeric(opt$processors)

print("ARGUMENTS:")
print(c("study_accession: ", study_accession))
print(c("tree_file: ", tree_file))
print(c("clusters_file: ", clusters_file))
print(c("out_dir: ", out_dir))
print(c("mcmc_iter: ", mcmc_iter))
print(c("processors/cores: ", cores))



# study_accession <- "PAKISTAN_ALL"
# tree_file <- "beast_results/PAKISTAN_ALL.mcc.tree"
# clusters_file <- "metadata/PAKISTAN_ALL.clusters.csv"
# out_dir <- "transphylo_results/"
# mcmc_iter <- 1000


# Other files ----

# Trees 
trees_file <- paste0(out_dir, study_accession, ".ptree_list.")

# Plots
transphylo_res_plot_file <- paste0(out_dir, study_accession, ".transphylo.pdf")
medoid_tree_plot_file <- paste0(out_dir, study_accession, ".med_tree.pdf")
trans_tree_plot_file <- paste0(out_dir, study_accession, ".trans_tree.pdf")
mat_prob_plot_file <- paste0(out_dir, study_accession, ".mat_prob.pdf")
mat_pair_plot_file <- paste0(out_dir, study_accession, ".mat_pair.pdf")
cases_plot_file <- paste0(out_dir, study_accession, ".cases.pdf")
gen_time_plot_file <- paste0(out_dir, study_accession, ".gen_time.pdf")
samp_time_plot_file <- paste0(out_dir, study_accession, ".samp_time.pdf")
inf_time_plot_file <- paste0(out_dir, study_accession, ".inf_time.pdf")
offspring_plot_file <- paste0(out_dir, study_accession, ".offspring.pdf")

# Tables
es_table_file <- paste0(out_dir, study_accession, ".es_table.csv")


# TransPhylo parameters ----

# Parameters of Gamma distr representing generation time
w.shape <- 2.2
w.scale <- 2.1


# Load in data ---- 

# Load and parse MCC tree from BEAST
phy <- ape::read.nexus(tree_file)

# Load in clusters
clusters <- read.csv(clusters_file, header = F, col.names = c("clust", "id"))


# Prep data for transphylo ----

# Get last date from tree
last_date_all <- as.numeric(max(unlist(lapply(strsplit(phy$tip.label, "_"), function(x){x[length(x)]}))))

# Split clusters table to get clusters to loop over
clusters_split <- split(clusters, clusters[, 1])

# Split the tree into those clusters
tree_list <- list()
for (i in seq(clusters_split)){
  tree_list[[i]] <- keep.tip(phy, clusters_split[[i]][,2])
}

# Get last dates from tree list
last_date_list <- list()
for (i in seq(tree_list)){
  last_date_list[[i]] <- as.numeric(max(unlist(lapply(strsplit(tree_list[[i]]$tip.label, "_"), function(x){x[length(x)]}))))
}

# Convert ape phylo object into phylogenetic tree. Pass last date.
# ptree <- ptreeFromPhylo(phy, dateLastSample=last_date)
ptree_list <- list()
for (i in seq(tree_list)){
  ptree_list[[i]] <- ptreeFromPhylo(tree_list[[i]], dateLastSample=last_date_list[[i]])
}



# Transphylo ----
# https://xavierdidelot.github.io/TransPhylo/articles/infer.html

# Run transphylo model ----

# Infer transmission tree given phylogenetic tree
# n.b. dateT parameter - Time at which observation of cases stopped
# Need to add small number greater than 1e-10 - bug in the inferTTree() function 
# See - https://github.com/xavierdidelot/TransPhylo/blob/master/R/inferTTree.R
# res <- inferTTree(ptree, mcmcIterations = mcmc_iter, w.shape=w.shape, w.scale=w.scale, dateT=dateT)

# cores <- parallel::detectCores()
cl <- parallel::makeCluster(cores) 
doParallel::registerDoParallel(cl)

res_list <- foreach::foreach(i = 16) %dopar% {
  TransPhylo::inferTTree(ptree_list[[i]], 
                         # Shape parameter of the Gamma distribution representing the generation time
                         w.shape = 2.2,
                         # Scale parameter of the Gamma distribution representing the generation time
                         w.scale = 2.1,
                         # Mean of the Gamma distribution representing the generation time
                         w.mean = 3.9,
                         # Std of the Gamma distribution representing the generation time
                         w.std = NA,
                         # Shape parameter of the Gamma distribution representing the sampling time
                         ws.shape = 2.2,
                         # Scale parameter of the Gamma distribution representing the sampling time
                         ws.scale = 2.1,
                         # Mean of the Gamma distribution representing the sampling time
                         ws.mean = 3.9,
                         # Std of the Gamma distribution representing the sampling time
                         ws.std = NA,
                         # Starting value of within-host coalescent parameter Ne*g
                         startNeg = 1.48,
                         # Whether of not to update the parameter Ne*g
                         updateNeg = TRUE,
                         # Starting value of parameter off.r
                         startOff.r = 1,
                         # Whether or not to update the parameter off.r
                         updateOff.r = TRUE,
                         # Starting value of parameter off.p
                         startOff.p = 0.5,
                         # Whether or not to update the parameter off.p
                         updateOff.p = FALSE,
                         # Starting value of sampling proportion pi
                         startPi = 0.01,
                         # Whether or not to update the parameter pi
                         updatePi = TRUE,
                         # Optional combined tree to start from  
                         startCTree = NA,
                         # Whether or not to update the transmission tree
                         updateTTree = TRUE,
                         # Type of optimisation to apply to MCMC start point (0=none, 1=slow, 2=fast)
                         optiStart = 2,
                         # Date when process stops (this can be Inf for fully simulated outbreaks)
                         dateT = last_date_list[[i]]+runif(1)*1e-08,
                         # Number of MCMC iterations to run the algorithm for
                         mcmcIterations = mcmc_iter,
                         # MCMC thinning interval between two sampled iterations   
                         thinning = mcmc_iter * 0.01,
                         # Whether or not to use verbose mode (default is false)
                         verbose = T)
}
#stop cluster
parallel::stopCluster(cl)


save(res_list, file = paste0(out_dir, "TP_cluster_16_test.Rdata"))


















