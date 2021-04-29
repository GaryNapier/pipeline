#!/usr/bin/env Rscript


# RUN:
# transphylo.R -s <study_accession> -t <tree_file> -o <out_dir> -c <clusters_file>
# transphylo.R

# Setup ----

rm(list=ls())

library(optparse)
library(ape)
library(TransPhylo)
library(coda)

heaD <- function(x,...){
  head(x, ...)
}

len_str <- function(string){
  length(unlist(strsplit(string, split = "")))
}

hs <- function(x, ...){
  print(head(x, ...))
  print("---")
  str(x, ...)
}


# Variables - setup ----

# Plots
units <- "px"
wth <- 3000
ht <- 2000
resolution <- 300
plot_text_sz <- 0.5


# Arguments ----

option_list = list(
  # EXAMPLE:
  # make_option(c("-t", "--template_file_name"), type="character", default=NULL,
  #             help="input template xml file", metavar="character"),

  make_option(c("-s", "--study_accession"), type="character", default=NULL,
              help="enter study accession code e.g. PAKISTAN", metavar="character"),
  # make_option(c("-p", "--plots_path"), type="character", default="./",
  #             help="path to save plots to; default is ./", metavar="character"),
  make_option(c("-t", "--tree_file"), type="character", default=NULL,
              help="input (dated) tree file - i.e. MCC tree from BEAST and TreeAnnotator", metavar="character"),
  make_option(c("-c", "--clusters_file"), type="character", default=NULL,
              help="input file of sample names and their clusters from cut_tree.py", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL,
              help="name of transphylo output directory - outputs files and plots", metavar="character"),
  make_option(c("-m", "--mcmc_iter"), type="character", default=1000,
              help="number of mcmc iterations of TransPhylo", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Input args ----

study_accession <- opt$study_accession
tree_file <- opt$tree_file
clusters_file <- opt$clusters_file
out_dir <- opt$out_dir
mcmc_iter <- opt$mcmc_iter

print("ARGUMENTS:")
print("study_accession: ", study_accession)
print("tree_file: ", tree_file)
print("clusters_file: ", clusters_file)
print("out_dir: ", out_dir)
print("mcmc_iter: ", mcmc_iter)

# study_accession <- "PAKISTAN_ALL"
# tree_file <- "beast_results/PAKISTAN_ALL.mcc.tree"
# clusters_file <- "metadata/PAKISTAN_ALL.clusters.csv"
# out_dir <- "transphylo_results/"
# mcmc_iter <- 1000


# Other files ----

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
res_list <- list()
for(i in seq(ptree_list)){
  res_list[[i]] <- inferTTree(ptree_list[[i]], 
                              mcmcIterations = mcmc_iter, 
                              w.shape = w.shape, w.scale = w.scale, 
                              dateT = last_date_list[[i]]+runif(1)*1e-08)
}

# Plot to check convergence
# png(transphylo_res_plot_file, units=units, width=wth, height=ht, res = resolution)
# plot(res)
# dev.off()

pdf(transphylo_res_plot_file)
for (i in seq(res_list)){
  plot(res_list[[i]])
}
dev.off()

# ESS
# FRom XavierDidelot gihub page  - "Further assessment of the MCMC convergence and mixing can be obtained using the CODA package,
# for example to obtain the effective sample size (ESS) of paramaters as follows, 
# making sure that the ESS of each parameter is at least 100":
# mcmc <- convertToCoda(res)
# effectiveSize(mcmc)

mcmc_list <- list()
es_list <- list()
for (i in seq(res_list)){
  mcmc_list[[i]] <- convertToCoda(res_list[[i]])
  es_list[[i]] <- effectiveSize(mcmc_list[[i]])
}


# NEED TO DO SOMETHING WITH THIS


# Interpretation of output ----

# Find the most representative (aka medoid) colored tree
# med <- medTTree(res)
# plot(med, cex = plot_text_sz)

med_list <- list()
pdf(medoid_tree_plot_file)
for (i in seq(res_list)){
  med_list[[i]] <- medTTree(res_list[[i]])
  plot(med_list[[i]], cex = plot_text_sz, maxTime = last_date_list[[i]]+1)
}
dev.off()

# Plot corresponding transmission tree:
# ttree <- extractTTree(med)
# plot(ttree,type='detailed',w.shape,w.scale, cex = plot_text_sz)

ttree_list <- list()
pdf(trans_tree_plot_file)
for(i in seq(med_list)){
  ttree_list[[i]] <- extractTTree(med_list[[i]])
  plot(ttree_list[[i]], cex = plot_text_sz)
}
dev.off()

# Matrix of probability of direct transmission for all pairs of individuals 
# mat_prob <- computeMatWIW(res)
# lattice::levelplot(mat_prob, xlab='', ylab='', cex = plot_text_sz, scales=list(x=list(rot=90)))

mat_prob_list <- list()
pdf(mat_prob_plot_file)
for(i in seq(res_list)){
  mat_prob_list[[i]] <- computeMatWIW(res_list[[i]])
  print(lattice::levelplot(mat_prob_list[[i]], xlab='', ylab='', cex = plot_text_sz, scales=list(x=list(rot=90))))
}
dev.off()

# Plot matrix indicating for each pair of individuals how many intermediates there are in the transmission chain
# mat_pair <- computeMatTDist(res)
# lattice::levelplot(mat_pair, xlab='', ylab='', cex = plot_text_sz, scales=list(x=list(rot=90)))

mat_pair_list <- list()
pdf(mat_pair_plot_file)
for(i in seq(res_list)){
  mat_pair_list[[i]] <- computeMatTDist(res_list[[i]])
  print(lattice::levelplot(mat_pair_list[[i]], xlab='', ylab='', cex = plot_text_sz, scales=list(x=list(rot=90))))
}
dev.off()

# Plot of sampled and unsampled cases over time:
# a <- getIncidentCases(res, show.plot = T)

cases_list <- list()
pdf(cases_plot_file)
for(i in seq(res_list)){
  cases_list[[i]] <- getIncidentCases(res_list[[i]], show.plot = T)
}
dev.off()

# Distribution of realised generation times:
# a <- getGenerationTimeDist(res, show.plot = T)

gen_time_list <- list()
pdf(gen_time_plot_file)
for(i in seq(res_list)){
  gen_time_list[[i]] <- getGenerationTimeDist(res_list[[i]], show.plot = T)
}
dev.off()

# Distribution of realised sampling times:
# a <- getSamplingTimeDist(res, show.plot = T)

samp_time_list <- list()
pdf(samp_time_plot_file)
for(i in seq(res_list)){
  samp_time_list[[i]] <- getSamplingTimeDist(res_list[[i]], show.plot = T)
}
dev.off()





# Distribution of infection time for the individuals labelled ‘1’ and ‘2’:
# a <- getInfectionTimeDist(res, k=c('1','2'), show.plot = T)

inf_time_list <- list()
pdf(inf_time_plot_file)
for(i in seq(res_list)){
  inf_time_list[[i]] <- getInfectionTimeDist(res_list[[i]], k = clusters_split[[i]]$id, show.plot = T)
}
dev.off()





# Offspring distribution for the individuals labelled ‘1’ and ‘2’:
# a <- getOffspringDist(res,k=c('1','2'),show.plot = T)

offspring_list <- list()
pdf(offspring_plot_file)
for(i in seq(res_list)){
  offspring_list[[i]] <- getOffspringDist(res_list[[i]], k = clusters_split[[i]]$id, show.plot = T)
}
dev.off()




  















