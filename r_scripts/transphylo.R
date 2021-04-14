#!/usr/bin/env Rscript

# <file_name>.R


# Description
# 

# Arguments to script:

# Input:

# Main steps:

# Output:

# RUN:
# Rscript r_scripts/<file_name>.R
# Rscript r_scripts/<file_name>.R

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

# Variables 

# Plots
units <- "px"
wth <- 3000
ht <- 2000
resolution <- 300
plot_text_sz <- 0.5

# Arguments ----

# option_list = list(
#   # EXAMPLE:
#   # make_option(c("-t", "--template_file_name"), type="character", default=NULL,
#   #             help="input template xml file", metavar="character"),
#   
#   make_option(c("-s", "--study_accession"), type="character", default=NULL,
#               help="enter study accession code e.g. PAKISTAN", metavar="character"),
#   make_option(c("-p", "--plots_path"), type="character", default="./", 
#               help="path to save plots to; default is ./", metavar="character"),
#   make_option(c("-t", "--tree_file"), type="character", default=NULL, 
#               help="input (dated) tree file - i.e. MCC tree from BEAST and TreeAnnotator", metavar="character"), 
#   make_option(c("-o", "--out_file"), type="character", default=NULL, 
#               help="name of transphylo output file", metavar="character")
# ); 

# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);

# print("ARGUMENTS:")
# print(opt)
# print("---")
# print(str(opt))


# Paths ----

# plots_path <- opt$plots_path
plots_path <- "plots/"

# Files ----

# study_accession <- opt$study_accession
# tree_file <- opt$tree_file
# out_file <- opt$out_file
study_accession <- "THAILAND_TEST"
tree_file <- "beast_results/THAILAND_TEST.manual.mcc_tree"
out_file <- "transphylo_results/THAILAND_TEST.transphylo"

mcc_plot_file <- paste0(plots_path, basename(tree_file), ".png")
transphylo_res_plot_file <- paste0(plots_path, basename(tree_file), ".transphylo.png")

print("ARGS:")
print(c("study_accession: ", study_accession))
print(c("tree_file: ", tree_file))
print(c("out_file: ", out_file))
print(c("mcc_plot_file", mcc_plot_file))


# Transphylo ----

# https://xavierdidelot.github.io/TransPhylo/articles/infer.html


# Load and parse MCC tree from BEAST
# phy <- ape::read.tree(tree_file)
phy <- ape::read.nexus(tree_file)

# Get last date from tree
last_date <- as.numeric(max(unlist(lapply(strsplit(phy$tip.label, "_"), function(x){x[length(x)]}))))

# Convert ape phylo object into phylogenetic tree. Pass last date.
ptree <- ptreeFromPhylo(phy, dateLastSample=last_date)

# Plot and save
png(mcc_plot_file, units=units, width=wth, height=ht, res = resolution)
plot(ptree, cex = plot_text_sz)
dev.off()

# Parameters ----

# Parameters of Gamma distr representing generation time
w.shape <- 2.2
w.scale <- 2.1
# Time at which observation of cases stopped
# Need to add small number greater than 1e-10 - see https://github.com/xavierdidelot/TransPhylo/blob/master/R/inferTTree.R
dateT <- last_date+runif(1)*1e-09 
mcmc_iter <- 10000

# Run transphylo model ----

# Infer transmission tree given phylogenetic tree
res <- inferTTree(ptree, mcmcIterations = mcmc_iter, w.shape=w.shape, w.scale=w.scale, dateT=dateT)

# Plot to check convergence
png(transphylo_res_plot_file, units=units, width=wth, height=ht, res = resolution)
plot(res)
dev.off()

# Further assessment MCMC convergence and mixing 
# e.g. to obtain effective sample size (ESS) of paramaters, making sure that the ESS of each parameter 
# is at least 100:
mcmc <- convertToCoda(res)
effectiveSize(mcmc)


# Interpretation of output ----

# Find the most representative (aka medoid) colored tree
med <- medTTree(res)
plot(med, cex = plot_text_sz)

# Plot corresponding transmission tree:
ttree <- extractTTree(med)
plot(ttree,type='detailed',w.shape,w.scale, cex = plot_text_sz)

# Matrix of probability of direct transmission for all pairs of individuals 
mat_prob <- computeMatWIW(res)
lattice::levelplot(mat_prob, xlab='', ylab='', cex = plot_text_sz, scales=list(x=list(rot=90)))

# Plot matrix indicating for each pair of individuals how many intermediates there are in the transmission chain
mat_pair <- computeMatTDist(res)
lattice::levelplot(mat_pair, xlab='', ylab='', cex = plot_text_sz, scales=list(x=list(rot=90)))

# Plot of sampled and unsampled cases over time:
a <- getIncidentCases(res, show.plot = T)

# Distribution of realised generation times:
a <- getGenerationTimeDist(res,show.plot = T)

# Distribution of realised sampling times:
a <- getSamplingTimeDist(res,show.plot = T)

# Distribution of infection time for the individuals labelled ‘1’ and ‘2’:
a <- getInfectionTimeDist(res, k=c('1','2'), show.plot = T)

# Offspring distribution for the individuals labelled ‘1’ and ‘2’:
a <- getOffspringDist(res,k=c('1','2'),show.plot = T)





  















