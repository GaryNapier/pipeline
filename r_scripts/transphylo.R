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


# Customise Transphylo plot functions

plotCTree_x <- function (tree, showLabels = TRUE, showStars = TRUE, cols = NA, maxTime = NA, cex = 1) {
  nam = tree$nam
  tree = tree$ctree
  nsam <- sum(tree[, 2] + tree[, 3] == 0)
  nh <- nrow(tree) - 3 * nsam + 1
  ntot <- nsam + nh
  oldpar <- par("yaxt", "bty", "xpd")
  on.exit(par(oldpar))
  # par(yaxt = "n", bty = "n", xpd = T)
  par(yaxt = "n", bty = "n", xpd = T, cex.axis = 0.75)
  
  # NEW
  time_range <- c(min(tree[, 1]), max(tree[, 1]))
  max_x <- time_range[2] + ceiling(length(seq(time_range[1], time_range[2]))*0.1)
  max_x_maxTime <- maxTime + ceiling(length(seq(time_range[1], maxTime))*0.1)
  
  # plot(0, 0, type = "l", 
  #      xlim = c(min(tree[, 1])-1, 
  #               ifelse(is.na(maxTime), max(tree[, 1]), maxTime)), 
  #      ylim = c(0, nsam + 1), xlab = "", ylab = "")
  
  plot(0, 0, type = "l", 
       xlim = c(min(time_range)-1, 
                ifelse(is.na(maxTime), 
                       max_x, 
                       max_x_maxTime)), 
       ylim = c(0, nsam + 1), xlab = "", ylab = "")
  
  host <- tree[, 4]
  
  if (ntot > 1) {
    if (is.na(cols[1])) 
      grDevices::palette(grDevices::rainbow(min(1024, ntot)))
    else grDevices::palette(cols)
  }
  
  root <- which(host == 0)
  ys <- matrix(0, nsam, 1)
  todo <- cbind(root, 0, 0.5, 1)
  while (nrow(todo) > 0) {
    w <- todo[1, 1]
    x <- todo[1, 2]
    y <- todo[1, 3]
    scale <- todo[1, 4]
    if (tree[w, 2] == 0 && tree[w, 3] == 0) {
      ys[w] <- y
    }
    else if (tree[w, 3] == 0) {
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1], 
                                y, scale, deparse.level = 0))
    }
    else {
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1], 
                                y + scale/2, scale/2, deparse.level = 0), 
                    cbind(tree[w, 3], tree[w, 1], y - scale/2, scale/2, deparse.level = 0))
    }
    todo <- rbind(todo[-1, ])
  }
  
  ys <- rank(ys)
  for (i in ((nsam + 1):nrow(tree))) {
    children <- c()
    todo <- i
    while (length(todo) > 0) {
      children = c(children, todo[1])
      todo = c(todo[-1], setdiff(tree[todo[1], 2:3], 0))
    }
    ys[i] <- mean(ys[children[which(children <= nsam)]])
  }
  
  todo <- cbind(root, tree[root, 1])
  while (nrow(todo) > 0) {
    w <- todo[1, 1]
    x <- todo[1, 2]
    y <- ys[w]
    col = host[w]
    if (tree[w, 2] == 0 && tree[w, 3] == 0) {
      lines(c(x, tree[w, 1]), c(y, y), col = col, lwd = 2)
      if (showLabels) 
        text(tree[w, 1], y, nam[w], cex = cex, pos = 4)
    }
    else if (tree[w, 3] == 0) {
      lines(c(x, tree[w, 1]), c(y, y), col = col, lwd = 2)
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1]))
    }
    else {
      lines(c(x, tree[w, 1]), c(y, y), col = col, lwd = 2)
      lines(c(tree[w, 1], tree[w, 1]), cbind(ys[tree[w, 2]], ys[tree[w, 3]]), col = col, lwd = 2)
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1]), 
                    cbind(tree[w, 3], tree[w, 1]))
    }
    todo <- rbind(todo[-1, ])
  }
  todo <- cbind(root, tree[root, 1])
  while (nrow(todo) > 0 && showStars) {
    w <- todo[1, 1]
    x <- todo[1, 2]
    y <- ys[w]
    col = host[w]
    if (tree[w, 2] == 0 && tree[w, 3] == 0) {
    }
    else if (tree[w, 3] == 0) {
      points(tree[w, 1], y, col = "red", pch = 8)
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1]))
    }
    else {
      todo <- rbind(todo, cbind(tree[w, 2], tree[w, 1]), 
                    cbind(tree[w, 3], tree[w, 1]))
    }
    todo <- rbind(todo[-1, ])
  }
  return(invisible(tree))
}

plotTTree2_x <- function (ttree, showLabels = TRUE, showMissingLinks = 0, cex = 1) {
  nam = ttree$nam
  ttree = ttree$ttree
  ttree = cbind(ttree, rep(1, nrow(ttree)))
  if (showMissingLinks > 0) {
    i = which(is.na(ttree[, 2]))[1]
    while (i < nrow(ttree)) {
      w = which(ttree[, 3] == i)
      if (length(w) == 1) {
        ttree[w, 3] = ttree[i, 3]
        ttree[w, 4] = ttree[w, 4] + ttree[i, 4]
        ttree = ttree[-i, ]
        ttree[which(ttree[, 3] > i), 3] = ttree[which(ttree[, 3] > i), 3] - 1
      }
      else i = i + 1
    }
  }
  if (showMissingLinks == 2) {
    ttree[which(ttree[, 4] >= 2), 4] = 2
  }
  n = nrow(ttree)
  ys <- rep(0, n)
  scale <- rep(1, n)
  todo = c(which(ttree[, 3] == 0))
  while (length(todo) > 0) {
    f = which(ttree[, 3] == todo[1])
    o = rank(-ttree[f, 1])
    f[o] = f
    for (i in f) {
      ys[i] = ys[todo[1]] + scale[todo[1]] * which(f == i)/(length(f) + 1)
      scale[i] = scale[todo[1]]/(length(f) + 1)
      todo = c(todo, i)
    }
    todo = todo[-1]
  }
  ys = rank(ys)
  oldpar <- par("yaxt", "bty")
  on.exit(par(oldpar))
  par(yaxt = "n", bty = "n", cex.axis = 0.75)
  mi = min(ttree[which(!is.na(ttree[, 1])), 1])-1
  ma = max(ttree[which(!is.na(ttree[, 1])), 1])+2
  plot(c(), c(), 
       # xlim = c(mi - (ma - mi) * 0.05, ma + (ma - mi) * 0.05), 
       xlim = c( (mi - (ma - mi) * 0.05), 
                 (ma + (ma - mi) * 0.15) ), 
       ylim = c(0, n + 1), 
       xlab = "", ylab = "")
  pal = grDevices::gray.colors(max(ttree[, 4]))
  for (i in 1:n) {
    if (ttree[i, 3] != 0) {
      dircol = pal[ttree[i, 4]]
      arrows(ttree[ttree[i, 3], 1], ys[ttree[i, 3]], ttree[i, 
                                                           1], ys[i], length = 0, col = dircol)
    }
    if (showLabels && !is.na(ttree[i, 2])) 
      text(ttree[i, 1], ys[i], nam[i], pos = 4, cex = cex)
  }
  for (i in 1:n) {
    points(ttree[i, 1], ys[i], pch = 21, bg = ifelse(is.na(ttree[i, 
                                                                 2]), "white", "black"), cex = cex)
  }
  if (length(pal) > 2) 
    legend("topleft", legend = 0:(length(pal) - 1), col = pal, 
           lty = 1, cex = cex, title = "Missing links")
  return(invisible(ttree))
}

getInfectionTimeDist_x <- function(record, burnin = 0.5, k, numBins = 10, show.plot = F) {
  record = record[max(1, round(length(record) * burnin)):length(record)]
  for (i in 1:length(record)) record[[i]] = extractTTree(record[[i]]$ctree)
  times = matrix(NA, length(k), length(record))
  for (i in 1:length(k)) for (j in 1:length(record)) {
    ii = which(record[[j]]$nam == k[i])
    times[i, j] = record[[j]]$ttree[ii, 1]
  }
  if (show.plot) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    # par(mfrow = c(length(k), 1))
    scale <- 0.75
    par(mfrow = n2mfrow(length(k)), 
        cex.axis = scale, 
        cex.main = scale, 
        cex.lab = scale)
    xlim = c(min(times), max(times))
    uni = length(unique(as.vector(times)))
    numBins = min(numBins, uni)
    br = seq(xlim[1], xlim[2], length.out = numBins + 1)
    for (i in 1:length(k)) {
      h = hist(times[i, ], breaks = br, plot = F)$counts/length(record)
      # barplot(h, main = "", xlab = "", ylab = sprintf("Infection time of %s", k[i]))
      barplot(h, main = k[i], xlab = "", ylab = "Infection time")
      if (xlim[1] == xlim[2]) 
        axis(1, at = 0.7, labels = xlim[1])
      else {
        labs = pretty(xlim, 6)
        axis(1, at = (labs - xlim[1])/(xlim[2] - xlim[1]) * 1.2 * numBins, labels = labs)
      }
    }
  }
  if (length(k) == 1) 
    times = as.vector(times)
  return(times)
}

getOffspringDist_x <- function (record, burnin = 0.5, k, show.plot = F) {
  record = record[max(1, round(length(record) * burnin)):length(record)]
  for (i in 1:length(record)) record[[i]] = extractTTree(record[[i]]$ctree)
  offspring = matrix(NA, length(k), length(record))
  for (i in 1:length(k)) for (j in 1:length(record)) {
    ii = which(record[[j]]$nam == k[i])
    offspring[i, j] = as.numeric(length(which(record[[j]]$ttree[, 3] == ii)))
  }
  if (show.plot) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    # par(mfrow = c(length(k), 1))
    scale <- 0.75
    par(mfrow = n2mfrow(length(k)), 
        cex.axis = scale, 
        cex.main = scale, 
        cex.lab = scale)
    xlim = c(-0.5, max(offspring) + 0.5)
    br = seq(xlim[1], xlim[2])
    for (i in 1:length(k)) {
      h = hist(offspring[i, ], breaks = br, plot = F)$counts/length(record)
      # barplot(h, main = "", xlab = "", ylab = sprintf("Offspring of %s", 
      #                                                 k[i]), names.arg = 0:max(offspring))
      barplot(h, main = k[i], xlab = "", 
              ylab = "Offspring", names.arg = 0:max(offspring))
    }
  }
  if (length(k) == 1)
    offspring = as.vector(offspring)
  return(offspring)
}



# Variables - setup ----

# Plots
plot_text_sz <- 0.5


# Arguments ----

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
              help="number of mcmc iterations of TransPhylo", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Input args ----

study_accession <- opt$study_accession
tree_file <- opt$tree_file
clusters_file <- opt$clusters_file
out_dir <- opt$out_dir
mcmc_iter <- as.numeric(opt$mcmc_iter)

print("ARGUMENTS:")
print(c("study_accession: ", study_accession))
print(c("tree_file: ", tree_file))
print(c("clusters_file: ", clusters_file))
print(c("out_dir: ", out_dir))
print(c("mcmc_iter: ", mcmc_iter))

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

# # Save ptrees
# for (i in seq(ptree_list)){
#   save(trees_file, paste0(trees_file, i, ".Rdata") )
# }









# Transphylo ----
# https://xavierdidelot.github.io/TransPhylo/articles/infer.html

# Run transphylo model ----

# Infer transmission tree given phylogenetic tree
# n.b. dateT parameter - Time at which observation of cases stopped
# Need to add small number greater than 1e-10 - bug in the inferTTree() function 
# See - https://github.com/xavierdidelot/TransPhylo/blob/master/R/inferTTree.R
# res <- inferTTree(ptree, mcmcIterations = mcmc_iter, w.shape=w.shape, w.scale=w.scale, dateT=dateT)

cores <- parallel::detectCores()
cl <- parallel::makeCluster(cores[1]-1) #not to overload your computer
doParallel::registerDoParallel(cl)

# res_list <- list()
# for(i in seq(ptree_list)){
#   res_list[[i]] <- inferTTree(ptree_list[[i]], 
#                               mcmcIterations = mcmc_iter, 
#                               w.shape = w.shape, w.scale = w.scale, 
#                               dateT = last_date_list[[i]]+runif(1)*1e-08)
# }

res_list <- foreach::foreach(i = seq(ptree_list)) %dopar% {
  inferTTree(ptree_list[[i]], mcmcIterations = mcmc_iter, 
             w.shape = w.shape, w.scale = w.scale, 
             dateT = last_date_list[[i]]+runif(1)*1e-08)
}
#stop cluster
parallel::stopCluster(cl)







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

es_table <- do.call("rbind", es_list)
write.csv(es_table, es_table_file, header = T)


# Interpretation of output ----

# Find the most representative (aka medoid) colored tree
# med <- medTTree(res)
# plot(med, cex = plot_text_sz)

med_list <- list()
pdf(medoid_tree_plot_file)
for (i in seq(res_list)){
  med_list[[i]] <- medTTree(res_list[[i]])
  plotCTree_x(med_list[[i]], cex = plot_text_sz, maxTime = last_date_list[[i]]+1)
}
dev.off()

# Plot corresponding transmission tree:
# ttree <- extractTTree(med)
# plot(ttree,type='detailed',w.shape,w.scale, cex = plot_text_sz)

ttree_list <- list()
pdf(trans_tree_plot_file)
for(i in seq(med_list)){
  ttree_list[[i]] <- extractTTree(med_list[[i]])
  plotTTree2_x(ttree_list[[i]], cex = plot_text_sz)
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
  inf_time_list[[i]] <- getInfectionTimeDist_x(res_list[[i]], k = clusters_split[[i]]$id, show.plot = T)
}
dev.off()

# Offspring distribution for the individuals labelled ‘1’ and ‘2’:
# a <- getOffspringDist(res,k=c('1','2'),show.plot = T)

offspring_list <- list()
pdf(offspring_plot_file)
for(i in seq(res_list)){
  offspring_list[[i]] <- getOffspringDist_x(res_list[[i]], k = clusters_split[[i]]$id, show.plot = T)
}
dev.off()


















