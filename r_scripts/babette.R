


# Checklist for Babette
# See XML file in (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7725332/)
# Upload XML file to Beauti to get parameters:

# Check if Babette is parsing years (Tip Dates tab > Use Tip Dates = T > Auto-configure > after last "_")
#  - Tip date file and create inference model
# create_inference_model
  # tipdates_filename 
#  - MRCA prior
#  - Node at the crown age

# CLOCK MODEL
# - strict clock
# - clock rate = 0.0000001 (1.0E-7)

# PRIORS
# - Tree.t - Coalescent Constant Population
#   - Pop Size = 100 

# - popSize.t
#   - Log Normal
#   - Lower = 0
#   - Upper = 200
#   - Value = 100

# - Add Prior
#   - Taxon set label = TreePrior
#   - Distribution = Laplace Distribution
#   - Monophyletic = True
#   - Mu = [DATE] - HOW IS THIS DATE OBTAINED??

# SITE MODELS
# - GTR

# MANUAL ADDITIONS:
# beast2 correction for ascertainment bias, specifying the number of invariant A, C, G and T sites as 758511 1449901 1444524 758336. 
# <taxa id="TaxonSet.snps_london_351_dna" spec="TaxonSet">
#   <alignment id="snps_london_351_dna" spec="FilteredAlignment" filter="-">
#     <data idref="snps_london_351_dna_original"/>
#     <constantSiteWeights id="IntegerParameter.0" spec="parameter.IntegerParameter" dimension="4" lower="0" upper="0">758511 1449901 1444524 758336</constantSiteWeights>
#   </alignment>
# </taxa>

# MCMC
# Chain Length = 10000000
# tracelog
# - File Name = <accession><clust>.log
# - Log Every = 10000
# screenlog
# - File Name = <accession><clust>.screen
# - Log Every = 10000
# treelog
# - <accession><clust>.trees
# - Log Every = 10000


# POST
# Tracer

# Densitree

# TreeAnnotator


# -----------

remotes::install_github("ropensci/beautier")

library(babette)
library(seqinr)

remove_tail <- function(x, sep = "_", del = 1){
  sapply(strsplit(x, split = sep, fixed = TRUE),
         function(i) paste(head(i, -del), collapse = sep))
}

# DIRECTORIES

setwd("~/Documents/transmission/")

fasta_dir <- "fasta/"
metadata_local_dir <- "metadata/"

# FILES
fasta_file <- paste0(fasta_dir, "THAILAND_TEST.clust_1.dated.fa")
fasta_id_date_df_outfile <- paste0(metadata_local_dir, remove_tail(basename(fasta_file), sep = "."), ".txt")


# READ IN FILES
fasta <- read.fasta(file = fasta_file, forceDNAtolower = F)

# Get fasta sample names and parse into name and date
fasta_names <- names(fasta)

fasta_id_date_df <- do.call(rbind, lapply(strsplit(fasta_names ,"_"), function(x){
  data.frame(id = paste(x[1:(length(x))], collapse = "_"), year = x[length(x)])
  }))

# Save df of ids and year as csv for model
write.table(fasta_id_date_df, file = fasta_id_date_df_outfile, quote = F, row.names = F, col.names = F, sep = "\t")


# CLOCK MODEL
clock_rate <- 0.0000001
clock_model <- create_strict_clock_model(clock_rate_param = create_clock_rate_param(value = clock_rate),
                                         clock_rate_distr = create_log_normal_distr(value = clock_rate, m = 1, s = 1.25))

# MCMC
every <- 10000

mcmc <- create_mcmc(
  chain_length = 1e+08,
  tracelog = beautier::create_tracelog(log_every = every),
  screenlog = beautier::create_screenlog(log_every = every),
  treelog = beautier::create_treelog(log_every = every)
)

# TREE PRIOR
tree_prior <- create_ccp_tree_prior(
  pop_size_distr = create_log_normal_distr(
    m = 1,
    s = 1.25,
    value = 100.0,
    lower = 0.0,
    upper = 200.0
  ))

# MRCA PRIOR
mrca_prior <- create_mrca_prior(
  is_monophyletic = TRUE, 
  mrca_distr = create_laplace_distr(mu = 1990))


# MAKE XML
create_beast2_input_file(
  fasta_file,
  "beast_xml/THAILAND_TEST.babette.xml",
  site_model = create_gtr_site_model(),
  clock_model = clock_model,
  tree_prior = tree_prior,
  mrca_prior = mrca_prior,
  mcmc = mcmc,
  tipdates_filename = fasta_id_date_df_outfile
)







