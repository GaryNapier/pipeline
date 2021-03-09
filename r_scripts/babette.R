

library(babette)


# SITE MODEL ----
site_model <- create_site_model_gtr()




# CREATE INFERENCE MODEL

inference_model <- create_test_inference_model(site_model = site_model)

# MCMC ----

# Chain length
# 100,000,000
# 100000000
chain_length <- 100000000

# Store/log/print every.. 
every <- 10000

store_every <- every
log_every <- every

tracelog <- create_tracelog(log_every = every)

screenlog_every <- every

?create_screenlog

treelog_every <- every


create_mcmc(
  chain_length = chain_length,
  store_every = store_every,
  pre_burnin = 0,
  n_init_attempts = 10,
  sample_from_prior = FALSE,
  tracelog = tracelog,
  screenlog = create_screenlog(),
  treelog = create_treelog()
)







# SITE MODELS ----


beast_path <- "/Users/garynapier/miniconda3/bin/beast"

beast_options <- create_beast2_options(beast2_path = beast_path)

fasta_filename <- beautier::get_babette_path("anthus_aco_sub.fas")
library(testthat)
expect_true(file.exists(fasta_filename))

out <- babette::bbt_run_from_model(beast2_options = beast_options, 
                          fasta_filename = fasta_filename)



inference_model <- create_inference_model()

out <- bbt_run_from_model(
  beast2_options = beast_options,
  fasta_filename = fasta_filename,
  inference_model = inference_model
)

# Saves these files to ./
# anthus_aco_sub.log	
# anthus_aco_sub.trees		
# beast2_b1017ffde289.xml.state

# Saves xml file to:
# ~/Library/Caches/beast2_b1017ffde289.xml
# !!!



# --------------------------------------------------------------------------------------------------------------





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


# "Test-driven development"


# -----------


# FILES
fasta_file <- "fasta/THAILAND_TEST.clust_1.dated.fa"


# GENERAL
id <- get_alignment_id(fasta_file, capitalize_first_char_id = FALSE)
check_alignment_id(id)


# CLOCK MODEL
clock_rate <- 0.0000001
clock_model <- create_clock_model_strict(create_clock_rate_param(value = clock_rate, id = NA))

# -------------------------
# -------------------------
# create_clock_rate_param:
# Note
# It cannot be estimated (as a hyper parameter) yet.
# -------------------------
# -------------------------



# PRIORS
# - Tree.t - Coalescent Constant Population
#   - Pop Size = 100 

# - popSize.t
#   - Log Normal
#   - Lower = 0
#   - Upper = 200
#   - Value = 100

create_ccp_tree_prior(
  id = NA,
  pop_size_distr = beautier::create_log_normal_distr(id = NA, m = 1, s = 1.25)
)


?create_param




library(XML)

london_xml_file <- "~/Documents/transmission/beast_xml/BEAST2_GTR_Model_London_TB_Outbreak.xml"

data <- xmlParse(london_xml_file)

xml_data <- xmlToList(data)

str(xml_data)














