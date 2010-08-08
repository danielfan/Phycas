################################################################################
# The code hear is a model-configuration portion of a full Phycas analysis.
#
# It assumes that the data have been read and the character subsets first_rag1,
#   second_rag1, third_rag1, first_rag2, second_rag2, and third_rag2 have been
#   defined.
#
# It sets up a model with 6 subsets in the partition, and sets the MCMC log file
#   names to convey the fact that we have partitioned by gene and codon position.
############

################################################################################
# First we create the 6 models. 
############
m1 = model()
m2 = model() 
m3 = model()
m4 = model()
m5 = model()
m6 = model()
 
################################################################################
# Now we partition the data using a function that takes 3 arguments:
#   1. the character subset
#   2. the model
#   3. A name for subset as a python string (hence we'll use quotes).
############
partition.addSubset(first_rag1, m1, "Rag1 First codon positions")
partition.addSubset(second_rag1, m2, "Rag1 Second codon positions")
partition.addSubset(third_rag1, m3, "Rag1 Third codon positions")
partition.addSubset(first_rag2, m4, "Rag2 First codon positions")
partition.addSubset(second_rag2, m5, "Rag2 Second codon positions")
partition.addSubset(third_rag2, m6, "Rag2 Third codon positions")
partition()


################################################################################
# Now we'll set up a tag that will be used in the steppingstoneDoMCMC.py script
#   to identify the output of this analysis.
############
analysis_tag = 'by_gene_and_codon_pos'

