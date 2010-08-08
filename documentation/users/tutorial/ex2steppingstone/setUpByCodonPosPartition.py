################################################################################
# The code hear is a model-configuration portion of a full Phycas analysis.
#
# It assumes that the data have been read and the character subsets first_rag1,
#   second_rag1, third_rag1, first_rag2, second_rag2, and third_rag2 have been
#   defined.
#
# It sets up a model with 3 subsets in the partition, and sets the MCMC log file
#   names to convey the fact that we have partitioned by codon position.
############

################################################################################
# First we create the 3 models. 
############
m1 = model()
m2 = model() 
m3 = model()


################################################################################
# Now we partition the data using a function that takes 3 arguments:
#   1. the character subset
#   2. the model
#   3. A name for subset as a python string (hence we'll use quotes).
############
partition.addSubset(first_rag1 + first_rag2, m1, "First codon positions")
partition.addSubset(second_rag1 + second_rag2, m2, "Second codon positions")
partition.addSubset(third_rag1 + third_rag2, m3, "Third codon positions")
partition()


################################################################################
# Now we'll set up a tag that will be used in the steppingstoneDoMCMC.py script
#   to identify the output of this analysis.
############
analysis_tag = 'by_codon_pos'
