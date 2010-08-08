################################################################################
# This code is the initial portion of a phycas analysis. It reads in the data.
#
# If you know that there are some settings that you'll use in other steps of
#   every analysis then you can put those instructions here.
#
# In this case, we'll read in the data, define some character subsets, 
############

################################################################################
# We need to tell python to import the functionality of the Phycas package...
################################################################################
from phycas import * 

################################################################################
# Now we'll tell Phycas what data file it should read...
################################################################################
mcmc.data_source = 'murphy29.rag1rag2.nex'


################################################################################
# Now we'll define 6 subsets of the characters...
################################################################################
first_rag1 = subset(1, 774, 3)
second_rag1 = subset(2, 774, 3)
third_rag1 = subset(3, 774, 3)
first_rag2 = subset(777, 1218, 3)
second_rag2 = subset(775, 1218, 3)
third_rag2 = subset(776, 1218, 3)
 


################################################################################
# We are going to use a different GTR+Gamma model.
################################################################################
model.type = "gtr"
model.num_rates = 4
model.gamma_shape = 0.5
model.gamma_shape_prior = Exponential(1.0)
model.update_freqs_separately = False

################################################################################
# Currently, the steppingstone method works whith fixed trees only.  Below
#   is the newick representation that I obtained for this dataset using a very
#   quick-and-dirty ML tree search with PAUP.
#
# Note that the trees are specified in terms of taxon number rather than name!
#
# We can use the print statement (with quoted strings of characters) to write
#   out status messages to remind us of what is happening.
################################################################################
print("Setting  tree that was obtained by a quick paup search...")
mcmc.starting_tree_source = TreeCollection(newick='(1:0.175599,(((((2:0.033518,4:0.038657):0.006866,3:0.041663):0.018831,((((5:0.064793,19:0.057497):0.003357,(24:0.034658,(25:0.019916,26:0.024097):0.007055):0.010581):0.001613,((18:0.053747,(27:0.040266,28:0.030639):0.016796):0.003339,((((20:0.016474,21:0.037134):0.015437,23:0.047395):0.005919,22:0.043220):0.014677,29:0.038711):0.002928):0.006034):0.011142,((((11:0.024795,12:0.039909):0.072141,(16:0.039810,17:0.035348):0.003100):0.003750,(13:0.030470,14:0.055289):0.019128):0.001971,15:0.037751):0.004012):0.005981):0.010440,((7:0.027685,8:0.023827):0.010063,10:0.040980):0.003835):0.017996,9:0.068549):0.024578,6:0.105696)')
