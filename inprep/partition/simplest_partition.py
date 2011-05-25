# Partitioned:
#   PAUP*  -29.056571 = (-2.349281) + (-8.170907) + (-6.466723) + (-3.567874) + (-8.501786)
#                       <------------- -16.986911 ------------>   <------ -12.06966 ------>
#                                          JC                       F81 (0.1 0.2 0.3 0.4)
#  taxon1 AAA AA
#  taxon2 ACA AC
#  taxon3 AGC AG
#  taxon4 ATT AT
#
#  tree only = [&U] (taxon1:0.2, taxon2:0.2, (taxon3:0.2, taxon4:0.2):0.2);
  
from phycas import *

if partitioning:
    print '==> Setting up partition model...'
    
    if False:
        model.type = 'jc'
        model.edgelen_prior = Exponential(10.0)
        model.edgelen_hyperprior = None

    else:
        model.type = 'jc'
        model.edgelen_prior = Exponential(10.0)
        model.edgelen_hyperprior = None
        m1 = model()
        
        model.type = 'hky'
        model.kappa = 1.0
        model.state_freqs = [0.1, 0.2, 0.3, 0.4]
        model.update_freqs_separately = False
        model.edgelen_prior = Exponential(10.0)
        model.edgelen_hyperprior = None
        m2 = model()
        
        partition.addSubset([1,2,3], m1, 'first')
        partition.addSubset([4,5],   m2, 'second')
        partition()
else:
    print '==> Setting up model...'
    
    model.type = 'jc'
    model.edgelen_prior = Exponential(10.0)
    model.edgelen_hyperprior = None

setMasterSeed(13579)
filename = 'simplest.nex'
blob = readFile(filename)

like.data_source = blob.characters
like.tree_source = TreeCollection(newick='(1:0.2,2:0.2,(3:0.2,4:0.2):0.2);') 
like.starting_edgelen_dist = None
like.store_site_likes = True
print '==> calling like() <=='
lnL = like()
print 'log-likelihood =',lnL
