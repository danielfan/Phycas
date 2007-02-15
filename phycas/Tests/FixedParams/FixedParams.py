from phycas import *

phycas = Phycas()
phycas.r.setSeed(13579)

# set up HKY model
phycas.model = Likelihood.HKYModel()

phycas.model.setKappa(4.0)
phycas.model.fixKappa()
phycas.model.setKappaPrior(ProbDist.ExponentialDist(1.0))

phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
phycas.model.fixStateFreqs()
phycas.model.setBaseFreqParamPrior(ProbDist.ExponentialDist(1.0))

phycas.model.setNGammaRates(4)
phycas.model.setShape(0.14)
phycas.model.fixShape()
phycas.model.setDiscreteGammaShapePrior(ProbDist.ExponentialDist(1.0))

phycas.model.setPinvarModel()
phycas.model.setPinvar(0.27)
phycas.model.fixPinvar()
phycas.model.setPinvarPrior(ProbDist.BetaDist(1.0, 1.0))

phycas.model.fixEdgeLengths()
phycas.model.setEdgeLenPrior(ProbDist.ExponentialDist(10.0))

phycas.model.fixEdgeLenHyperprior()
phycas.model.setEdgeLenHyperPrior(ProbDist.InverseGammaDist(2.1, 0.9090909))

# create a likelihood object and prepare the tree
phycas.setupLikelihood()

# read data file
phycas.data_source = 'file'
phycas.readNexusFile('../Data/nyldna4.nex')

# create a random starting tree
phycas.setupTree('random')

# TODO move both of these inside run()
phycas.likelihood.recalcRelativeRates()
phycas.likelihood.prepareForLikelihood(phycas.tree)

# Open new parameter and tree files
# TODO move paramFileOpen and treeFileOpen calls inside run()
phycas.param_file_name = 'params.p'
phycas.paramFileOpen()
phycas.tree_file_name = 'trees.t'
phycas.taxon_labels = phycas.reader.getTaxLabels() 
phycas.treeFileOpen()

phycas.ncycles = 2500
phycas.run()
