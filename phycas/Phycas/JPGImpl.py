import os,sys,math,random
from phycas import *
from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from LikelihoodCore import LikelihoodCore

def check(msg = 'check'):
    raw_input(msg)

class JPGImpl(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Needs to be written.
    
    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        
        """
        CommonFunctions.__init__(self, opts)
        self.opts   = opts

    def getStartingTree(self):
        """
            
        """
        return None

    def _loadData(self):
        """
            
        """
        ds = self.opts.data_source
        super_matrix = ds.getMatrix()
        self.phycassert(super_matrix is not None, 'Tried to load data from a non-existant matrix')
        self.data_matrix = super_matrix
        self.taxon_labels = super_matrix.getTaxLabels()
        self.datatype = super_matrix[0].getDatatype()
        self.ntax = self.data_matrix.getNTax()
        self.nchar = self.data_matrix.getNChar()
        self.phycassert(self.ntax > 0, 'Number of taxa in data matrix was 0')
        self.phycassert(len(self.taxon_labels) == self.ntax, "Number of taxon labels does not match number of taxa.")
    
    def _calcAvgLike(self, like_list):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the log of the average likelihood, factoring out the largest
        likelihood in the supplied list to avoid numerical underflow.
            
        """
        maxLnL = max(like_list)
        n = float(len(like_list))
        tmp = [math.exp(x - maxLnL) for x in like_list]
        return maxLnL - math.log(n) + math.log(sum(tmp)) 
        
    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            
        """
        self._loadData()
        
        core = LikelihoodCore(self)
        core.setupCore()

        tr_source = self.opts.tree_source
        if tr_source is not None:
            m = core.partition_model.getModel(0)
            d = m.getScalingFactorPrior()
            tree_index = 0
            lnLarr = []
            try:
                #tr_source.setActiveTaxonLabels(self.taxon_labels)
                it = iter(tr_source)
                while True:
                    try:
                        t = it.next()
                    except StopIteration:
                        break
                    tree_index += 1;
                    if tree_index >= self.opts.fromtree:
                        if tree_index > self.opts.totree:
                            break
                        self.output('tree %d' % tree_index)
                        core.setTree(t)
                        core.prepareForLikelihood()
                        for rep in range(self.opts.nreps):
                            phi = d.sample()
                            self.output('  rep = %g, phi = %g' % (rep,phi))
                            m.setScalingFactor(phi)
                            lnL = core.calcLnLikelihood()
                            lnLarr.append(lnL)
                self.output('log avg. likelihood = %g' % (self._calcAvgLike(lnLarr),))
            except:
                self.stdout.error("Trees could not be obtained from tree_source")
                raise

        