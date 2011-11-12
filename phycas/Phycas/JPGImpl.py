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
        self.opts       = opts
        self.detailsf   = None

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

    def detailsFileOpen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the sampling details file and writes a header line.
        
        """
        details_file_spec = self.opts.out.details
        self.detailsf = None
        try:
            self.detailsf = details_file_spec.open(self.stdout)
        except:
            print '*** Attempt to open details file (%s) failed.' % self.opts.out.details.filename
        
        if self.detailsf:
            self.detailsf.write('tree\trep\tscaling factor\tkappa\tlog-likelihood\n')

    def detailsFileClose(self):
        if self.detailsf is not None:
            self.detailsf.close()

    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            
        """
        self._loadData()
        partition.resetPartition()
        partition.validate(self.nchar)
        
        core = LikelihoodCore(self)
        core.setupCore()

        model_specs = partition.getModels()
        if len(model_specs) != 1:
            self.stdout.error("Expecting just one model, but found %d models" % (len(model_specs),))
            return
        mspec = model_specs[0]
        m = core.partition_model.getModel(0)

        kdist = None
        sfdist = None
        if m.__class__.__name__ == Likelihood.IrreversibleModel.__name__:
            m.setScalingFactorPrior(mspec.scaling_factor_prior.cloneAndSetLot(self.opts.rng))
            sfdist = m.getScalingFactorPrior()
        elif m.__class__.__name__ == Likelihood.BinaryModel.__name__:
            m.setScalingFactorPrior(mspec.scaling_factor_prior.cloneAndSetLot(self.opts.rng))
            sfdist = m.getScalingFactorPrior()
            m.setKappaPrior(mspec.kappa_prior.cloneAndSetLot(self.opts.rng))
            kdist = m.getKappaPrior()
        else:
            self.stdout.error("Expecting model to be type 'gain', 'loss' or 'binary', but instead found type %s" % (mspec.type,))
            return

        tr_source = self.opts.tree_source
        if tr_source is None:
            self.stdout.error("Tree source does not appear to have been specified")
            return

        self.detailsFileOpen()
        self.output('Model name: %s' % (m.getModelName(),))
                              
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
                    core.setTree(t)
                    core.prepareForLikelihood()
                    for rep in range(self.opts.nreps):
                        phi = sfdist.sample()
                        m.setScalingFactor(phi)
                        if kdist is not None:
                            k = kdist.sample()
                            m.setKappa(k)
                        lnL = core.calcLnLikelihood()
                        lnLarr.append(lnL)
                        if self.detailsf:
                            if kdist is not None:
                                self.detailsf.write('%d\t%d\t%g\t%g\t%g\n' % (tree_index,rep+1,phi,k,lnL))
                            else:
                                self.detailsf.write('%d\t%d\t%g\t\t%g\n' % (tree_index,rep+1,phi,lnL))
            if len(lnLarr) > 0:
                self.output('log avg. likelihood = %g' % (self._calcAvgLike(lnLarr),))
            else:
                self.stdout.error("Sample array empty, possibly because from and to settings are incorrect for the tree file being processed")
        except:
            self.stdout.error("Analysis could not be completed (sorry, I can't tell you why)")
            raise

        self.detailsFileClose()
        