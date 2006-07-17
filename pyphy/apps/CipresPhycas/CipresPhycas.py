from Phycas import Phycas
from Likelihood import PhycasIDLishMatrix
import ProbDist
import os, sys
import re
import threading
from PIPRes.util.server_import import *
from PIPRes.wrap.external_program import *
#from PIPRes.cipres_types import * 
from PIPRes.util.cipres import getCipresProperty
from PIPRes.util.io import cipresGetLogger
_LOG = cipresGetLogger('pipres.service_impl.cipres_phycas.py')
from cStringIO import StringIO

class CipresPhycas(CipresIDL_api1__POA.AsyncTreeInfer, SimpleServer, Phycas):
    def __init__(self, registry):
        _LOG.info('Initializing Phycas')
        Phycas.__init__(self)
        self.data_source = 'memory'
        #if sys.platform == 'darwin':
        #    self.data_file_name = r'/Users/mholder/Documents/projects/phycas_svn/phycasdev/pyphy/pyphy/nyldna4.nex'
        #else:
        #    self.data_file_name = r'C:\Synchronized\Projects\pyphydev\pyphy\pyphy\nyldna4.nex'
        SimpleServer.__init__(self, registry)
        self.registry = registry
        self.trees = []
        self.treesLock = threading.Lock()
        self.curr_tree_pos = None
        self.phycasIDLishMatrix = None
        self.ntax = 0
        self.leafSet = []
        self.serverRefCount = 1
        self.settings = {}
        #self.tree = CipresTree('(1,2,(3,4))')
        #self.nTrees = 10

    def remove(self):
        self.serverRefCount -= 1
        _LOG.warn('not removing')
        if self.serverRefCount is None: #< 1:
            SimpleServer.remove(self)

    def toCipresIDLTree(self, t):
        cycle, lnL, newick = t
        name = str(cycle)
        return CipresIDL_api1.Tree(
            newick,
            CipresIDL_api1.TreeScore(CipresIDL_api1.DOUBLE_SCORE_TYPE, lnL),
            self.leafSet,
            name)

    def returnTree(self, ind):
        self.curr_tree_pos = ind
        self.treesLock.acquire()
        try:
            t = self.trees[ind]
        finally:
            self.treesLock.release()
        return self.toCipresIDLTree(t)

    # Scriptable interface    
    def execute(self, cmd):
        _LOG.debug('CipresPhycas.execute cmd=%s' % cmd)
        parts = cmd.split(';')
        for s in parts:
            if s:
                splitOnEqual = s.split('=')
                if len(splitOnEqual) == 1:
                    return False, "No '=' character seen.  Expecting command in the form of <name> = <value>"
                name = splitOnEqual[0].strip()
                value = '='.join(splitOnEqual[1:]).strip()
                _LOG.debug('splitOnEqual[0]=%s, splitOnEqual[1]=%s, name=%s, value=%s' % (splitOnEqual[0], splitOnEqual[1], name, value))
                self.settings[name] = value
        return True, ''

    def applySettings(self):
        if self.settings.has_key('ncycles'):
            self.ncycles = int(self.settings['ncycles'])
            _LOG.debug('setting self.ncycles = %d' % self.ncycles)
        if self.settings.has_key('sample_every'):
            self.sample_every = int(self.settings['sample_every'])
            _LOG.debug('setting self.sample_every = %d' % self.sample_every)
        if self.settings.has_key('ls_move_weight'):
            self.ls_move_weight = int(self.settings['ls_move_weight'])
            _LOG.debug('setting self.ls_move_weight = %d' % self.ls_move_weight)
        if self.settings.has_key('random_seed'):
            self.random_seed = int(self.settings['random_seed'])
            _LOG.debug('setting self.random_seed = %d' % self.random_seed)
        if self.settings.has_key('using_hyperprior'):
            if int(self.settings['using_hyperprior']):
                self.using_hyperprior = True
                _LOG.debug('setting self.using_hyperprior = True')
                inversegamma_shape = 2.1
                if self.settings.has_key('inversegamma_shape'):
                    inversegamma_shape = float(self.settings['inversegamma_shape'])
                inversegamma_scale = 1.0/1.1
                if self.settings.has_key('inversegamma_scale'):
                    inversegamma_scale = float(self.settings['inversegamma_scale'])
                self.edgelen_hyperprior = ProbDist.InverseGammaDist(inversegamma_shape, inversegamma_scale)
            else:
                self.using_hyperprior = False
                _LOG.debug('setting self.using_hyperprior = False')

    # AsyncTreeInfer (will be called before AsyncTreeIterator functions)
    def setMatrix(self, mat):
        _LOG.debug('CipresPhycas.setMatrix')
        self.ntax = len(mat.m_matrix)
        if self.ntax == 0:
            self.nchar = 0
            self.matrix = None
            return
        self.nchar = len(mat.m_matrix[0])
        # TEMP hard coding datatype (need to look up how to deal with IDL enums).
        self.phycasIDLishMatrix = PhycasIDLishMatrix(self.ntax, self.nchar, mat.m_symbols, 0, mat.m_numStates)
        for n, stateList in enumerate(mat.m_charStateLookup):
            self.phycasIDLishMatrix.setCharStateLookup(n, stateList)
        for n, row in enumerate(mat.m_matrix):
            self.phycasIDLishMatrix.replaceRow(n, row)

    def inferTrees(self, outStream):
        _LOG.debug('CipresPhycas.inferTrees')
        if self.phycasIDLishMatrix is None:
            raise RuntimeError, 'Failed to call setMatrix before inferTree'
        self.applySettings()
        _LOG.debug('about to call setup')
        self.leafSet = range(1, self.ntax + 1)
        self.setup()
        self.runThread = threading.Thread(target=self.run, name='MCMC')
        _LOG.debug('about to launch mcmc thread')
        self.runThread.start()
        _LOG.debug('back from launching thread')
        self.serverRefCount += 1
        return self

    def copyDataMatrix(self):
        self.likelihood.copyDataFromIDLMatrix(self.phycasIDLishMatrix)
    # AsyncTreeIterator

    def getFirstTree(self):
        _LOG.debug('CipresPhycas.getFirstTree')
        if self.getNumTrees() > 0:
            return self.returnTree(0)
        raise RuntimeError, 'getFirstTree failed because there are currently no trees to supply'

    def getLastTree(self):
        _LOG.debug('CipresPhycas.getLastTree')
        n = self.getNumTrees()
        if n > 0:
            return self.returnTree(n-1)
        raise RuntimeError, 'getLastTree failed because there are currently no trees to supply'

    def getNextTree(self):
        _LOG.debug('CipresPhycas.getNextTree')
        n = self.getNumTrees()
        if (self.curr_tree_pos is not None) and self.curr_tree_pos + 1 < n:
            return self.returnTree(self.curr_tree_pos + 1)
        if self.curr_tree_pos is not None:
            errmsg = 'getNextTree failed because there are currently no trees beyond the current tree (%d)' % self.curr_tree_pos
        else:
            errmsg = 'getNextTree failed because there are currently no trees to supply'
        raise RuntimeError, errmsg
    def getPreviousTree(self):
        _LOG.debug('CipresPhycas.getPreviousTree')
        if (self.curr_tree_pos is not None) and (self.curr_tree_pos > 0):
            return self.returnTree(self.curr_tree_pos - 1)
        if self.curr_tree_pos is not None:
            errmsg = 'getPreviousTree failed because the current tree is the first tree'
        else:
            errmsg = 'getPreviousTree failed because there are currently no trees to supply'
        raise RuntimeError, errmsg
    def getTreeByIndex(self, whichTree):
        _LOG.debug('CipresPhycas.getTreeByIndex(%d)' % whichTree)
        n = self.getNumTrees()
        if (self.curr_tree_pos is not None) and (whichTree >= 0) and (whichTree < n):
            return self.returnTree(whichTree)
        if self.curr_tree_pos is None:
            errmsg = 'getTreeByIndex failed because there are currently no trees to supply'
        elif whichTree >= n:
            errmsg = 'getTreeByIndex failed because the requested tree index (%d) is greater than the current number of trees (%d)' % (whichTree, n)
        else:
            errmsg = 'getTreeByIndex failed because the requested tree index (%d) was less than 0' % whichTree
        raise RuntimeError, errmsg
    def pop(self):
        _LOG.debug('CipresPhycas.pop')
        n = self.getNumTrees()
        if self.curr_tree_pos is not None:
            t = self.returnTrees(self.curr_tree_pos)
            self.trees.pop()
            if n > 0:
                self.curr_tree_pos = n - 1
            else:
                self.curr_tree_pos = None
            return t
        raise RuntimeError, 'pop failed because there are currently no trees to supply'
    def push(self, treeParam):
        raise RuntimeError, 'Phycas does not accept trees via the push method!'
    def shift(self):
        _LOG.debug('CipresPhycas.shift')
        n = self.getNumTrees()
        if self.curr_tree_pos is not None:
            t = self.returnTrees(0)
            self.trees.pop(0)
            if n > 0:
                self.curr_tree_pos = 0
            else:
                self.curr_tree_pos = None
            return t
        raise RuntimeError, 'shift failed because there are currently no trees to supply'
    def unshift(self, treeParam):
        raise RuntimeError, 'Phycas does not accept trees via the unshift method!'
    def getNumTrees(self):
        _LOG.debug('CipresPhycas.getNumTrees')
        self.treesLock.acquire()
        try:
            n = len(self.trees)
        finally:
            self.treesLock.release()
        return n
    def exists(self, treeInd):
        _LOG.debug('CipresPhycas.exists %d ' % treeInd)
        return treeInd < self.getNumTrees()
    # Phycas overrides
    def recordSample(self, cycle, lnL = 0.0):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls base class version, then stores tree description in list so that
        trees can be supplied when requested.
        """
        Phycas.recordSample(self, cycle, lnL)
        
        t = (cycle + 1, lnL, self.tree.makeNewick())
        _LOG.debug(str(t))
        self.treesLock.acquire()
        try:
            self.trees.append(t)
            if not self.curr_tree_pos:
                self.curr_tree_pos = 0
        finally:
            self.treesLock.release()

if __name__=='__main__':
    argv = sys.argv
    try:
        interfaceDict = {}
        for k in ['AsyncTreeInfer']:
            interfaceDict[k] = CipresPhycas
        cipresServe(argv, interfaceDict)
    except:
        logException(_LOG)

