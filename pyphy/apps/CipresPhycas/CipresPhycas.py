from Phycas import Phycas
import os, sys
import re
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
        SimpleServer.__init__(self, registry)
        self.registry = registry
        self.trees = []
        self.curr_tree_pos = None
        #self.nTax = 0
        #self.tree = CipresTree('(1,2,(3,4))')
        #self.nTrees = 10

    # Scriptable interface    
    def execute(self, cmd):
        _LOG.debug('CipresPhycas.execute cmd=%s' % cmd)
        return True, ''

    # AsyncTreeInfer (will be called before AsyncTreeIterator functions)
    def setMatrix(self, mat):
        _LOG.debug('CipresPhycas.setMatrix')
        self.nTax = len(mat.m_matrix)

    def inferTrees(self, outStream):
        _LOG.debug('CipresPhycas.inferTrees')
        if self.nTax == 0:
            raise RuntimeError, 'Failed to call setMatrix before inferTree'
        return self
    
    # AsyncTreeIterator
    def getFirstTree(self):
        _LOG.debug('CipresPhycas.getFirstTree')
        if self.getNumTrees() > 0:
            self.curr_tree_pos = 0
            return self.trees[0]
        else:
            RuntimeError, 'getFirstTree failed because there are currently no trees to supply'
    def getLastTree(self):
        _LOG.debug('CipresPhycas.getLastTree')
        n = self.getNumTrees()
        if n > 0:
            self.curr_tree_pos = n - 1
            return self.trees[self.curr_tree_pos]
        else:
            RuntimeError, 'getLastTree failed because there are currently no trees to supply'
    def getNextTree(self):
        _LOG.debug('CipresPhycas.getNextTree')
        if self.curr_tree_pos and self.curr_tree_pos + 1 < self.getNumTrees():
            self.curr_tree_pos += 1
            return self.trees[self.curr_tree_pos]
        else:
            if self.curr_tree_pos:
                errmsg = 'getNextTree failed because there are currently no trees beyond the current tree (%d)' % self.curr_tree_pos
            else:
                errmsg = 'getNextTree failed because there are currently no trees to supply'
            RuntimeError, errmsg
    def getPreviousTree(self):
        _LOG.debug('CipresPhycas.getPreviousTree')
        if self.curr_tree_pos and (self.curr_tree_pos > 0):
            self.curr_tree_pos -= 1
            return self.trees[self.curr_tree_pos]
        else:
            if self.curr_tree_pos:
                errmsg = 'getPreviousTree failed because the current tree is the first tree'
            else:
                errmsg = 'getPreviousTree failed because there are currently no trees to supply'
            RuntimeError, errmsg
    def getTreeByIndex(self, whichTree):
        _LOG.debug('CipresPhycas.getTreeByIndex(%d)' % whichTree)
        if self.curr_tree_pos and (whichTree >= 0) and (whichTree < self.getNumTrees()):
            self.curr_tree_pos = whichTree
            return self.trees[self.curr_tree_pos]
        else:
            if not self.curr_tree_pos:
                errmsg = 'getTreeByIndex failed because there are currently no trees to supply'
            elif whichTree >= self.getNumTrees():
                errmsg = 'getTreeByIndex failed because the requested tree index (%d) is greater than the current number of trees (%d)' % (whichTree, len(self.trees))
            else:
                errmsg = 'getTreeByIndex failed because the requested tree index (%d) was less than 0' % whichTree
            RuntimeError, errmsg
    def pop(self):
        _LOG.debug('CipresPhycas.pop')
        if self.curr_tree_pos:
            t = self.trees.pop()
            n = self.getNumTrees()
            if n > 0:
                self.curr_tree_pos = n - 1
            else:
                self.curr_tree_pos = None
            return t
        else:
            RuntimeError, 'pop failed because there are currently no trees to supply'
    def push(self, treeParam):
        raise RuntimeError, 'Phycas does not accept trees via the push method!'
    def shift(self):
        _LOG.debug('CipresPhycas.shift')
        if self.curr_tree_pos:
            t = self.trees.pop(0)
            n = self.getNumTrees()
            if n > 0:
                self.curr_tree_pos = 0
            else:
                self.curr_tree_pos = None
            return t
        else:
            RuntimeError, 'shift failed because there are currently no trees to supply'
    def unshift(self, treeParam):
        raise RuntimeError, 'Phycas does not accept trees via the unshift method!'
    def getNumTrees(self):
        _LOG.debug('CipresPhycas.getNumTrees')
        return len(self.trees)
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
        self.trees.append(self.tree.makeNewick())
        if not self.curr_tree_pos:
            self.curr_tree_pos = 0

if __name__=='__main__':
    argv = sys.argv
    try:
        interfaceDict = {}
        for k in ['AsyncTreeInfer']:
            interfaceDict[k] = CipresPhycas
        cipresServe(argv, interfaceDict)
    except:
        logException(_LOG)

