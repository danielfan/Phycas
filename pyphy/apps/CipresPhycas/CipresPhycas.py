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
        self.nTax = 0
        self.tree = CipresTree('(1,2,(3,4))')
        self.nTrees = 10

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
        return self.tree
    def getLastTree(self):
        _LOG.debug('CipresPhycas.getLastTree')
        return self.tree
    def getNextTree(self):
        _LOG.debug('CipresPhycas.getNextTree')
        return self.tree
    def getPreviousTree(self):
        _LOG.debug('CipresPhycas.getPreviousTree')
        return self.tree
    def getTreeByIndex(self, whichTree):
        _LOG.debug('CipresPhycas.getTreeByIndex(%d)' % whichTree)
        return self.tree
    def pop(self):
        _LOG.debug('CipresPhycas.pop')
        return self.tree
    def push(self, treeParam):
        raise RuntimeError, 'Phycas does not accept trees via the push method!'
    def shift(self):
        _LOG.debug('CipresPhycas.shift')
        return self.tree
    def unshift(self, treeParam):
        raise RuntimeError, 'Phycas does not accept trees via the unshift method!'
    def getNumTrees(self):
        _LOG.debug('CipresPhycas.getNumTrees')
        return self.nTrees
    def exists(self, treeInd):
        _LOG.debug('CipresPhycas.exists %d ' % treeInd)
        return treeInd < self.nTrees



if __name__=='__main__':
    argv = sys.argv
    try:
        interfaceDict = {}
        for k in ['AsyncTreeInfer']:
            interfaceDict[k] = CipresPhycas
        cipresServe(argv, interfaceDict)
    except:
        logException(_LOG)

