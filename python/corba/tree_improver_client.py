#!/usr/bin/python
import sys
from omniORB import CORBA
import CosNaming, CipresIDL
import re
from corba_util import *



nameContext = initOrbAndGetNameService(sys.argv)
treeImproveObjRef = nameContext.getServiceOrExit(CipresIDL.TreeImprove)

from cipres_tree import CipresTree 
from cipres_matrix import CipresMatrix

cipresTree = CipresTree('(1:0.089964,(2:0.050264,3:0.034803):0.009528,4:0.055375)') 
print 'startingTree = ',  cipresTree
corbaTree = cipresTree.makeCorbaTree()

rMat = CipresMatrix([	'ACGT', 
				'ACAT', 
				'ACGT', 
				'ACAT',
				])
print 'about to set matrix'
treeImproveObjRef.setMatrix(rMat.getCorbaMatrix())
print 'about to set tree', corbaTree.m_newick
treeImproveObjRef.setTree(corbaTree)
print 'about to inferTree'
returnedTree = treeImproveObjRef.inferTree()
print 'returnedTree = ',  returnedTree.m_newick

