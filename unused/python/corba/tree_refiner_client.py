#!/usr/bin/python
import sys
from omniORB import CORBA
import CosNaming, CipresIDL
import re
from corba_util import *

def makeCorbaCharacters(matRow, symbols):
	return [symbols.index(i) for i in matRow]

def makeCorbaRawMatrix(mat, symbols):
	return [makeCorbaCharacters(row, symbols) for row in mat]

def makeCorbaDataMatrix(mat, symbols = 'ACGT'):
	nChars = len(mat)> 0 and len(mat[0]) or 0
	rawMat = makeCorbaRawMatrix(mat, symbols)
	print rawMat;
	return CipresIDL.DataMatrix(symbols, 4, 1, nChars , rawMat)



nameContext = initOrbAndGetNameService(sys.argv)
treeRefiner = nameContext.getServiceOrExit(CipresIDL.TreeRefine)

from cipres_tree import CipresTree 
cipresTree = CipresTree('(1:0.089964,2:0.050264,3:0.034803,4:0.055375)') 
print 'startingTree = ',  cipresTree
corbaTree = cipresTree.makeCorbaTree()

goodMatrix = [	'ACTT', 
				'ACTT', 
				'ACAT', 
				'ACAT',
				]
badMatrix = [	'ACGT', 
				'ACAT', 
				'ACGT', 
				'ACAT',
				]
rMat = makeCorbaDataMatrix(badMatrix)
print 'about to set matrix'
treeRefiner.setMatrix(rMat)

print 'about to set tree', corbaTree.m_newick
treeRefiner.setTree(corbaTree)
print 'about to refineTree'
for i in range(1000):
	returnedTree = treeRefiner.refineTree()
	print 'returnedTree = ',  returnedTree.m_newick

