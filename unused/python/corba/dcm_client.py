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

from cipres_tree import CipresTree 
cipresTree = CipresTree('(1:0.089964,(2:0.050264,3:0.034803):0.009528,4:0.055375)') 
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
treeImproveObjRef.setMatrix(rMat)
print 'about to set tree', corbaTree.m_newick
treeImproveObjRef.setTree(corbaTree)
print 'about to inferTree'
returnedTree = treeImproveObjRef.inferTree()
print 'returnedTree = ',  returnedTree.m_newick

def getPrunedCopy(corbaTree, leafSet):
	t = Tree(corbaTree.m_newick)
	t.pruneTo(leafSet)
	return t
class RecIDCM():
	def __init__(self, nameContext):
		self.treeImprover = nameContext.getServiceOrExit(CipresIDL.TreeImprove)
		self.treeRefiner = nameContext.getServiceOrExit(CipresIDL.TreeRefine)
		self.treeDecomposer = nameContext.getServiceOrExit(CipresIDL.TreeDecompose)
		self.mergeTrees = nameContext.getServiceOrExit(CipresIDL.TreeMerge)
	def smallEnough(self, tree):
		return True; ########NO RECURSION
	def done(self, tree):
		self.nIters += 1
		return self.nIters > 1;
	def recDCM(self, subTree):
		if self.smallEnough(subTree):
			treeImprover.setTree(subTree)
			return self.treeImprover.improveCurrentTree()
		else:
			leafSets = self.treeDecomposer.leafSetDecompose(subTree)
			tinyTreesList = []
			for leafSet in leafSets:
				getPrunedCopy(subTree, leafSet)
				tinyTreesList.append(self.recDCM(prunedTree))
			tree = self.treeMerger.mergeTrees(tinyTreesList)
			self.treeRefiner.setTree(tree)
			return self.treeRefiner.refineTree()

	def improveTree(self, matrix, tree):
		return self.runDCM(matrix, tree)
	def runDCM(self, matrix, tree):
		#if not tree:
		#	tree = fastTreeInferer.inferTree(matrix)
		self.nIters = 0
		self.treeImprover.setMatrix(matrix)
		self.treeRefiner.setMatrix(matrix)
		while not done(tree):
			tree = self.recDCM(tree)
