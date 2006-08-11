from pyphy.ReadNexus import *
from pyphy.DataMatrix import *
from pyphy.Likelihood import *
import sys, os
from PIPRes.tree import *
import subprocess 
import os, sys
import new

isWindows = sys.platform.upper().startswith('WIN')
	
defaultFileName = 'nyldna4.nex'

def simpleRunProcess(exe, cmdLineArgs, commands):
	invocation  = [exe]
	invocation.extend(cmdLineArgs)
	p = subprocess.Popen( invocation, # cmdLineArgs, executable=exe,
						  stdin = subprocess.PIPE,
						  stdout = subprocess.PIPE,
						  stderr = None,
						  close_fds = True)
	stdout_text, stderr_text = p.communicate(commands)
	return stdout_text.rstrip()

def askokcancel(title,msg, automaticMode= True):
    print title
    if automaticMode:
        print msg
        return True
    else:
        resp = raw_input(msg)
        return bool(resp)

def prepareForLikelihood(self, dataMatrix):
	model = JCModel()
	self.model = model
	treeLikeInfo = allocateTreeLike(model, dataMatrix)
	self.treeLikeInfo = treeLikeInfo
	for nd in iterTips(self.root):
		nd.tipData = allocateTipData(model, dataMatrix, nd.index, treeLikeInfo)
	
	internalIterator = iterPreOrder(self.root, Node.isInternal)
	numRates = 1
	numPatterns = dataMatrix.getNChar()
	numStates = dataMatrix.getNStates()
	for num, nd in enumerate(internalIterator):
		nd.cla = allocateCondLike(model, dataMatrix, treeLikeInfo)
	
def calcLnLikelihood(self):
	internalIterator = iterPostOrder(self.root, Node.isInternal)
	treeLikeInfo = self.treeLikeInfo
	for nd in internalIterator:
		lChild = nd.lChild
		secChild = lChild.rSib
		if lChild.isLeaf():
			calcPMatTranspose(treeLikeInfo, lChild.tipData, lChild.edgeLength)
			if secChild.isLeaf():
				calcPMatTranspose(treeLikeInfo, secChild.tipData, secChild.edgeLength)
				calcCondLikeArrayTwoTips(nd.cla, lChild.tipData, secChild.tipData, treeLikeInfo)
			else:
				calcPMat(treeLikeInfo, secChild.cla, secChild.edgeLength)
				calcCondLikeArrayOneTip(nd.cla, lChild.tipData, secChild.cla, treeLikeInfo)
		else:
			calcPMat(treeLikeInfo, lChild.cla, lChild.edgeLength)
			if secChild.isLeaf():
				calcPMatTranspose(treeLikeInfo, secChild.tipData, secChild.edgeLength)
				calcCondLikeArrayOneTip(nd.cla, secChild.tipData, lChild.cla, treeLikeInfo)
			else:
				calcPMat(treeLikeInfo, secChild.cla, secChild.edgeLength)
				calcCondLikeArrayNoTips(nd.cla, lChild.cla, secChild.cla, treeLikeInfo)
		currNd = secChild.rSib
		while currNd:
			if currNd.isLeaf():
				calcPMatTranspose(treeLikeInfo, currNd.tipData, currNd.edgeLength)
				conditionOnAdditionaTip(nd.cla, currNd.tipData, treeLikeInfo)
			else:
				calcPMat(treeLikeInfo, currNd.cla, currNd.edgeLength)
				conditionOnAdditionaInternal(nd.cla, currNd.cla, treeLikeInfo)
			currNd = currNd.rSib
	return calcTotalLnLikelihood(self.root.cla, treeLikeInfo)

def getScoresFromPaup(fn):
	paupCmd = '''set storebr noWarnRoot;
	execute %s;
	deroot;
	lset nst=1 basefreq=equal rates=equal pinv=0 userbr;
	lscore / ScoreFile=paup.score replace;
	''' % fn
	print simpleRunProcess('paup.OSX', ['-n'], paupCmd)

if __name__ == '__main__':
	Tree.prepareForLikelihood = prepareForLikelihood
	Tree.calcLnLikelihood = calcLnLikelihood
	
	fn = len(sys.argv) < 2 and defaultFileName or sys.argv[1]
	getScoresFromPaup(fn)
	r = NexusReader(-1) # -1 means read all types of blocks
	try:
		r.ReadFile(fn)
	except ValueError, e:
		sys.exit('Error:\n\t%s' % e)
	except int:
		sys.exit('Unknown exception!')
	
	dataMatrix = getDiscreteMatrix(r, 0)
	nTaxa = dataMatrix.getNTax()
	treeList = r.GetTrees()
	for n, treeDesc in enumerate(treeList):
		print 'Tree', n, ':', treeDesc.newick
		tree = Tree(treeDesc.newick, rooted = False, taxaNamingStyle = TaxaNamingEnum.indicesOnly)
		if tree.hasEdgeLengths:
			print str(tree)
			# from time import sleep
			# print 'night, night'
			# sleep(2)
			tree.prepareForLikelihood(dataMatrix)
			lnLike = tree.calcLnLikelihood()
			print 'The log likelihood is', lnLike
		else:
			print 'Tree does not have edge lengths, likelihood will not be calculated.'
	
	os.system('%s paup.score' % (isWindows and 'type' or 'cat'))
	