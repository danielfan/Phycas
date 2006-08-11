#!/usr/bin/python
import sys, os
from nexus_public_blocks import *
from tree import *
usingBBEdit = True # True to use bbedit to highlight the lines of Nexus Errors

dirPrefix = 'boot'
	
def getInferredTree():
	''' Reads in the tree inferred from the non-bootstrapped data.'''
	inferredTree = getTreesFromNexusFileName(os.path.join('orig', 'unweighted', 'all.tre'))[0]
	inferredTree.unroot()
	return inferredTree

def countBootDirs():
	global dirPrefix
	nBootReps = 0
	nextDir = dirPrefix + str(nBootReps)
	while os.path.exists(nextDir) and os.path.isdir(nextDir):
		nBootReps += 1
		nextDir = dirPrefix + str(nBootReps)
	if nBootReps == 0:
		sys.exit('No boot# directories found.')
	return nBootReps
	
def classifyBootReps():
	'''Creates a directory for every split in the inferred tree. 
creates boot# subdirectories in each of these split directories for every 
bootstrap replicate that does not have that split under the boot weights.
analyses of the line search for the boundary point occurs in each of these 
directories.  Finally the boundary characters weights are written in each boot 
subdir.'''
	global dirPrefix
	nBootReps = countBootDirs()
	inferredTree = getInferredTree()
	taxMask = inferredTree.taxonMask
	inferredSplits = inferredTree.getSplits(incTrivials = False, orient = -1)
	conflictingRepList = [[] for i in xrange(len(inferredSplits))]
	for n in range(nBootReps):
		p = os.path.join(dirPrefix + str(n), 'bootWts', 'strictCon.tre')	
		bootConTree = getTreesFromNexusFileName(p)[0]
		bootConTree.unroot()
		for splitN, s in enumerate(inferredSplits):
			if not bootConTree.hasSplit(s):
				conflictingRepList[splitN].append(n)
	print conflictingRepList

def writeLScoreCommands():
	'''Adds command files for lscores of all inferred trees on the linear weighting schemes associated with each boot rep.'''
	nBootReps = countBootDirs()
	bootTrees = {}
	for n in range(nBootReps):
		bootNDir = dirPrefix + str(n)
		p = os.path.join(bootNDir, 'bootWts', 'all.tre')	
		bootConTree = getTreesFromNexusFileName(p)
		for t in bootConTree:
			t.unroot()
			bootTrees[t] = -1
		linSearchN = 0
		while True:
			p = os.path.join(bootNDir, 'linSearch' + str(linSearchN), 'all.tre')
			if not os.path.exists(p):
				break
			lsTre = getTreesFromNexusFileName(p)
			for newTree in lsTre:
				treeIsNew = True
				newTree.unroot()
				if not newTree in bootTrees:
					bootTrees[newTree] = linSearchN
		for k, v in bootTrees.iteritems():
			print v, k
		
class NoCorrespondingKeyError(ValueError): pass
if __name__ == '__main__':
	if len(sys.argv) != 2:
		sys.exit('Usage: ' + sys.argv[0] + ' <data set name>')
	dataSetName = sys.argv[1]
	if not os.path.exists(dataSetName) or not os.path.isdir(dataSetName):
		sys.exit('The directory ' + dataSetName + ' does not exist.')
	os.chdir(dataSetName)
	try:
		writeLScoreCommands()
		sys.exit(0)
		classifyBootReps()
	except NexusError, x:
		import os
		print >>sys.stderr,  x
		if usingBBEdit and NexusParsing.currentFile != '':
			if x.startPos != None:
				s = 'bbedit +%d %s' % (x.startPos, NexusParsing.currentFile)
				#s = 'bbedit +%d %s' % (x.startPos.line, NexusParsing.currentFile)
			else:
				s = 'bbedit %s' % (NexusParsing.currentFile)
			os.system(s)
			sys.exit(1)
		NexusParsing.statusStream = None
	sys.exit(0)
