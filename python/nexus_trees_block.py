#!/usr/bin/python
from nexus_parser import NexusBlock, NexusBlockStream, NexusParsing
from nexus_primitives import *
from nexus_command_reader import NexusCommandReader
from sets import Set
from tree import Tree
class NexusTreesTranslateReader:
	def __init__(self, treesBlock):
		self.treesBlock = treesBlock
	def readCommand(self, translateToken, cStream, obj = None, blockObj = None):
		if blockObj != None: 
			self.treesBlock = blockObj
		if self.treesBlock.__dict__.has_key('labelTranslation') and len(self.treesBlock.labelTranslation) > 0:
			raise NexusAfterTokenError(translateToken, 'Translation table already constructed:  Translate command (if used) must occur before any trees (and there can only be on translate command).')
		tokStream = cStream.getTokenStream()
		translateDict = {}
		k, v = str(tokStream.next()), str(tokStream.next())
		transKeyOrder = []
		while True:
			if translateDict.has_key(k):
				raise NexusAfterTokenError(k, 'Translate keys cannot occur more than once (%s was duplicated)' % str(k))
			translateDict[k] = v
			transKeyOrder.append(k)
			k = tokStream.next()
			if str(k) == ';':
				break
			if str(k) != ',':
				raise NexusMissingTokenError(', or ;', k)
			k, v = str(tokStream.next()), str(tokStream.next())
		self.treesBlock.addTranslateTable(translateDict, transKeyOrder)
		return True
class NexusTreeReader:
	import re
	treeNamePattern = re.compile(r'.*[a-zA-Z0-9]+.*') 
	def __init__(self, treesBlock):
		self.treesBlock = treesBlock
	def isValidTreeName(s): return NexusTreeReader.treeNamePattern.match(s)
	isValidTreeName = staticmethod(isValidTreeName)
	def readCommand(self, treeCommandNameToken, cStream, obj = None, blockObj = None):
		if blockObj != None: 
			self.treesBlock = blockObj
		if self.treesBlock:
			self.treesBlock.chooseTaxaManager()
		tokStream = cStream.getTokenStream()
		n = tokStream.next()
		treeName = str(n)
		if not NexusTreeReader.isValidTreeName(treeName):
			raise NexusIllegalName('Tree', n)
		tokStream.demandToken('=')
		t = Tree(name = treeName, taxManager = self.treesBlock)
		t.readTokenStream(tokStream)
		self.treesBlock.addTree(t)
		return True
class NexusTreesBlock(NexusBlock, NexusTaxaManager):
	treatUknownNumbersAsLabels = False
	treatUknownNumbersAsIndices = False
	maxAnonymousNumericLabelAllowed = 1000L
	cmdHandlers = [	NexusCommandReader('Translate', readerToCreate = NexusTreesTranslateReader), 
					NexusCommandReader('Tree', readerToCreate = NexusTreeReader) ]
	def __init__(self, beginCmd = None, commandStream = None, previousBlocks = None):
		self.labelTranslation = {}
		if commandStream == None:
			self.prepareToRead(previousBlocks or [])
		self.tolerateNewLabels = True
		self.allowDifferingLeafSets = False
		NexusBlock.__init__(self, beginCmd, commandStream, previousBlocks or [])
	def chooseTaxaManager(self):
		''' called right before a tree is read.  Hook to set the taxa block (if one is unambiguous'''
		if len(self.labelTranslation) > 0 or ('taxManager' in self.__dict__ and self.taxManager): 
			return
		if len(self.potentialTaxaBlocks) > 0:
			firstTaxaBlock = self.potentialTaxaBlocks[0]
			for p in self.potentialTaxaBlocks[1:]:
				if firstTaxaBlock != p:
					return # ambiguous taxa block don't use any of the previous blocks
			self.addNewLabels(firstTaxaBlock.getTaxLabels())
			self.taxManager = firstTaxaBlock
	def prepareToRead(self, previousBlocks):
		''' prepare to read a block.'''
		self.potentialTaxaBlocks = []
		for p in previousBlocks:
			try:
				tl = p.getTaxLabels()
				self.potentialTaxaBlocks.append(p)
			except AttributeError: pass
		self.labelTranslation = {}
		self.trees = []
	def getTrees(self):
		return self.trees
	def addNewLabels(self, labels):
		n = len(self.getTaxLabels())
		for i, l in enumerate(labels):
			ind = n + i
			self.labelTranslation[str(ind + 1)] = ind
			self.labelTranslation[l.upper()] = ind
		if self.__dict__.has_key('taxLabels'):
			self.taxLabels.extend(labels)
		else:
			self.taxLabels = labels
	def addTree(self, t):
		if len(self.__dict__.get('newLabels', []))> 0:
			self.addNewLabels(self.newLabels)
			del self.newLabels
		self.trees.append(t)
	def translateTaxLabel(self, tLabel):
		capLabel = str(tLabel).upper()
		if 'labelTranslation' in self.__dict__:
			ind = self.labelTranslation.get(capLabel)
			if ind != None:
				return ind
		try:
			tl = self.getTaxLabels()
			try:
				if tl:
					return findNexusIndex(tLabel, tl)
			except NexusAfterTokenError:
				if not self.allowDifferingLeafSets:
					raise NexusUnknownTaxonError(tLabel, tl)		
		except AttributeError:
			if not self.allowDifferingLeafSets:
				raise NexusUnknownTaxonError(tLabel, tl)
			tl = []
		if self.tolerateNewLabels:
			if not self.__dict__.has_key('newLabels'):
				self.newLabels = []
			else:
				try:
					return index(capLabel, self.newLabels, lambda c, u: c == u.upper())
				except ValueError: pass
			if not capLabel.isdigit() or NexusTreesBlock.treatUknownNumbersAsLabels:
				self.newLabels.append(str(tLabel))
				return len(tl) + len(self.newLabels) - 1
			if NexusTreesBlock.treatUknownNumbersAsIndices:
				intLabel = int(capLabel)
				if intLabel > 0 and intLabel <= NexusTreesBlock.maxAnonymousNumericLabelAllowed:
					intLabel -= 1
					self.maxAnonIndex = max(self.__dict__.get('maxAnonIndex', 0L), intLabel)
					self.supportAnonNumbering = True
					return intLabel
		raise NexusUnknownTaxonError(tLabel, tl)
	
		
	def addTranslateTable(self, translateDict, translateKeyOrder):
		if self.__dict__.has_key('taxLabels'):
			NexusTreesBlock.validateTranslateTable(translateDict, self, translateKeyOrder)
		else:
			for t in self.potentialTaxaBlocks:
				try:
					if NexusTreesBlock.validateTranslateTable(translateDict, t, translateKeyOrder):
						self.taxManager = t
						break
				except NexusUnknownTaxonError: pass
			if not self.__dict__.has_key('taxLabels'):
				self.taxLabels = [translateDict[k] for k in translateKeyOrder]
				if len(Set(self.taxLabels)) != len(self.taxLabels):
					self.taxLabels = stableUnique(self.taxLabels)
		self.createLabelTranslationDict(translateDict, translateKeyOrder)
	def createLabelTranslationDict(self, translateDict, translateKeyOrder):
		self.labelTranslation = {}
		for k in translateKeyOrder:
			v = translateDict[k]
			inds = self.translateTaxLabel(v)
			self.labelTranslation[k.upper()] = inds
			self.labelTranslation[v.upper()] = inds
		for i in range(self.getNTax()):
			self.labelTranslation[str(i + 1)] = i
			lab = self.getTaxLabel(i).upper()
			self.labelTranslation[lab] = i
	def validateTranslateTable(translateDict, taxaContext, translateKeyOrder):
		for k in translateKeyOrder:
			taxaContext.translateTaxLabel(translateDict[k])
		return True
	validateTranslateTable = staticmethod(validateTranslateTable)
	
def getTreesFromNexus(inF, getTaxaFromAllPublic = True):
	tList = []
	if getTaxaFromAllPublic:
		import nexus_public_blocks
		handlerDict = nexus_public_blocks.ALL_PUBLIC_BLOCKS
	else:
		handlerDict = {'TREES': NexusTreesBlock}
	for b in NexusBlockStream(inF, handlerDict, True, []):
		try:
			tList.extend(b.getTrees())
		except AttributeError: pass
	return tList

def getTreesFromNexusString(s, getTaxaFromAllPublic = True):
	import cStringIO
	return getTreesFromNexus(cStringIO.StringIO(s), getTaxaFromAllPublic)

def getTreesFromNexusFileName(inFilename, getTaxaFromAllPublic = True):
	import os
	if not os.path.exists(inFilename):
		m = 'The tree file %s does not exist' % inFilename
		raise ValueError, m
	NexusParsing.pushCurrentFile(inFilename)
	inF = open(inFilename)
	tList = getTreesFromNexus(inF, getTaxaFromAllPublic)
	NexusParsing.popCurrentFile()
	inF.close()
	return tList

