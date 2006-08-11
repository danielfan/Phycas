#!/usr/bin/python
from nexus_primitives import *
from nexus_command_reader import *
import re

class NexusTaxLabelsReader:
	taxNamePattern = re.compile(r'.*[a-zA-Z0-9]+.*') 
	def __init__(self, effectiveTaxaBlock, supportUnknownDimension = True):
		self.effTaxaBlock = effectiveTaxaBlock
		self.allowUnknownDimension = supportUnknownDimension
	def validateTaxName(newLabel, prevLabels = []):
		nUpper = str(newLabel).upper()
		if len(filter(lambda x: x== nUpper, prevLabels)) > 0:
			raise NexusAfterTokenError(newLabel, 'the taxon label %s was repeated' % newLabel)
		if nUpper.isdigit():
			if int(nUpper) != (len(prevLabels) + 1):
				raise NexusIllegalName('Taxon', newLabel)
		elif not NexusTaxLabelsReader.taxNamePattern.match(nUpper):
			raise NexusIllegalName('Taxon', newLabel)
		return True
	validateTaxName = staticmethod(validateTaxName)
	def readCommand(self, taxLabelsToken, cStream, obj = None, blockObj = None):
		if blockObj != None: 
			self.effTaxaBlock = blockObj
		self.taxLabels = []
		if self.effTaxaBlock and ('ntax' not in self.effTaxaBlock.__dict__ or self.effTaxaBlock.ntax < 1):
			if not self.allowUnknownDimension:
				raise NexusAfterTokenError(taxLabelsToken, 'The number of taxa must be specified (in a DIMENSIONS command) before taxlabels can be read.')
			nLabels = None
		else:
			nLabels = self.effTaxaBlock.ntax
		tokStream = cStream.getTokenStream()
		lab = tokStream.next()
		while str(lab) != ';':
			NexusTaxLabelsReader.validateTaxName(lab, self.taxLabels)
			if nLabels and len(self.taxLabels) >= nLabels:
				raise NexusAfterTokenError(lab, 'Expecting a ; (Check that the number of taxa specified is correct).')				
			self.taxLabels.append(str(lab))
			lab = tokStream.next()
		if nLabels and len(self.taxLabels) < nLabels:
			raise NexusAfterTokenError(lab, 'Unexpected ; (Check that the number of taxa specified is correct).')				
		if self.effTaxaBlock:
			self.effTaxaBlock.addTaxa(self.taxLabels)
		return True
	
class NexusTaxaBlock(NexusBlock, NexusTaxaManager):
	cmdHandlers = [	NexusCommandReader('Dimensions', [NexusIntSubCommandReader('ntax', 0, NexusSubCommandReader.isPositiveNumber)]), 
					NexusCommandReader('TaxLabels', readerToCreate = NexusTaxLabelsReader) ]
	def endBlockEncountered(self):
		sd = self.__dict__
		if 'taxLabels' not in sd and 'ntax' in sd and 'taxManager' not in sd :
			self.taxLabels = [str(i+1) for i in range(self.ntax)]

def getTaxaFromNexus(inF, getTaxaFromAllPublic = True):
	'''returns first taxa list from the stream.'''
	if getTaxaFromAllPublic:
		import nexus_public_blocks
		handlerDict = nexus_public_blocks.ALL_PUBLIC_BLOCKS
	else:
		handlerDict = {'TAXA': NexusTaxaBlock}
	for b in NexusBlockStream(inF, handlerDict, True, []):
		try:
			tl = b.getTaxLabels()
			if len(tl) > 0:
				return tl 
		except AttributeError: pass
	return []

def getTaxaFromNexusString(s, getTaxaFromAllPublic = True):
	import cStringIO
	return getTaxaFromNexus(cStringIO.StringIO(s), getTaxaFromAllPublic)

def getTaxaFromNexusFileName(inFilename, getTaxaFromAllPublic = True):
	NexusParsing.pushCurrentFile(inFilename)
	inF = open(inFilename)
	tList = getTaxaFromNexusFileName(inF, getTaxaFromAllPublic)
	NexusParsing.popCurrentFile()
	return tList
