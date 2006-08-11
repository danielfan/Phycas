#!/usr/bin/python
import string, copy, sys
class MthSeq:
	def findSubSequenceStartPoint(self, potentialSubSet):
		 # TEMP hideous algorithm
		def matchFunc(x, s, p):
			c = s.seqList[x]
			d = p[x]
			inter = c&d
			return inter == c or inter == d # Note we may want "return inter != 0" here and dispense with all of the determining which is a subset and just look for an intersection at every site
		if potentialSubSet._isMatchingSubset(self, matchFunc): return 0
		shorterCopy = copy.copy(self)
		maxOffSet  = len(self.seqList) - len(potentialSubSet)
		if maxOffSet < 1:
			return None
		for i in xrange(1, maxOffSet + 1):
			shorterCopy.pop(0)
			if potentialSubSet._isMatchingSubset(shorterCopy, matchFunc):
				return i
		return None
		
	def isMoreResolvedSubset(self, potentialSuperSet):
		return self._isMatchingSubset(potentialSuperSet, lambda x, s , p : s.seqList[x] & p[x] == s.seqList[x])
	def _isMatchingSubset(self, potentialSuperSet, matchFunc):
		ls = len(self.seqList)
		if ls > len(potentialSuperSet): 
			if ls > len(potentialSuperSet) + 1 or self.seqList[ls-1] != self.stopSeqCode: 
				return False
			ls -= 1
		for i in range(ls):
			if not matchFunc(i, self, potentialSuperSet):
				print >>sys.stderr, self.seqList[i], 'does not match', potentialSuperSet[i], 'in _isMatchingSubset()'
				return i == ls - 1 and potentialSuperSet[i] == self.stopSeqCode
		return True
	def isValidSymbol(self, n):
		return n <= self.missingCode
	def appendCode(self, n):
		if n > self.missingCode:
			s = '%s is not a valid numeric code for %s' %  (n, self.__class__.__name__) 
			raise ValueError
		self.seqList.append(n)
	def insertCode(self, i, n):
		if n > self.missingCode:
			s = '%s is not a valid numeric code for %s' % (n, self.__class__.__name__) 
			raise ValueError, s
		self.seqList.insert(i, n)
	def removeGaps(self):
		gapN = self.symbolDict['-']
		tmp = self.seqList
		self.seqList = []
		for i in tmp:
			if i & gapN:
				if i != gapN:
					self.seqList.append(i - gapN)
			else:
				self.seqList.append(i)
	def pop(self, i): self.seqList.pop(i)
	def __delitem__(self, k): del self.seqList[k]
	def __str__(self): return ''.join([self.translateCode(i) for i in self.seqList])
	def __len__(self): return len(self.seqList)
	def __getitem__(self, k): return self.seqList[k]
	def __contains__(self, potentialSubSet): # TEMP hideous algorithm
		if potentialSubSet.isMoreResolvedSubset(self): return True
		shorterCopy = copy.copy(self)
		maxOffSet  = len(self.seqList) - len(potentialSubSet)
		if maxOffSet < 1: return False
		for i in xrange(maxOffSet):
			shorterCopy.pop(0)
			if potentialSubSet.isMoreResolvedSubset(shorterCopy):
				return True
		return False
	def __iter__(self): return iter(self.seqList)
class DNASeq(MthSeq):
	symbols = 'ACGT-'
	equates = [	('R', '{AG}'),
					('Y', '{CT}'),
					('M', '{AC}'),
					('K', '{GT}'),
					('S', '{CG}'),
					('W', '{AT}'),
					('H', '{ACT}'),
					('B', '{CGT}'),
					('V', '{ACG}'),
					('D', '{AGT}'),
					('N', '{ACGT}'),
					('X', '{ACGT}')]
	listBitLists = [	[],
						[1],
						[2],
						[1, 2],
						[4],
						[1, 4],
						[2, 4],
						[1, 2, 4],
						[8],
						[1, 8],
						[2, 8],
						[1, 2, 8],
						[4, 8],
						[1, 4, 8],
						[2, 4, 8],
						[1, 2, 4, 8]]
	def __init__(self, transSeqList, stringForm = ''):
		self.seqList = [i for i in transSeqList]
	def __copy__(self):
		return DNASeq(copy.copy(self.seqList))
	def translateCode(self, n):
		return DNASeq.codeToSymbolDict.get(n, '?')
	def toCodon(codonString):
		assert(len(codonString) == 3)
		nList = [DNASeq.symbolDict[i] for i in codonString[:3]]
		return DNASeq.numericCodesToCodon(nList)
	toCodon = staticmethod(toCodon)
	def numericCodesToCodon(nList):
		codonCode = 0
		for n in nList[0:3]:
			if n == DNASeq.gapCode: raise ValueError, 'Gap code found in ...ToCodon() method.'
			if n > DNASeq.gapCode:
				n -= DNASeq.gapCode # we'll tolerate possible gaps and just mask off the gap
			codonCode <<= 4
			codonCode += n
		return codonCode
	numericCodesToCodon = staticmethod(numericCodesToCodon)
	def codonCodeToStr(codonNumberCode):
		print hex(codonNumberCode)
		tCode = codonNumberCode & 15
		sCode = (codonNumberCode >> 4) & 15
		fCode = (codonNumberCode >> 8) & 15
		return ''.join([DNASeq.codeToSymbolDict[fCode], DNASeq.codeToSymbolDict[sCode], DNASeq.codeToSymbolDict[tCode]])
	codonCodeToStr = staticmethod(codonCodeToStr)
	def iter_next_nucleotide(self):
		'''skips gaps and masks off the gap code of positions that have the gap bit'''
		for i in self.seqList:
			if i >= DNASeq.gapCode:
				if i != DNASeq.gapCode:
					yield i - DNASeq.gapCode
			else:
				yield i
	def iter_codons(self):
		'''skips gap sites, returns triples. Pads the end with 'N' codes'''
		retList = [0, 0]
		for i in self.iter_next_nucleotide():
			if retList[0] == 0:
				retList[0] = i
			elif retList[1] == 0:
				retList[1] = i
			else:
				yield retList[0], retList[1], i
				retList = [0, 0]
		if retList[0] != 0:
			if retList[1] == 0:
				yield retList[0], DNASeq.anyNonGapCode, DNASeq.anyNonGapCode
			else:
				yield retList[0], retList[1], DNASeq.anyNonGapCode
	def getSubsetCodons(f, s, t):
		all_unambig = []
		for f_unambig in DNASeq.listBitLists[f]:
			for s_unambig in DNASeq.listBitLists[s]:
				for t_unambig in DNASeq.listBitLists[t]:
					all_unambig.append(DNASeq.numericCodesToCodon((f_unambig, s_unambig, t_unambig)))
		return all_unambig
	getSubsetCodons = staticmethod(getSubsetCodons)	
	def iter_ambig_codons(self):
		for f, s, t in self.iter_codons():
			yield DNASeq.getSubsetCodons(f, s, t)
class AASeq(MthSeq):
	symbols = 'ACDEFGHIKLMNPQRSTVWY*-'
	equates = [	('B', '{DN}'),
					('Z', '{EQ}'),
					('X', 'ACDEFGHIKLMNPQRSTVWY*')]
	def __init__(self, nucList):
		self.seqList = nucList
	def __copy__(self):
		return AASeq(copy.copy(self.seqList))
	def translateCode(self, n):
		c = AASeq.codeToSymbolDict.get(n)
		if c == None:
			return n > self.gapCode and '?' or 'X'
		return c		
def addLowerCaseEquateKeys(eList):
	lcEquates = [(i[0].lower(), i[1]) for i in eList]
	eList.extend(lcEquates)

def constructNexusSymbolsTranslation(blockObj, respectCase, missing, useEquatesInOutput, tok = None):
	identityTrans = string.maketrans('', '')
	blockObj.charShift = 1
	blockObj.symbolDict = {}
	blockObj.codeToSymbolDict = {}
	for s in blockObj.symbols:
		if respectCase: 
			blockObj.symbolDict[s] = blockObj.charShift
			blockObj.codeToSymbolDict[blockObj.charShift] = s
		else: # assumes that the symbols list will not have both 'a' and 'A'
			blockObj.symbolDict[s.upper()] = blockObj.charShift
			blockObj.symbolDict[s.lower()] = blockObj.charShift
			blockObj.codeToSymbolDict[blockObj.charShift] = s.upper()
		blockObj.charShift = blockObj.charShift << 1
	blockObj.missingCode = blockObj.charShift - 1
	blockObj.symbolDict[missing] = blockObj.missingCode
	for e, v in blockObj.equates:
		if v.count('(') + v.count(')') > 0:
			raise NexusAfterTokenError(tok, 'Polymorphic Equates are not supported yet')
		eqSyms = v.translate(identityTrans, '{}()')
		n = 0
		try:
			for es in eqSyms:
				n += blockObj.symbolDict[es]
			blockObj.symbolDict[e] = n
			if useEquatesInOutput:
				blockObj.codeToSymbolDict[n] = e.upper() # bio python IUPAC for protein is slightly different than NEXUS (X is {SC}), so we do not use protein equates
		except KeyError, x: 
			raise NexusAfterTokenError(tok, 'Symbol %s in expansion of EQUATE macro "%s => %s" is invalid' % (x, e, v))
	# for k,v in blockObj.symbolDict.iteritems(): 
	#	print k, '= %20s' % hex(v)
		# make sure when we are outputting we use IUPAC in preference to whatever missing character was assigned (if there is a clash)
	if blockObj.symbols.startswith(DNASeq.symbols[:-1]) and blockObj.symbolDict.has_key('N'):
		blockObj.codeToSymbolDict[blockObj.symbolDict['N']] = 'N'
	elif blockObj.symbols.startswith(AASeq.symbols[:-1]) and blockObj.symbolDict.has_key('X'):
		blockObj.codeToSymbolDict[blockObj.symbolDict['X']] = 'X'

addLowerCaseEquateKeys(DNASeq.equates)
constructNexusSymbolsTranslation(DNASeq, False, '?', True)
DNASeq.missingCode = DNASeq.symbolDict['?']
DNASeq.gapCode = DNASeq.symbolDict['-']
DNASeq.anyNonGapCode = DNASeq.symbolDict['N']
DNASeq.stopSeqCode = 0
addLowerCaseEquateKeys(AASeq.equates)
constructNexusSymbolsTranslation(AASeq, False, '?', False)
AASeq.missingCode = AASeq.symbolDict['?']
AASeq.gapCode = AASeq.symbolDict['-']
AASeq.anyNonGapCode = AASeq.symbolDict['X']
AASeq.stopSeqCode = AASeq.symbolDict['*']
