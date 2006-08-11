#!/usr/bin/python
from nexus_primitives import *
from nexus_command_reader import *
from mth_seqs import *
from sets import Set
import nexus_taxa_block
	
class NexusMatrixReader:
	def __init__(self, containingBlock): pass
	def verifyBlockField(blockObj, s, tok, v = None):
		if not blockObj.__dict__.has_key(s):
			#print blockObj.__dict__
			raise NexusAfterTokenError(tok, 'FORMAT command must precede Matrix (%s field missing)' % s)
		if v != None and blockObj.__dict__[s] != v:
			raise NexusAfterTokenError(tok, 'Currently FORMAT %s must equal %s' % (s, str(v)))
		return True
	verifyBlockField = staticmethod(verifyBlockField)
	
	def translateSeqList(sList, cDict, missingCode):
		strForm = []
		for c in sList:
			strForm.append(cDict.get(c, missingCode))
		return ''.join(strForm)
	translateSeqList = staticmethod(translateSeqList)
	
	def readCommand(self, cName, cStream, obj, blockObj):
		NexusMatrixReader.verifyBlockField(blockObj, 'interleave', cName, False)
		NexusMatrixReader.verifyBlockField(blockObj, 'labels', cName)
		NexusMatrixReader.verifyBlockField(blockObj, 'symbolDict', cName)
		NexusMatrixReader.verifyBlockField(blockObj, 'symbols', cName)
		NexusMatrixReader.verifyBlockField(blockObj, 'codeToSymbolDict', cName)
		NexusMatrixReader.verifyBlockField(blockObj, 'nchar', cName)
		NexusMatrixReader.verifyBlockField(blockObj, 'ntax', cName)
		NexusMatrixReader.verifyBlockField(blockObj, 'taxLabels', cName)
		# print 'appear to have all fields'
		self.readNonInterleavedMatrix(cStream.getTokenStream(), obj, blockObj)
		cStream.demandCommandEnd()
		return True
	def readDataCell(self, tokStream, blockObj):
		c = tokStream.nextChar()
		try:
			return blockObj.symbolDict[c]
		except KeyError:
			if c == '(':
				raise NexusAfterTokenError(tokStream.cbcTok(), 'Polymorphic data cells are not supported yet')
			if c == '{':
				n = 0
				try:
					while True:
						c = tokStream.nextChar()
						if c == '}':
							break
						n += blockObj.symbolDict[c]
					return n
				except KeyError:pass
			raise NexusAfterTokenError(tokStream.cbcTok(), 'Unexpected character %s (at' % c)
	def read_n_chars(self, n, tokStream, blockObj):
		nc = []
		for nCharRead in range(n):
			try:
				nc.append(self.readDataCell(tokStream, blockObj))
			except NexusError, n:
				n.message  = n.message + ' (at character number %d)' % (nCharRead+1)
				raise n
		if tokStream.cbcMoreCharacters():
			tok = tokStream.cbcTok()
			raise NexusAfterTokenError(tok, 'Expecting the end of the character array, found %s' % tokStream.nextChar())
		return nc
	def readNonInterleavedMatrix(self, tokStream, obj, blockObj):
		addLabels = len(blockObj.taxLabels) == 0
		
		if addLabels and not blockObj.labels:
			addLabels = False
			blockObj.taxLabels = range(i, blockObj.ntax + 1)
		dataForTaxInd = ['' for i in range(blockObj.ntax)]
		taxIndSet = Set()
		for nTaxRead in range(blockObj.ntax):
			if blockObj.labels:
				nameTok = tokStream.next()
				if addLabels and nexus_taxa_block.NexusTaxLabelsReader.validateTaxName(nameTok, blockObj.taxLabels):
					blockObj.taxLabels.append(str(nameTok))
					taxInd = nTaxRead
				else:
					taxInd = findNexusIndex(nameTok, blockObj.taxLabels, 'taxon')
			if taxInd in taxIndSet:
				raise NexusAfterTokenError(nameTok, 'Data for taxon number %d (name = %s) have already been read' % (taxInd + 1, blockObj.taxLabels[taxInd].upper()))
			taxIndSet.add(taxInd)
			NexusParsing.statusMessage('Reading data for %s...\n' % blockObj.taxLabels[taxInd])
			dataForTaxInd[taxInd] = self.read_n_chars(blockObj.nchar, tokStream, blockObj)
		 	
		if 0:
				#now well create a biomatrix of bioseq object
			if blockObj.dataType != 'DNA' and blockObj.dataType != 'PROTEIN':
				return True
			from Bio import Alphabet
			from Bio.Alphabet import IUPAC
			if blockObj.dataType == 'DNA':
				alfa = Alphabet.Gapped(IUPAC.ambiguous_dna)
				missingIUPAC = 'N'
			elif blockObj.dataType == 'PROTEIN':
				NexusParsing.statusMessage('Converting protein sequence to IUPAC code all sites with any ambiguity will become "X".\n')
				alfa = Alphabet.Gapped(IUPAC.protein)
				missingIUPAC = 'X'
			bmat = []
			from Bio.Seq import Seq
			NexusParsing.statusMessage('Converting to IUPAC code all missing sites will become %s.\n' % missingIUPAC)
			for nTaxRead in range(blockObj.ntax):
				NexusParsing.statusMessage('Translating data for %s to Bio.Seq() ...\n' % blockObj.taxLabels[taxInd])
				s = NexusMatrixReader.translateSeqList(dataForTaxInd[nTaxRead], blockObj.codeToSymbolDict, missingIUPAC)
				bmat.append(Seq(s, alfa))
			obj.bioMatrix = bmat
		if blockObj.dataType == 'DNA':
			obj.matrix = [DNASeq(i) for i in dataForTaxInd]
		elif blockObj.dataType == 'PROTEIN':
			obj.matrix = [AASeq(i) for i in dataForTaxInd]
		else:
			obj.matrix = dataForTaxInd
		NexusParsing.statusMessage('Data matrix read...\n')
		return True
class EquateSubCommandReader(NexusUnsupportedSubCommand):
	def __init__(self):
		super(EquateSubCommandReader, self).__init__('Equate', '')
	
def verifyNotInSymbols(c, sym, name, tok):
	if c != None and sym.count(c) > 0:
		raise NexusAfterTokenError(tok, 'The "%s" character (%s) may not be included in the SYMBOLS list' % (name, c))
	return True
			
class NexusCharactersBlock(NexusBlock, NexusTaxaManager):
	cmdHandlers = []
	illegalSpecialChars = '''()[]{}/\\,;:=*'"`<>^'''
	def isValidGapChar(subCmd, c, tok, obj):
		if len(c) > 1 or len(c.strip()) == 0 or NexusCharactersBlock.illegalSpecialChars.find(c) != -1:
			NexusError(tok.startPos, tok.endPos, '%s is not a valid GAP character' % c)
		return True
	isValidGapChar = staticmethod(isValidGapChar)
	def isValidMatchChar(subCmd, c, tok, obj):
		if len(c) > 1 or len(c.strip()) == 0 or NexusCharactersBlock.illegalSpecialChars.find(c) != -1:
			NexusError(tok.startPos, tok.endPos, '%s is not a valid MATCH character' % c)
		return True
	isValidMatchChar = staticmethod(isValidMatchChar)
	def isValidSymbolsList(subCmd, c, tok, obj):
		s =''.join(c.split()) 
		s = list(Set([i for i in s]))
		s.sort()
		obj['Symbols'] = ''.join(s)
		return True
	isValidSymbolsList = staticmethod(isValidSymbolsList)
	def isValidMissingChar(subCmd, c, tok, obj):
		if len(c) > 1 or len(c.strip()) == 0 or NexusCharactersBlock.illegalSpecialChars.find(c) != -1:
			NexusError(tok.startPos, tok.endPos, '%s is not a valid MISSING character' % c)
		return True
	isValidMissingChar = staticmethod(isValidMissingChar)
	
	def formatCommandIsValid(self, cmd, obj, tok):
		_npb_identityTrans = string.maketrans('', '')
		pref = 'The FORMAT Command\'s'
		notEqualOrNexusError(obj,'MatchChar', 'Gap', pref, tok)
		notEqualOrNexusError(obj,'MatchChar', 'Missing', pref, tok)
		notEqualOrNexusError(obj,'Gap', 'Missing', pref, tok)
		missing = obj.get('Missing')
		matchChar = obj.get('MatchChar')
		gap = obj.get('Gap')
		dt = obj['DataType'].upper()
		self.dataType = dt
		sym = obj['Symbols']
		respectCase = obj['RespectCase']
		assert sym != None
		equates = obj.get('Equate')
		if equates == None: equates = [] # key might exist in dict, but still be at the default None
		if dt == 'CONTINUOUS':
			raise NexusAfterTokenError(tok, 'CONTINUOUS datatype is not supported yet')
		if dt == 'NUCLEOTIDE':
			sym = sym.translate(_npb_identityTrans, 'Uu') # U is dealt with using an equate
			equates.insert(0, ('U', 'T'))							
		if dt == 'STANDARD':
			native = sym == '' and '01' or ''
			nativeEquates = []
		elif dt == 'DNA' or dt == 'Nucleotide':
			native = DNASeq.symbols[:-1] # we don't always use the -
			nativeEquates = DNASeq.equates
		elif dt == 'RNA':
			native = 'ACGU'
			nativeEquates = NexusCharactersBlock.rnaEquates
		elif dt == 'PROTEIN':
			native = AASeq.symbols[:-1] # we don't always add the -
			nativeEquates = AASeq.equates
		nativeEquates.extend(equates)
		equates = nativeEquates
		verifyNotInSymbols(gap, native, 'Gap', tok)
		verifyNotInSymbols(missing, native, 'Missing', tok)
		verifyNotInSymbols(matchChar, native, 'MatchChar', tok)
		if gap != None:
			native = native + gap
		if len(sym) > 0:
			self.dataType = 'STANDARD' # once you introduce new symbols you are in no-man's land in terms of data type
			raise NexusAfterTokenError(tok, 'User defined symbols are not supported yet')
		sym = sym.translate(_npb_identityTrans, native)
		if not respectCase:
			sym = sym.translate(_npb_identityTrans, native.lower())
		sym = native + sym
		self.symbols = sym
		self.gap = gap
		self.missing = missing
		self.matchChar = matchChar
		self.equates = equates
		self.respectCase = respectCase
		self.taxLabels = [] # ['a', 'b', 'c', 'D'] #//@ need to support taxa block/taxLabels command
		constructNexusSymbolsTranslation(self, self.respectCase, self.missing, self.dataType != 'PROTEIN', tok)
		obj['Symbols'] = sym
		if self.matchChar!= None: raise  NexusAfterTokenError(tok, 'MatchChar is not supported yet.')
		return True
	formatCommandIsValid = staticmethod(formatCommandIsValid)
class NexusDataBlock(NexusCharactersBlock):
	def writeNexusBlock(self, out):
		ntax = len(self.taxLabels)
		assert(ntax == self.ntax and ntax == len(self.matrix))
		nchar = self.nchar
		assert(nchar == len(self.matrix[0]))
		out.write('BEGIN DATA;\n\tDimensions ntax = %d nchar = %d;\n\tFormat datatype = %s gap = -;\n\tmatrix\n' % (ntax, nchar, self.dataType))
		tokenizedTaxLabels = [NexusToken.escapeString(i) for i in self.taxLabels]
		maxLabelLen = len(tokenizedTaxLabels[0])
		for i in range(1,ntax):
			maxLabelLen = max(maxLabelLen, len(tokenizedTaxLabels[i]))
		formatStr = '%%-%ds %%s\n' % maxLabelLen
		for i in range(0, ntax):
			out.write(formatStr % (tokenizedTaxLabels[i], self.matrix[i]))
		out.write(';\nEND;\n')
		
def initializeCharDataBlock():
	if len(NexusCharactersBlock.cmdHandlers) == 0:
		NexusCharactersBlock.dimensionsCommand = NexusCommandReader('Dimensions', [
							NexusBoolSubCommandReader('NewTaxa', False),
							NexusIntSubCommandReader('ntax', 0, NexusSubCommandReader.isPositiveNumber),
							NexusIntSubCommandReader('nchar', 0, NexusSubCommandReader.isPositiveNumber)])
		NexusCharactersBlock.formatCommand = NexusCommandReader('Format', [
							NexusChoiceSubCommandReader('DataType', ['Standard', 'DNA', 'RNA', 'Nucleotide', 'Protein', 'Continuous'], 0, None),
							NexusChoiceSubCommandReader('Items', ['Min', 'Max','Median', 'Average', 'Var', 'SampleSize', 'States'], -1, None),
							NexusChoiceSubCommandReader('StateFormat', ['StatesPresent', 'Individuals','Count', 'Frequency'], -1, None),
							NexusBoolSubCommandReader('Tokens', False),
							NexusBoolSubCommandReader('RespectCase', False),
							NexusBoolSubCommandReader('Transpose', False),
							NexusBoolSubCommandReader('interleave', False),
							NexusBoolSubCommandReader('labels', True),
							NexusStringSubCommandReader('Gap', None, NexusCharactersBlock.isValidGapChar),
							NexusStringSubCommandReader('MatchChar', None, NexusCharactersBlock.isValidMatchChar),
							NexusStringSubCommandReader('Symbols', '', NexusCharactersBlock.isValidSymbolsList),
							EquateSubCommandReader(),
							NexusStringSubCommandReader('Missing', '?', NexusCharactersBlock.isValidMissingChar),],
							NexusCharactersBlock.formatCommandIsValid)
							
		NexusCharactersBlock.eliminateCommand = NexusUnsupportedCommand('Eliminate')
		NexusCharactersBlock.taxLabelsCommand = NexusCommandReader('TaxLabels', readerToCreate = nexus_taxa_block.NexusTaxLabelsReader)
		NexusCharactersBlock.charStateLabelsCommand = NexusUnsupportedCommand('CharStateLabels')
		NexusCharactersBlock.charLabelsCommand  = NexusUnsupportedCommand('CharLabels')
		NexusCharactersBlock.stateLabelsCommand  = NexusUnsupportedCommand('StateLabels')
		NexusCharactersBlock.matrixCommand  = NexusCommandReader('Matrix', readerToCreate = NexusMatrixReader)
		NexusCharactersBlock.cmdHandlers = [	NexusCharactersBlock.dimensionsCommand,
												NexusCharactersBlock.formatCommand,
												NexusCharactersBlock.eliminateCommand,
												NexusCharactersBlock.taxLabelsCommand,
												NexusCharactersBlock.charStateLabelsCommand,
												NexusCharactersBlock.charLabelsCommand,
												NexusCharactersBlock.stateLabelsCommand,
												NexusCharactersBlock.matrixCommand
												]
		NexusDataBlock.dimensionsCommand = NexusCommandReader('Dimensions', 
								[NexusIntSubCommandReader('ntax', 0, NexusSubCommandReader.isPositiveNumber),
								NexusIntSubCommandReader('nchar', 0, NexusSubCommandReader.isPositiveNumber)])
		NexusDataBlock.cmdHandlers = [NexusDataBlock.dimensionsCommand]
		NexusDataBlock.cmdHandlers.extend(NexusCharactersBlock.cmdHandlers[1:])
		t_to_u = string.maketrans('tT', 'uU')
		NexusCharactersBlock.rnaEquates = [(i[0].translate(t_to_u, ''), i[1].translate(t_to_u, '')) for i in DNASeq.equates]
initializeCharDataBlock()	
