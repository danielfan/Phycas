#!/usr/bin/python
'''
Supplies commonly used classes (utility and exceptions) and functions used through
nexus parsing code (may be imported as *).
'''
import copy
class NexusTaxaManager:
	def __init__(self, taxLabels = None, taxSets = None):
		if taxLabels:
			self.taxLabels = taxLabels
		if taxSets:
			self.taxSets = taxSets
	def translateTaxLabel(self, tLabel):
		'''Checks for label in taxManager field, thentaxLabels, then taxSets (if those fields exist)'''
		if 'taxManager' in self.__dict__:
			return self.taxManager.translateTaxLabel(tLabel)
		knownLabels = self.__dict__.has_key('taxLabels') and self.taxLabels or []
		try:
			return findNexusIndex(tLabel, knownLabels, 'taxon')
		except NexusAfterTokenError: pass
		if self.__dict__.has_key('taxSets') and self.taxSets.has_key(tLabel):
			return self.taxSets[tLabel]
		if self.__dict__.get('supportAnonNumbering') and tLabel.isdigit() and long(tLabel) > 0:
			return long(tLabel) - 1
		raise NexusUnknownTaxonError(tLabel, knownLabels)
	def getNTax(self):
		if 'taxManager' in self.__dict__:
			return self.taxManager.getNTax()
		if 'taxLabels' in self.__dict__:
			return len(self.taxLabels)
		if self.__dict__.get('supportAnonNumbering'):
			return self.__dict__.get('maxAnonIndex', 0)
		return 0
	def getTaxLabel(self, ind):
		if 'taxManager' in self.__dict__:
			return self.taxManager.getTaxLabel(ind)
		if 'taxLabels' in self.__dict__:
			return self.taxLabels[ind]
		if self.__dict__.get('supportAnonNumbering'):
			return str(ind + 1)
		raise IndexError, 'indexing empty NexusTaxaManager'
	def addTaxa(self, newLabels):
		if 'taxManager' in self.__dict__:
			return self.taxManager.addTaxa(newLabels)
		if 'taxLabels' in self.__dict__:
			self.taxLabels.extend(newLabels)
		else:
			self.taxLabels = copy.copy(newLabels)
	def getTaxLabels(self):
		if 'taxManager' in self.__dict__:
			return self.taxManager.getTaxLabels
		if 'taxLabels' in self.__dict__:
			return self.taxLabels
		return []
	def __eq__(self, other):
		if other == None: return False
		if id(self) == id(self): return True
		if 'taxManager' in self.__dict__:
			return self.taxManager == other
		if 'taxManager' in other.__dict__:
			return self == other.taxManager
		return False
		
class NexusError(ValueError):
	def __init__(self, s = None, e = None, m = ''):
		self.startPos = s
		self.endPos = e
		self.message = m
	def __str__(self):
		if self.startPos == None:
			posInfo = 'position unknown'
		else:
			if self.endPos == None or self.endPos == self.startPos:
				s = self.endPos == None and 'starting ' or ''
				posInfo =  '%sat line %s' % (s, self.startPos)
			else:
				posInfo = 'from line %s to line %s' % (self.startPos, self.endPos)
		return 'Nexus Error:  %s (%s)' % (self.message, posInfo)
class NexusAfterTokenError(NexusError):
	def __init__(self, tok, m):
		if tok.__class__.__name__ != 'str':
			sP = tok.__dict__.has_key('startPos') and tok.startPos or None
			eP = tok.__dict__.has_key('endPos') and tok.endPos or None
		else:
			sP, eP = None, None
		NexusError.__init__(self, sP, eP, m)
		
class NexusOpenCommentError(NexusError):
	def __init__(self, startPos):
		NexusError.__init__(self, startPos, None, 'Unterminated comment')
class NexusOpenQuoteError(NexusError):
	def __init__(self, startPos):
		NexusError.__init__(self, startPos, None, 'Unterminated quoted string')
class NexusOpenCommandError(NexusError):
	def __init__(self, name, startPos):
		NexusError.__init__(self, startPos, None, 'Expecting ; to end the %s command' % name)
class NexusOpenBlockError(NexusError):
	def __init__(self, name, startPos):
		NexusError.__init__(self, startPos, None, 'Expecting "END;" to end the %s block' % name)
class NexusBareCommandError(NexusError):
	def __init__(self, c):
		NexusError.__init__(self, c.startPos, c.endPos, 'Expecting BEGIN command to start a new block (found %s)' % c.name)
		self.command = c
class NexusUnexpectedTokenError(NexusError):
	def __init__(self, c, tok):
		NexusError.__init__(self, tok.startPos, tok.endPos, 'Unexpected "%s" in %s' % (str(tok), c))
class NexusMissingTokenError(NexusError):
	def __init__(self, expected, tok):
		NexusError.__init__(self, tok.startPos, tok.endPos, 'Expecting %s but found %s' % (expected, str(tok)))
class NexusUnsupportedError(NexusError):
	def __init__(self, n, tok):
		NexusError.__init__(self, tok.startPos, tok.endPos, '%s is not currently supported' % (n))
class NexusIllegalName(NexusError):
	def __init__(self, n, tok):
		NexusError.__init__(self, tok.startPos, tok.endPos, '%s is not a legal %s label' % (str(tok), n))
class NexusUnknownLabelError(ValueError):
	def __init__(self, lab, type = '', knownLabels = None):
		m = 'Unknown %s label "%s"' % (type, lab)
		if len(knownLabels or []) > 0:
			m = m + ' (possible names = %s)' % str(knownLabels)
		ValueError.__init__(self, m)
class NexusUnknownTaxonError(NexusUnknownLabelError):
	def __init__(self, lab, knownLabels = None):
		NexusUnknownLabelError.__init__(self, lab, 'Taxon', knownLabels)
class NexusBadTreeDefError(NexusAfterTokenError):
	def __init__(self, tok, m):
		NexusAfterTokenError.__init__(self, tok, 'Invalid tree definition:  ' + m)

def notEqualOrNexusError(obj, s1, s2, pref, tok):
	if obj.get(s1) != None:
		if obj[s1] == obj.get(s2): raise NexusAfterTokenError(tok, '%s %s cannot equal %s' %(pref, s1, s2))

firstOnBit15Tuple = (0, 1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1)
def getFirstOnBitAsNumber(a):
	'''Return first bit set to 1 (0->0, 1->1, 2->2, 3->1, 4->3).  horrendous algorithm'''
	global firstOnBit15Tuple
	c = long(a)
	if c != a:
		raise ValueError, 'non-integer arg to getFirstOnBitAsNumber'
	if c < 1L:
		if c == 0L:
			return 0
		raise ValueError, 'negative arg to getFirstOnBitAsNumber' 
	fBit = 0
	while True:
		if c & 0xFL:
			return fBit + firstOnBit15Tuple[c & 0xFL]
		fBit += 4
		c >>= 4
def getFirstBitOnly(a):
	'''returns a long with only the first bit set to 1. (0->0, 1->1, 2->2, 3->1, 4->4, 5->1)
	equivalent to pow(2, getFirstOnBitAsNumber(a) - 1) for a > 0'''
	c = long(a)
	if c != a:
		raise ValueError, 'non-integer arg to getFirstOnBitAsNumber'
	if c < 1L:
		if c == 0L:
			return 0
		raise ValueError, 'negative arg to getFirstOnBitAsNumber'
	if c&1:
		return 1
	return c&(c^(c-1))
def index(item, searchList, op):
	i = 0
	for el in searchList:
		if op(item, el):
			return i
		i += 1
	raise ValueError, 'index(item, searchList, op): item is not in searchList'

def findNexusIndex(labelTok, labels, type = ''):
	searchName = str(labelTok).upper()
	try:
		return index(searchName, labels, lambda c, u: c == u.upper())
	except ValueError:
		if searchName.isdigit():
			i = int(searchName)
			if i <= len(labels):
				return i - 1
		raise NexusAfterTokenError(labelTok, '%s is not a known %s label' % (labelTok, type))


def stableUnique(seq, idfun = None):
	''' returns list with unique elements from seq.
	from http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52560
	and when you can't lose order..., Alex Martelli, 2001/10/13
	unique() systematically loses order, but if items are hashable it's not hard to keep order 
	intact -- only eliminate "later" insertions of already-present items. 
	In case items _aren't_ hashable, for a big enough N using their cPickle.dumps() might be 
	worth it... that's easily generalized to "uniqueness within equivalence classes" -- parameter 
	function idfun must return hashable objects which are == for and only for items that are duplicates.'''
	if idfun is None:
		def idfun(x): return x
	seen = Set()
	result = []
	for item in seq:
		marker = idfun(item)
		if not marker in seen:
			seen.add(marker)
			result.append(item)
	return result

if __name__ == '__main__':
	import doctest, sys
	doctest.testmod(sys.modules[__name__])
