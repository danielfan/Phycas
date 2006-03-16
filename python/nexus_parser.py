#!/usr/bin/python
'''

'''
import copy, cStringIO, sys
from nexus_primitives import NexusMissingTokenError, NexusBareCommandError
class NexusParsing:
	statusStream = None
	fileStack = []
	currentFile = ''
	def statusMessage(s):
		if NexusParsing.statusStream != None:
			NexusParsing.statusStream.write(' ...%s' % s)
			NexusParsing.statusStream.flush()
	statusMessage = staticmethod(statusMessage)
	def pushCurrentFile(f):
		if NexusParsing.currentFile != '':
			NexusParsing.fileStack.append(NexusParsing.currentFile)
		currentFile = f
	pushCurrentFile = staticmethod(pushCurrentFile)
	def popCurrentFile():
		NexusParsing.currentFile = len(NexusParsing.fileStack) > 0 and NexusParsing.fileStack.pop(-1) or ''
	popCurrentFile = staticmethod(popCurrentFile)

class PosTriple:
	tabWidth = 4
	def __init__(self, p, l, c):
		self.pos = p
		self.line = l
		self.column = c
	def nextChar(self, c, prev):
		'''returns length of previous line (if c is a newline)'''
		self.pos += 1
		if c == '\n' or c =='\r':
			ret = self.column
			self.column = 0
			if c == '\n' and prev == '\r':
				return -1
			self.line += 1
			return ret
		if c == '\t':
			self.column += PosTriple.tabWidth
		else:
			self.column += 1
		return -1
	def __str__(self): return 'pos=%d line=%d col=%d' % (self.pos, self.line, self.column)

class CharPosStream:
	'''raw class that iterates through a file object returning the character and PosTriple for the characters location'''
	def __init__(self, fileObj):
		#fullPosSupport self.currPos = PosTriple(0, 1, 0)
		self.currPos = 1 
		#fullPosSupport
		self.inputStream = fileObj
		self.prev = ''
	def readNextChar(self):
		c = self.inputStream.read(1)
		if c == '': raise StopIteration
		return c
	def __iter__(self):
		while True:
			c = self.readNextChar()
			#fullPosSupport self.currPos.nextChar(c, self.prev)
			if c == '\n' or c =='\r' and (c == '\r' or prev == '\n'):
				self.currPos += 1
			#fullPosSupport
			self.prev = c
			yield c, self.currPos

class NexusCharStream:
	'''Wraps CharPosStream to allow for:
		- translation of all line-endings (\r, \r\n, and \n) to \n (note: NexusCharStream.peek() may return \r)
		- peek capabilities
	'''
	def __init__(self, fileObj):
		self.charIter = iter(CharPosStream(fileObj))
		#fullPosSupport self.nextCharPos = PosTriple(0, 1, 0)
		self.nextCharPos = 1
		#fullPosSupport 
		
		self._advanceStoredNextChar()
	def _advanceStoredNextChar(self):
		try: 
			self.nextChar, self.nextCharPos = self.charIter.next()
		except StopIteration:
			self.nextChar = ''
	def peek(self):
		return self.nextChar
	def __iter__(self): 
		while True: 
			yield self.next()
	def next(self):
		c, p = self.nextChar, copy.copy(self.nextCharPos)
		if c == '': raise StopIteration
		self._advanceStoredNextChar()
		if c == '\r':	# returns \n not \r
			if self.nextChar == '\n': #only return on \n for dos endings
				self._advanceStoredNextChar()
			c = '\n'
		return c, p
import string
class NexusToken:
	tokenBreakers = r';\'()]{}/\\,:=*"`+<>-'
	singleTokenChars = r'(){}"]/\\,;:=*`+<>-'
	identityTrans = string.maketrans('', '')
	def __init__(self, nexusCharStream = None):
		self.comments = []
		if nexusCharStream:
			c = ''
			while c.strip() =='':
				c, sPos = nexusCharStream.next()
				if c == '[':
					self.comments.append(self._skipComment(nexusCharStream))
					c = ''
			self.startPos = sPos
			if c == "'": 
				self._readRestOfSingleQuoted(nexusCharStream)
			elif NexusToken.singleTokenChars.find(c) == -1:
				self._readRestOfToken(c, nexusCharStream)
			else:
				self.chars = c
				self.endPos = sPos
	def _readRestOfSingleQuoted(self, nexusCharStream):
		tok = cStringIO.StringIO()
		try:
			while True:
				c, self.endPos = nexusCharStream.next()
				if c == "'":
					if nexusCharStream.peek() == "'":
						c, self.endPos = nexusCharStream.next()
					else:
						break
				tok.write(c)
		except StopIteration:
			raise NexusOpenQuoteError, NexusOpenQuoteError(self.startPos)
		self.chars = tok.getvalue()
	def _skipComment(self, nexusCharStream):
		cmt = cStringIO.StringIO()
		try:
			while True:
				c, self.endPos = nexusCharStream.next()
				if c == ']':
					break
				else:
					if c == '[':
						cmt.write('[%s]' % self._skipComment(nexusCharStream))
					else:
						cmt.write(c)
		except StopIteration:
			raise NexusOpenCommentError, NexusOpenCommentError(self.startPos)
		return cmt.getvalue()
	def _readRestOfToken(self, c, nexusCharStream):
		sPos = self.startPos
		tok = cStringIO.StringIO()
		try:
			while True:
				if c == '[':
					self.comments.append(self._skipComment(nexusCharStream))
			 	else:
			 		tc = c == '_' and ' ' or c
					tok.write(tc)
					self.endPos = sPos
				n = nexusCharStream.peek()
				if n.strip() =='' or NexusToken.tokenBreakers.find(n) != -1:
					break
				c, sPos = nexusCharStream.next()
		except StopIteration: pass
		self.chars = tok.getvalue()
	def __str__(self):
		return self.chars
	def debugStr(self):
		return '->%s<- (from %s to %s)' % (self.chars, self.startPos, str(self.endPos))
	def __eq__(self, other):
		if other.__class__ == NexusToken:
			return (self.chars, self.startPos, str(self.endPos)) == (other.chars, other.startPos, str(other.endPos))
		return str(self).upper() == str(other).upper()
	def __ne__(self, other):return not self == other
	def escapeString(s):
		withoutSpecial = s.translate(NexusToken.identityTrans, NexusToken.tokenBreakers + '_[')
		if withoutSpecial == s:
			wsSplit = s.split()
			if len(wsSplit) == 1:
				return s
			singleSpaceSplit = s.split(' ')
			for i in singleSpaceSplit:
				if len(i.split()) > 1:
					return "'%s'" % s
			return '_'.join(singleSpaceSplit)
		return "'%s'" % s
	escapeString = staticmethod(escapeString) 
	def getSpecialComments(self):
		return [i for i in self.comments if (len(i)> 0 and (i[0] == '&' or i[0] == '!'))]
		
class NexusCommand:
	def __init__(self, nameToken, tokenStream):
		self.name = nameToken
		self.startPos = nameToken.startPos
		self.optionList = []
		for t in tokenStream:
			if t == ';':
				self.endPos = t.endPos
				return
			self.optionList.append(t)
		raise NexusOpenCommandError, NexusOpenCommandError(str(nameToken), nameToken.startPos)
	def __str__(self):
		s = ' '.join([str(o) for o in self.optionList])
		return '%s\n\t%s;' % (self.name, s)

class NexusCommandStream:
	def __init__(self, fileObj):
		self.nexusTokStream = NexusTokenStream(fileObj)
		self.firstCmd = True
		self.nextCommandName = None
	def next(self):
		cmdName = self.getNextCommandName()
		self.nextCommandName = None
		return NexusCommand(cmdName, self.nexusTokStream)
	def __iter__(self): 
		while True: 
			yield self.next()
	def getNextCommandName(self):
		if self.nextCommandName == None:
			self.nextCommandName = self.nexusTokStream.next()
			if self.firstCmd == True and self.nextCommandName == '#NEXUS':
				self.nextCommandName = self.nexusTokStream.next()
			self.firstCmd == False
			while self.nextCommandName == ';':
				self.nextCommandName = self.nexusTokStream.next()
		return self.nextCommandName
	def demandCommandEnd(self):
		self.nexusTokStream.demandToken(';')
		self.nextCommandName = None
	def getTokenStream(self):
		'''Asking for the token stream invalidates the nextCommandName field.'''
		self.nextCommandName = None
		return self.nexusTokStream
	def skipCommand(self):
		cmdName = self.getNextCommandName()
		self.nextCommandName = None
		for t in self.nexusTokStream:
			if t == ';':
				return

class NexusBlock:
	cmdHandlers = []
	def __init__(self, beginCmd = None, commandStream = None, previousBlocks = None):
		if beginCmd != None:
			self.name = str(beginCmd.optionList[0])
			self.startPos = beginCmd.name.startPos
		if commandStream != None:
			self.prepareToRead(previousBlocks or [])
			self.parseBlock(commandStream)
	def prepareToRead(self, previousBlocks):
		'''Called after the class of the block is determined.  Hook for initialization of derived classes).'''
		pass
	def endBlockEncountered(self): pass
	def skipBlock(self, commandStream):
		self.commands = []
		try:
			while True:
				nextName = commandStream.getNextCommandName()
				if str(nextName).upper() == 'END':
					c = commandStream.next()
					self.endPos = c.endPos
					break
				else:
					commandStream.skipCommand()
		except StopIteration:
			raise NexusOpenBlockError, NexusOpenBlockError(self.name, self.startPos)
	def parseBlock(self, commandStream):
		self.commands = []
		try:
			while True:
				nextName = commandStream.getNextCommandName()
				if str(nextName).upper() == 'END':
					c = commandStream.next()
					self.endPos = c.endPos
					self.endBlockEncountered()
					break
				self.interpretCommand(nextName, commandStream)
		except StopIteration:
			raise NexusOpenBlockError, NexusOpenBlockError(self.name, self.startPos)
	def interpretCommand(self, nextName, commandStream):
		for ch in self.__class__.cmdHandlers:
			if ch.attemptRead(nextName, commandStream, self, self):
				return
		self.commands.append(commandStream.next())

	def getCommandStream(self): return self.commandStream
	def getTokenStream(self): return self.commandStream.getTokenStream()			
	def __str__(self):
		s =  self.__dict__.has_key('commands') and '\n  '.join([str(c) for c in self.commands]) or ''
		return 'BEGIN %s;\n  %s\nEND;' % (self.name, s)

class NexusTokenStream:
	import re
	_mantissaPattern = re.compile(r'\A\d+\.?\d*[eE]?\d*\Z|\A\.\d+[eE]?\d*\Z')
	_fullFloatPattern = re.compile(r'\A-(\d+\.?\d*|\.\d+)[eE]-?\d+\Z|\A(\d+\.?\d*|\.\d+)[eE]-\d+\Z')
	def __init__(self, fileObj):
		self.charByCharToken = None
		self._nexusCharStream = NexusCharStream(fileObj)
	def next(self):
		return NexusToken(self._nexusCharStream)
	def nextChar(self):
		while self.charByCharToken == None or self.cbcIndex >= len(str(self.charByCharToken)):
			self.charByCharToken = self.next()
			self.cbcIndex = 0
		c = str(self.charByCharToken)[self.cbcIndex]
		self.cbcIndex += 1
		return c
	def cbcTok(self):
		return self.charByCharToken
	def cbcMoreCharacters(self):
		return self.charByCharToken != None and self.cbcIndex < len(str(self.charByCharToken))
	def demandToken(self, s):
		t = self.next()
		if str(t).upper() != s:
			raise NexusMissingTokenError(s, t)
	def readUntil(self, s):
		tList = []
		t = self.next()
		while str(t).upper() != s:
			tList.append(t)
			t = self.next()
		return tList
	def __iter__(self): 
		while True: 
			yield self.next()
	def nextFloat(self):
		''' '''
		t = self.next()
		sForm = str(t)
		makeNegative = sForm == '-'
		if makeNegative:
			breaker = self._nexusCharStream.peek()
			if not (breaker.isdigit() or breaker == '.'):  raise NexusMissingTokenError('a number', t)
			t = self.next()
			sForm = str(t)
		if not NexusTokenStream._mantissaPattern.match(sForm) and (makeNegative or sForm[0] != '-' or not NexusTokenStream._mantissaPattern.match(sForm[1:])):
			if NexusTokenStream._fullFloatPattern.match(sForm):
				return t
			raise NexusMissingTokenError('a number', t)
		if t.chars[-1].upper() == 'E':
			breaker = self._nexusCharStream.peek()
			if not (breaker.isdigit() or breaker == '-'):  raise NexusMissingTokenError('a number', t)
			ex = self.nextInteger()
			t.chars = t.chars + ex.chars
		if makeNegative:
			t.chars = '-' + t.chars	
		return t
	def nextInteger(self):
		''' returns a NexusToken with chars holding an integer from the next token in the stream (or two tokens if the number is negative)'''
		t = self.next()
		sForm = str(t)
		makeNegative = sForm == '-'
		if makeNegative:
			if not self._nexusCharStream.peek().isdigit():  raise NexusMissingTokenError('an integer', t)
			t = self.next()
			sForm = str(t)
		if not sForm.isdigit() and (makeNegative or sForm[0] != '-' or not sForm[1:].isdigit()):
			raise NexusMissingTokenError('an integer', t)
		if makeNegative:
			t.chars = '-' + t.chars
		return t

class NexusBlockStream:
	def __init__(self, fileObj, blockHandlers = {}, skipUnknownBlocks = True, previousBlocks = None):
		self.blockHandlers = blockHandlers
		self.commandStream = NexusCommandStream(fileObj)
		self.skipUnknownBlocks = skipUnknownBlocks
		self.previousBlocks = previousBlocks or []
	def next(self): 
		handlerClass = None
		while handlerClass == None:
			beginCmd = self.commandStream.next()
			if beginCmd.name != 'BEGIN':
				raise NexusBareCommandError, NexusBareCommandError(beginCmd)
			if len(beginCmd.optionList) != 1:
				raise NexusError, NexusError(beginCmd.name.endPos, None, 'Expecting block name and then a semi colon after the BEGIN command')
			blockName = str(beginCmd.optionList[0]).upper()
			handlerClass = self.blockHandlers.get(blockName)
			if handlerClass == None:
				if self.skipUnknownBlocks:
					NexusParsing.statusMessage('%s Block is not currently supported.  Skipping...' % blockName)
					self.skipBlock(self.commandStream)
				else:
					handlerClass = NexusBlock
		b = handlerClass(beginCmd, self.commandStream, copy.copy(self.previousBlocks))
		self.previousBlocks.append(b)
		return b
	def skipBlock(self, commandStream):
		while self.commandStream.next().name != 'END':
			pass
	def __iter__(self): 
		while True: 
			yield self.next()
	def getCommandStream(self):
		return self.commandStream
	def getTokenStream(self): 
		return self.commandStream.getTokenStream()
		
def iterNexusTokenizeStr(s):
	sStream = cStringIO.StringIO(str(s))
	return NexusTokenStream(sStream)
def getNexusTokenObjects(s):
	return [i for i in iterNexusTokenizeStr(s)]
def getNexusTokens(s):
	return [str(i) for i in iterNexusTokenizeStr(s)]

if __name__ == '__main__':
	import doctest, sys
	doctest.testmod(sys.modules[__name__])

if 0:
		#	s = '''#NEXUS''t[hi]e_s't[hiagain]'' 's\r\nf
		# a;b\d(e)f]g{h}i/k,l:m=n*o"p`q+r<s>t-ua a(a a)a a{a a}a a"a a]a a/a a\\a a,a a;a a:a a=a a*a a`a a+a a<a a>a a-a;'''
		s = '''begin a;
		b dga dga ad a;
		d;end;'''
		print s
		if 0:
			readString = cStringIO.StringIO(s)
			for i, p in CharPosStream(readString):
				print i, p
			readString = cStringIO.StringIO(s)
			for i, p in NexusCharStream(readString):
				print i, p
		readString = cStringIO.StringIO(s)
