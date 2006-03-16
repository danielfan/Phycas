#!/usr/bin/python2.3
# author mth
# some xml reading code taken from:
#	pyxml/test/test_sax.py.
# 	The header from that file
# 	regression test for SAX 2.0
# 
# $Id: sax_constructible.py,v 1.4 2005/07/20 04:05:18 mholder Exp $
import os
from xml.sax import make_parser, ContentHandler, \
					SAXException, SAXReaderNotAvailable, SAXParseException
try:
	make_parser()
except SAXReaderNotAvailable:
	raise ImportError("no XML parsers available")
	
from xml.sax.saxutils import XMLGenerator, escape, quoteattr, XMLFilterBase 
from xml.sax.expatreader import create_parser

import copy

	
def adict(**kwds):
	return kwds

def writeSAXConstructibleToXML(outStr, elName, atts, text, subElements):
	outStr.write('<%s' % elName)
	for k, v in atts.iteritems():
		outStr.write(' %s="%s"' % (k, v))
	outStr.write('>%s\n' % text)
	for element in subElements:
		element.writeXML(outStr)
	outStr.write('</%s>\n' % elName)
	
class XMLEchoHandler(ContentHandler):
	'''Prints element name and characters with indentation indicating element depth'''
	_nTabs = 0
	def startElement(self, name, attrs):
		for i in xrange(self._nTabs):
			print '\t',
		print name.encode('ascii')
		self._nTabs += 1
	def endElement(self, name):
		self._nTabs -= 1
	def characters(self, ch):
		print ch

class SkipSubTreeHandler(ContentHandler):
	'''Class to absorb SAX events for a subtree.'''
	def __init__(self, parser, parent, depth = 1):
		'''Sets self to the parser's ContentHandler until the subtree is read, and then installs parent as the ContentHandler.
		
		Pass depth or 0 if the object will receive a startElement for the element to ignore, and 1 if the startElement has 
		already been read.'''
		self._parser = parser
		self._parent = parent
		self._depth = depth
		self._parser.setContentHandler(self)
	def startElement(self, name, attrs):
		self._depth += 1
	def endElement(self, name):
		self._depth -= 1
		if self._depth == 0:
			self._parser.setContentHandler(self._parent)
	def characters(self, ch): pass

class ChildrenToHandle(object):
	'''LIsts attributes to read and maps children elements to SAXConstructible class to handle the element'''
	def __init__(self, **children):
		'''Take dictionary mapping attr, singleEl, and multiEl to a list and 2 dictionaries.'''
		self.attributesToRead = children.get('attr', [])
		self.singleElementDict = children.get('singleEl', {})
		self.multiElementDict = children.get('multiEl', {})
	def getChildElementsCTor(self, name):
		'''Returns the a sub-class of SAXConstructible to handle an element 'name' or None.'''
		childTypePair = self.singleElementDict.get(name)
		if childTypePair == None:
			return self.multiElementDict.get(name)
		return childTypePair
	def isSingleElement(self, name):
		'''Returns True if the 'name' is the name of an element that is expected as a child once'''
		return self.singleElementDict.get(name) != None
	def getUnion(self, other):
		a = copy.copy(self.attributesToRead)
		a.extend(other.attributesToRead)
		s = copy.copy(self.singleElementDict)
		s.update(other.singleElementDict)
		m = copy.copy(self.multiElementDict)
		m.update(other.multiElementDict)
		return ChildrenToHandle(attr = a, singleEl = s, multiEl = m)
	def __str__(self):
		a =  'attr = ' + str(self.attributesToRead)
		s = ', '.join([i for i in self.singleElementDict.iterkeys()])
		m = ', '.join([i for i in self.multiElementDict.iterkeys()])
		return a + '\nsingle elements = ' + s + '\nmulti elements = ' + m
class SAXConstructible(ContentHandler, object):
	''' Abstract class for Sax ContententHandlers which delegate parsing of elements to contained SAXConstructible objects.
	
		Methods startSelfElement, childDone are pure virtual'''
	initializeFromClassStatic = False # set to True to use the derived class's _attr, _singleEl, and _multiEl to intialize the _childrenToHandle field
	def __init__(self, elName, parseContext, childrenToHandle = None):
		''' Sets _elName to elName and calls startSelfElement if attrs were sent in parseContext.'''
		self._elName = elName
		if SAXConstructible.initializeFromClassStatic:
			self._childrenToHandle = self.__class__._CTH
		else:
			if childrenToHandle == None:	raise TypeError, 'must use childrenToHandle arg in SAXConstructible.__init__() when SAXConstructible.initializeFromClassStatic is False'
			self._childrenToHandle = childrenToHandle
		self.clear()
		self._parent = parseContext.get('parent')
		self._parser = parseContext.get('parser')
		if self._parser != None: self._parser.setContentHandler(self)
		if parseContext.get('attrs') != None:	
			self.callStartSelfElement(self._elName, parseContext.get('attrs'))

	def version(): return '0.1'
	version = staticmethod(version)
	
	def getSubParser(self, inParser, name, inAttrs):
		''' Creates a new ContentHandler of the class specified in _childrenToHandle.
		
			If an element of type name found in _childrenToHandle dictionary
			will be used to return a new object of the correct type. If it is not
			found a SkipSubTreeHandler will be returned.'''
		cTorArg = self._childrenToHandle.getChildElementsCTor(name)
		if cTorArg == None:	return SkipSubTreeHandler(inParser, self)
		parseContext = adict(parser = inParser, parent = self, attrs = inAttrs)
		return cTorArg(name, parseContext)
	
	def getElementName(self):
		''' Returns the elName passed in as to __init__.'''
		return self._elName
	
	def getObject(self):
		''' Returns self or None (if the element has not been read.'''
		if self._saxRead or self.hasChars():
			return self
		return None
	
	def hasChars(self):
		'''Returns True if self._rawChars is not empty.'''
		return len(self._rawChars) > 0
		
	def clear(self):
		''' Resets _rawChars, all elements and attributes and flags.'''
		self._saxRead = False
		for attTag in self._childrenToHandle.attributesToRead:
			self.__dict__[attTag] = ''
		for elTag in self._childrenToHandle.singleElementDict.iterkeys():
			self.__dict__[elTag] = None
		for elTag in self._childrenToHandle.multiElementDict.iterkeys():
			self.__dict__[elTag] = []
		self._allElements = []
		self._rawChars = ''
		
	def parseFile(self, xmlFile):
		''' Creates a new parser for fileName, using self as the initial ContentHandler.'''
		self._parent = None
		self._parser = create_parser()
		self._parser.setContentHandler(self)
		if isinstance(xmlFile, str):
			fileName = os.path.expanduser(xmlFile)
			xmlFile = open(fileName, 'r')
		self._parser.parse(xmlFile)
		return True

	
	def callStartSelfElement(self, name, attrs):
		''' Calls startSelfElement(), sets self._saxRead to True, and calls postStartSelfElement()'''
		self.startSelfElement(name, attrs)
		self._saxRead = True
		self.postStartSelfElement()
	
	def startSelfElement(self, name, attrs):
		''' Adds all attributes in _childrenToHandle.attributesToRead to self's dictionary '''
		for attTag in self._childrenToHandle.attributesToRead:
			self.__dict__[attTag] = str(attrs.get(attTag, ""))
	
	def postStartSelfElement(self): 
		''' Hook for post-startSelfElement() behavior'''
		pass
		
	def startElement(self, name, attrs):
		''' Spawns a subParser or calls startSelfElement if self has not parsed and name matches getElementName().'''
		if name == self._elName and self.getObject() == None:
			self.callStartSelfElement(name, attrs)
		else:
			self.subParser = self.getSubParser(self._parser, name, attrs)
			
	def characters(self, ch): 
		''' Strips and then append ch.'''
		stripped = ch.strip()
		if (len(stripped) > 0):
			self._rawChars += stripped.encode('ascii') 
			
	def childDone(self, obj, name): 
		''' Appends the obj to the list of _allElements and to the element data member 'name'.'''
		self._allElements.append(obj)
		if self._childrenToHandle.isSingleElement(name):
			self.__dict__[name] = obj
		elif self.__dict__.get(name) == None: 
			return
		else:
			self.__dict__[name].append(obj)
		
	def endSelfElement(self, name):
		''' Hook for post-endElement() behavior.'''
		pass
	
	def endElement(self, name): 
		''' Calls endSelfElement() then parent's childDone() method with getObject() as an argument, and sets the parent as the parser's ContentHandler '''
		if name != self._elName:
			raise ValueError, 'Unexpected name in endElement'
		self.endSelfElement(name)		
		if self._parent != None:
			self._parent.childDone(self.getObject(), self._elName)
			if self._parser != None:
				self._parser.setContentHandler(self._parent)
			self._parent = None
		elif self._parser != None:
			self._parser.setContentHandler(XMLEchoHandler())
			self._parser = None

	def getRawChars(self):
		''' Returns stripped concatenation of all characters.'''
		return self._rawChars
	
	def getEscapedChars(self):
		''' Returns stripped concatenation of all characters, with &, < and > escaped.'''
		return escape(self._rawChars)

	def __str__(self):
		''' Concatenates self._rawChars and string for of every member self._allElement.'''
		return self._rawChars + ''.join([str(x) for x in self._allElements])
	
	def writeXML(self, outStr):
		'''Writes the content of the object to XML preserving element order, but not order of characters relative to elements'''
		filledAttDict = {}
		for attTag in self._childrenToHandle.attributesToRead:
			if self.__dict__[attTag] != '':
				filledAttDict[attTag] = self.__dict__[attTag]
		writeSAXConstructibleToXML(outStr, self._elName, filledAttDict, self.getEscapedChars(), self._allElements)

	def ignorableWhitespace(self, content):
		'''Ignores content.'''
		pass

	def processingInstruction(self, target, data):
		'''Ignores the processing content.'''
		pass

class TextOnlyElement(SAXConstructible):
	'''Only handles character content.'''
	automaticallyConvertToString = False
	_CTH = ChildrenToHandle()
	def __init__(self, elName, parseContext, translateFunc = None):
		super(TextOnlyElement, self).__init__(elName, parseContext, TextOnlyElement._CTH)
		self._translateFunc = translateFunc
		
	def childDone(self, obj, n): 
		raise ValueError, 'Unexpected Child in TextOnlyElement'
	def endSelfElement(self, name):
		if self._translateFunc != None:
			self._rawChars = self._translateFunc(self._rawChars)
	def getObject(self):
		if TextOnlyElement.automaticallyConvertToString:
			return str(self)
		return self
		
class TranslateTextOnlyElement(SAXConstructible):
	'''Callable class that creates instances of TextOnlyElement with a stored translation function.'''
	def __init__(self, translateFunc):
		self._translateFunc = translateFunc
	def __call__(self, elName, parseContext): 
		return TextOnlyElement(elName, parseContext, self._translateFunc)
	
	
class IgnoreContentElement(SAXConstructible):
	'''Only handles character content.'''
	_CTH = ChildrenToHandle()
	def __init__(self, elName, parseContext):
		super(IgnoreContentElement, self).__init__(elName, parseContext, IgnoreContentElement._CTH)
	
