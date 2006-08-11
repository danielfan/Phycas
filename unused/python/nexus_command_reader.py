#!/usr/bin/python
import copy, sys
from nexus_parser import *
from nexus_primitives import *
class NexusSubCommandReader(object):
	def isPositiveNumber(nscr, n, tok, obj):
		if n > 0:
			return True
		raise NexusUnexpectedTokenError('%s (a positive number was expected)' % nscr.name, tok)
	isPositiveNumber = staticmethod(isPositiveNumber)
	def __init__(self, n, default = None, validation = None):
		self.name = n
		self.default = default
		self.validation = validation
		self.castFunc = str
	def canRead(self, capStr):
		return self.name.upper() == capStr
	def setToDefault(self, d):
		d[self.name] = self.default
		
	def readSubCommand(self, cList, i, dictToAddTo):
		if i+1 >= len(cList) or cList[i+1] != '=':
			raise NexusMissingTokenError('=', cList[i])
		i += 2
		if i >= len(cList):
			raise NexusMissingTokenError(self.getTypeDescrip(), cList[i-1])
		return self.readValue(cList, i, dictToAddTo)

	def readValue(self, cList, i, dictToAddTo):
		try:
			n = self.castFunc(str(cList[i]))
		except ValueError:
			raise NexusUnexpectedTokenError('%s (%s was expected)' %(str(cList[i]), self.getTypeDescrip()), cList[i])
		self.setAndValidate(n, cList[i], dictToAddTo)
		return i+1
	def setAndValidate(self, n, tok, dictToAddTo):
		backup = None
		hadKey = dictToAddTo.has_key(self.name)
		if hadKey: backup = n
		dictToAddTo[self.name] = n
		try:
			if self.validation != None:
				self.validation(self, n, tok, dictToAddTo)
		except: # catch all exceptions
			typ, val, trac = sys.exc_info()
			if hadKey: 
				dictToAddTo[self.name] = backup
			else:
				del dictToAddTo[self.name]
			raise typ, val, trac # reraise

class NexusBoolSubCommandReader(NexusSubCommandReader):
	def __init__(self, n, default, validation = None):
		super(NexusBoolSubCommandReader, self).__init__(n, default, validation)
		self.castFunc = bool
	def getTypeDescrip(self): return 'True/False'
	def canRead(self, capStr):
		nup = self.name.upper()
		if nup == capStr:
			self.converse = False
			return True
		if 'NO' + nup == capStr:
			self.converse = True
			return True
		return False
	def readSubCommand(self, cList, i, dictToAddTo):
		if self.converse:
			self.setAndValidate(False, cList[i], dictToAddTo)
			return i + 1
		if i+1 >= len(cList) or cList[i+1] != '=':
			self.setAndValidate(True, cList[i], dictToAddTo)
			return i + 1
		NexusSubCommandReader.readSubCommand(self, cList, i, dictToAddTo)
	def readValue(self, cList, i, dictToAddTo): 
		raise TypeError, 'NexusBoolSubCommandReader.readValue called'
		

class NexusNumericSubCommandReader(NexusSubCommandReader):
	def __init__(self, n, default, validation = None):
		super(NexusNumericSubCommandReader, self).__init__(n, default, validation)

class NexusIntSubCommandReader(NexusNumericSubCommandReader):
	def __init__(self, n, default, validation = None):
		super(NexusIntSubCommandReader, self).__init__(n, default, validation)
		self.castFunc = int
	def getTypeDescrip(self): return 'Integer'
	
class NexusFloatSubCommandReader(NexusNumericSubCommandReader):
	def __init__(self, n, default, validation = None):
		super(NexusFloatSubCommandReader, self).__init__(n, default, validation)
		self.castFunc = float
	def getTypeDescrip(self): return 'Real numbr'
	
class NexusStringSubCommandReader(NexusSubCommandReader):
	def __init__(self, n, default, validation = None):
		super(NexusStringSubCommandReader, self).__init__(n, default, validation)
	def getTypeDescrip(self): return 'String'
	
class NexusChoiceSubCommandReader(NexusSubCommandReader):
	def __init__(self, n, choices, default, validation = None):
		self.sec_validation = validation
		self.choices = choices
		super(NexusChoiceSubCommandReader, self).__init__(n, default, NexusChoiceSubCommandReader.choiceValidation)
	def choiceValidation(self, c, tok, dictToAddTo):
		capC = c.upper()
		for pc in self.choices:
			if capC == pc.upper():
				if self.sec_validation != None:
					self.sec_validation(self, c, tok, dictToAddTo)
				return True
		raise NexusMissingTokenError('|'.join(self.choices), tok)
	def getTypeDescrip(self): return '|'.join(self.choices)
	
class NexusUnsupportedSubCommand(NexusSubCommandReader):
	def __init__(self, n, default = None):
		super(NexusUnsupportedSubCommand, self).__init__(n, default, None)
	def readSubCommand(self, cList, i, dictToAddTo):
		raise NexusUnsupportedError('The ' + self.name + ' sub-command', cList[i])

class NexusCommandReader(object):
	def __init__(self, n, subcmds = None, postValidation = None, putVarsInCmdNamedDict = False, readerToCreate = None):
		self.name = n
		self.subCommands = subcmds or []
		self.putVarsInCmdNamedDict = putVarsInCmdNamedDict
		self.readerToCreate = readerToCreate
		if postValidation != None:
			self.postValidation = postValidation
	def attemptRead(self, cName, cStream, obj, blockObj):
		if str(cName).upper() != self.name.upper():
			return False
		if self.readerToCreate != None:
			readerObj = self.readerToCreate(blockObj)
			return readerObj.readCommand(cName, cStream, obj, blockObj)
		cList = cStream.next().optionList
		i = 0
		scList = copy.copy(self.subCommands)
		if self.putVarsInCmdNamedDict:
			obj.__dict__[self.name] = {}
			oToPass = obj.__dict__[self.name]
		else:
			oToPass = obj.__dict__
		for sc in self.subCommands:
			sc.setToDefault(oToPass)
		while len(cList) > i:
			i = self.readSubCommands(cList, i, scList, oToPass)
		if self.__dict__.has_key('postValidation'):
			self.postValidation(blockObj, self, oToPass, cName)
		return True
	def readSubCommands(self, cList, i, scList, dictToAddTo):
		capStr = str(cList[i]).upper()
		for scReader in scList:
			if scReader.canRead(capStr):
				scList.remove(scReader)
				return scReader.readSubCommand(cList, i, dictToAddTo)
		raise NexusUnexpectedTokenError(self.name, cList[i])
class NexusUnsupportedCommand(NexusCommandReader):
	def __init__(self, n):
		super(NexusUnsupportedCommand, self).__init__(n)
	def attemptRead(self, cName, cStream, obj, blockObj):
		if str(cName).upper() != self.name.upper():
			return False
		raise NexusUnsupportedError('The ' + self.name + ' command', cName)

