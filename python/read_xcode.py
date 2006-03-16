#!/usr/bin/python2.3
from __future__ import generators
import sys, re, copy, os
from sets import Set

versionString = '0.01'

autoCompilerFlags = [] # filled based on platform (with no connection to what is in the xcode project)
compilerTranslations = {} # dictionary used to convert options recognized by mac's g++ (used by xcode) to the appropriate flag for the targetted compiler/platform
headersMatch = re.compile(r'.+\.h(pp)?$')
sourceMatch = re.compile(r'.+\.c(pp)?$')
			
def xmlWrap(elName, s):
	'''Wraps s in elName tags (does NOT escape s for xml)'''
	return '<%s>%s</%s>' % (elName, s, elName)
def tokenizeFile(xcodeFileName):
	inF = open(xcodeFileName, 'r')
	for line in inF:
		words = line.split()
		if len(words) > 0 and not words[0].startswith('//'): #lines beginning with // appear to be comments
			prev = words[0] 
			inQuote = prev.startswith('"')
			for w in words[1:]:
				if not inQuote:
					yield prev
					inQuote = w.startswith('"') and not (w.endswith('"') or w.endswith('";'))
					prev = w
				else:
					prev = prev + ' ' + w
					inQuote = not (w.endswith('"') or w.endswith('";'))
			if prev.endswith(';'): 	#last word in every non opening brace line seems to be a ; We'll strip it off
				prev = prev[:-1]
			yield prev

idStringRE = re.compile(r'[0-9A-F]{24,24}')
def isIDString(s):
	return idStringRE.match(s)

idLookup = {}
knownBuildStyles = []
def changeXCodeObjType(obj, tttMapping):
	if tttMapping.get(obj.__class__) != None:
		obj.__class__ = tttMapping[obj.__class__]
def translateFromXCodeTypeIterable(iterable, tttMapping):
	for o in iterable:
		if isinstance(o, XCodeDOMObject) or isinstance(o, SourceFileTree):
			translateFromXCodeType(o, tttMapping)
		elif isinstance(o, list):
			translateFromXCodeTypeIterable(o, tttMapping)
		elif isinstance(o, dict):
			translateFromXCodeTypeIterable(o.iteritems(), tttMapping)
		else:
			changeXCodeObjType(o, tttMapping)
def translateFromXCodeType(o, tttMapping):
	if isinstance(o, SourceFileTree):
		o.translateFileTreeFromXCodeType(tttMapping)
	else:
		translateFromXCodeTypeIterable([c for c in o.__dict__.itervalues()], tttMapping)
		changeXCodeObjType(o, tttMapping)

class XCodeDOMObject:
	indentationLevel = 0
	def __init__(self, tokenStream):
		while(True):
			label = tokenStream.next()
			if label == '}': return
			e = tokenStream.next()
			if e != '=': 
				print e, ' !!!! '
				raise ValueError, 'Expecting ='
			if label == 'isa':
				isa = tokenStream.next()
				self.__class__ = stringToClassDict[isa]
				if isa == 'PBXBuildStyle':
					knownBuildStyles.append(self)
			else:
				o = readDOMObject(tokenStream)
				if isIDString(label):
					idLookup[label] = o
				self.__dict__[label] = o
	def resolveReferences(self):
		for label, val in self.__dict__.iteritems():
			if isinstance(val, str):
				if isIDString(val):
					self.__dict__[label] = idLookup[val]
			elif isinstance(val, list):
				for i in range(0,len(val)):
					o = val[i]
					if isinstance(o, str) and isIDString(o):
						val[i] = idLookup[o]
			else:
				val.resolveReferences()
	def resolvePaths(self, parPath):
		for label, val in self.__dict__.iteritems():
			if not isinstance(val, str):
				if isinstance(val, list):
					for v in val:
						if not isinstance(v, str):
							v.resolvePaths(parPath)
				else:
					val.resolvePaths(parPath)
	def __str__(self):
		tabStr = ''.join(['\t' for i in range(0, XCodeDOMObject.indentationLevel)])
		XCodeDOMObject.indentationLevel += 1
		objList = [''.join([' (__class__ = ', self.__class__.__name__, ')\n'])]
		for l, o in self.__dict__.iteritems():
			if isinstance(o, list):
				objList.append(''.join([tabStr, l, ' LIST\n']))
				XCodeDOMObject.indentationLevel += 1
				for subObj in o:
					objList.append(''.join([tabStr, '\t', str(subObj), '\n']))
				XCodeDOMObject.indentationLevel -= 1
			else:
				objList.append(''.join([tabStr, l, ' ', str(o), '\n']))
		XCodeDOMObject.indentationLevel -= 1
		return ''.join(objList)

class PBXBuildFile(XCodeDOMObject): pass
class PBXBuildStyle(XCodeDOMObject):
	def standardizeNames(self):
		if self.name == 'Development': self.name = 'Debug'
		elif self.name == 'Deployment': self.name = 'Release'
class PBXProject(XCodeDOMObject):
	'''the root object of an xcode project file.'''
	def expandTargets(self, buildStyleList):
		et = []
		for t in self.targets:
			et.extend(t.createResolvedCopies(buildStyleList))
		self.targets = et
class PBXSourcesBuildPhase(XCodeDOMObject): pass
class PBXHeadersBuildPhase(XCodeDOMObject): pass
class PBXFileReference(XCodeDOMObject):
	def resolvePaths(self, parPath):
		if parPath != '' and self.sourceTree == '"<group>"':
			self.path = parPath + self.path
			self.sourceTree = 'python_set' + self.sourceTree
class PBXFrameworksBuildPhase(XCodeDOMObject): pass
class PBXRezBuildPhase(XCodeDOMObject): pass
class PBXGroup(XCodeDOMObject):
	def resolvePaths(self, parPath):
		if self.sourceTree == '"<group>"':	
			if parPath == '': return	# we might recurse into a relative to parent group without the proper parent info.  If so, we bail out
		elif parPath.startswith('_python_set'): return  # we've already hit this node in the dom tree bail out
		else: parPath = self.sourceTree + '/' 	# we aren't relative to the group, so we'll use whatever is in the sourceTree field 
		pathAdd = self.__dict__.get('path', '') # path = "ame of this node" this info might not exist, or may be ""
		if pathAdd == '""': 
			pathAdd = ''
		  # append to the parent path and propogate this full path to children
		if len(pathAdd) > 0:
			newpath = parPath + pathAdd + '/'
		else:
			newpath = parPath
		if self.sourceTree == '"<group>"':
			self.path = newpath
			self.sourceTree = '_python_set' + self.sourceTree
		children = self.__dict__.get('children', [])
		for c in children:
			c.resolvePaths(copy.copy(newpath))
class PBXCopyFilesBuildPhase(XCodeDOMObject): pass
sourceTreeAlias = {}
def filterForTargetCompiler(s):
	global compilerTranslations
	inFlags = s.split()
	o = []
	for i in inFlags:
		if compilerTranslations.get(i) != None:
			o.append(compilerTranslations[i])
		else:
			o.append(i)
	return  ' '.join(o)

class BuildSettings: 
	'''xcode has two levels of build settings "build styles" and target-specific.  Each target is based on a build-style and tweaks those settings
		This class merges the two levels to make a target-specific set of flags
		buildStyleName = 'Debug|Release'
		PRODUCT_NAME = 'console_phycas|phycas|mth_dev_phycas|pol_dev_phycas'
		isDebugTarget = True|False
		headerSearchPaths = [header]
		optimizationLevel = integer
		'''
	def __init__(self, primaryBuildStyle, secondaryBS):
		global autoCompilerFlags
		self.__dict__ = copy.deepcopy(primaryBuildStyle.buildSettings.__dict__)
		for k,v in secondaryBS.__dict__.iteritems():
			if self.__dict__.has_key(k):
				if len(v) > 1 and v[0] == '"' and v[1] != '"':
					s = self.__dict__[k]
					self.__dict__[k] = s[0:-1] + ' ' + v[1:]
			else:
				self.__dict__[k] = copy.copy(v)
		self.buildStyleName = primaryBuildStyle.name
		self.isDebugTarget = self.buildStyleName == 'Debug'
		self.removeQuotes('GCC_PREFIX_HEADER')
		gphl = self.GCC_PREFIX_HEADER.split('/')
		if gphl[0] == 'PHYCAS_ROOT' or gphl[0] == '$PHYCAS_ROOT': 
			self.GCC_PREFIX_HEADER = '/'.join(gphl[1:])
		self.removeQuotes('HEADER_SEARCH_PATHS')
		self.headerSearchPaths = self.HEADER_SEARCH_PATHS.split()
		del self.__dict__['HEADER_SEARCH_PATHS']
		self.removeQuotes('OTHER_CPLUSPLUSFLAGS')
		self.adjustForCompiler('OTHER_CPLUSPLUSFLAGS')
		self.removeQuotes('WARNING_CFLAGS')
		self.adjustForCompiler('WARNING_CFLAGS')
		self.removeQuotes('OPTIMIZATION_CFLAGS')
		self.adjustForCompiler('OPTIMIZATION_CFLAGS')
		if len(self.OPTIMIZATION_CFLAGS) == 3 and self.OPTIMIZATION_CFLAGS.startswith('-O'):	#convert string of optimization level to number
			self.optimizationLevel = int(self.OPTIMIZATION_CFLAGS[2])
		else:
			if self.isDebugTarget: self.optimizationLevel = 0
			else: self.optimizationLevel = 3
		del self.__dict__['OPTIMIZATION_CFLAGS']
		self.OTHER_CPLUSPLUSFLAGS = self.OTHER_CPLUSPLUSFLAGS + ' ' + ' '.join(autoCompilerFlags)
		
	def removeQuotes(self, s):
		strForm = self.__dict__.get(s, '""')
		if strForm[0] == '"':
			if strForm[-1] == '"':
				self.__dict__[s] = strForm[1:-1]
			else:
				self.__dict__[s] = strForm[1:]
		else:
			if strForm[-1] == '"':
				self.__dict__[s] = strForm[:-1]
			else:
				self.__dict__[s] = strForm

	def adjustForCompiler(self, s):
		if self.__dict__.get(s) != None:
			self.__dict__[s] = filterForTargetCompiler(self.__dict__[s])


def verifyDirs(p):
	if not os.path.exists(p):
		os.makedirs(p)
def findCommonParent(l):
	cp = os.path.commonprefix(l)
	lastInd = cp.rfind('/')
	if lastInd < 0: return ''
	return cp[:lastInd+1]
	
class SourceFileTree:
	''' corresponds to a directory in the source file-directory tree.
	name holds the directory name.
	files is a set of files in this directory.
	subdirs maps a subdir name to the SourceFileTree object that describes that sub directory.
	sourceTreeRootName is the name of parent of this source tree.'''
	knownSourceTrees = {}
	def __init__(self, n = '', pathFromRoot = ''):
		self.files = Set()
		self.subdirs = {}
		self.name = n
		self.sourceTreeRootName = ''
		if pathFromRoot == '': self.pathFromSourceTreeRoot = n
		else: self.pathFromSourceTreeRoot = pathFromRoot + '/' + n
	def add(self, path):
		splitpath = path.split('/')
		nextName = splitpath[0]
		if len(splitpath) == 1:
			self.files.add(nextName)
		else:
			toAdd = '/'.join(splitpath[1:])
			if self.subdirs.get(nextName) == None:
				self.subdirs[nextName] = SourceFileTree(nextName, self.pathFromSourceTreeRoot)
			self.subdirs[nextName].add(toAdd)
	def getFilePaths(self, fileMatch = None, dirSep = '/', prefix = ''):
		''' returns full path to all files '''
		ret = []
		if len(prefix) > 0: prefix = prefix + dirSep
		for f in self.files:
			if fileMatch == None or fileMatch.match(f): 
				ret.append(prefix + f)
		for k, sd in self.subdirs.iteritems():
			ret.extend(sd.getFilePaths(fileMatch, dirSep, prefix + k))
		return ret
	def getSubDirs(self): return [d for d in self.subdirs.itervalues()]
	def getFileList(self): return [f for f in self.files]
	def getFileAndSubDirDict(self):
		'''returns (list, dict) holding file-list and dict of sub-directory to SourceFileTree.'''
		return self.getFileList(), self.subdirs
	def setSourceTreeRootName(self, sourceTreeRootName):
		for sDir in self.subdirs.itervalues():
			sDir.setSourceTreeRootName(sourceTreeRootName)
		if self.pathFromSourceTreeRoot.startswith(sourceTreeRootName + '/'):
			self.pathFromSourceTreeRoot = self.pathFromSourceTreeRoot[len(sourceTreeRootName) + 1:]
		self.sourceTreeRootName = sourceTreeRootName
	def translateFileTreeFromXCodeType(self, tttMapping):
		if tttMapping.get(self.__class__) != None:
			self.__class__ = tttMapping[self.__class__]
		f, sd = self.getFileAndSubDirDict()
		for v in sd.itervalues():
			v.translateFileTreeFromXCodeType(tttMapping)
	def catalogSourceTrees(self):
		SourceFileTree.knownSourceTrees = self.subdirs
		allSourceTreeSubdirs = []
		for sourceTreeRootName, sDir in self.subdirs.iteritems():
			sDir.setSourceTreeRootName(sourceTreeRootName)
			allSourceTreeSubdirs.extend(sDir.getSubDirs())
		return allSourceTreeSubdirs
	def getMWerksGroupedFileRef(self, f, aTargetName, indentLevel):
		nextIndent = indentLevel + '\t'
		return '%s<FILEREF>\n%s%s\n%s<PATHTYPE>PathRelative</PATHTYPE>\n%s%s\n%s<ACCESSPATH></ACCESSPATH>\n%s%s\n%s<PATHFORMAT>Unix</PATHFORMAT>\n%s</FILEREF>\n' % (
			indentLevel, 
			nextIndent, xmlWrap('TARGETNAME', aTargetName),
			nextIndent,
			nextIndent, xmlWrap('PATHROOT', self.sourceTreeRootName),
			nextIndent,
			nextIndent, xmlWrap('PATH', self.pathFromSourceTreeRoot + '/' + f),
			nextIndent,	indentLevel)
	def getPathWithSourceTreeRoot(self):
		return self.sourceTreeRootName + '/' + self.pathFromSourceTreeRoot
	def writeMWerksFileGroups(self, out, targets, indentLevel = '\t\t\t'):
		fileList, subDirDict = self.getFileAndSubDirDict()
		out.write('%s<GROUP>%s\n' % (indentLevel, xmlWrap('NAME', self.name)))
		fullPath = self.getPathWithSourceTreeRoot()
		for sd in subDirDict.itervalues():
			sd.writeMWerksFileGroups(out, targets, indentLevel+'\t')
		for f in fileList:
			for t in targets:
				if t.usesFile(fullPath + '/' + f):
					aTargetName = t.getMWerksTargetName(False)
					out.write(self.getMWerksGroupedFileRef(f, aTargetName, indentLevel + '\t'))
					break
		out.write('%s</GROUP>\n' % indentLevel)

class PBXNativeTarget(XCodeDOMObject):
	def createResolvedCopies(self, listOfBuildStyles):
		resCopies = []
		for bs in listOfBuildStyles:
			c = copy.copy(self)
			c.buildSettings = BuildSettings(bs, c.buildSettings) # note the higher level (non-target specific) build settings have precedence
			if c.buildSettings.isDebugTarget or (not c.buildSettings.PRODUCT_NAME.startswith('pol_') and not c.buildSettings.PRODUCT_NAME.startswith('mth_')):
				resCopies.append(c)
		return resCopies
	def getFilePaths(self): return self.headers + self.sources 	# consolidateFiles must be called first
	def getSources(self): 	return self.sources  				# consolidateFiles must be called first
	def getHeaders(self): 	return self.headers					# consolidateFiles must be called first
	def consolidateFiles(self, fileClassifier):
		for bp in self.buildPhases:
			if isinstance(bp, PBXSourcesBuildPhase): self.sources = [f.fileRef.path for f in bp.files]
			elif isinstance(bp, PBXHeadersBuildPhase): self.headers = [f.fileRef.path for f in bp.files]
		usedFilePaths = self.getFilePaths()
		for f in usedFilePaths:
			fileClassifier.add(f)
	def usesFile(self, fp):
		filePaths = self.getFilePaths()
		return filePaths.count(fp) > 0
		
stringToClassDict = {
	'PBXBuildFile' : PBXBuildFile, 
	'PBXBuildStyle' : PBXBuildStyle, 
	'PBXProject' : PBXProject, 
	'PBXSourcesBuildPhase' : PBXSourcesBuildPhase, 
	'PBXHeadersBuildPhase' : PBXHeadersBuildPhase, 
	'PBXFileReference' : PBXFileReference, 
	'PBXFrameworksBuildPhase' : PBXFrameworksBuildPhase, 
	'PBXRezBuildPhase' : PBXRezBuildPhase, 
	'PBXGroup' : PBXGroup,
	'PBXCopyFilesBuildPhase' : PBXCopyFilesBuildPhase,
	'PBXNativeTarget' : PBXNativeTarget
	}
	
def readDOMObject(tokenStream):
	for token in tokenStream:
		if token == '{':
			return XCodeDOMObject(tokenStream)
		if token == '(':
			group = []
			for parensToken in tokenStream:
				if parensToken == ')': return group
				if not parensToken.endswith(','): 
					print parensToken, ' !!!!'
					raise ValueError, 'Expecting ) or token ending with a ,'
				group.append(parensToken[:-1])
		return token
	return None
	
def readXCodeProject(xcodeFileName, additionalCompilerFlags = None, compilerTranslateDict = None, typeToTypeMap = None):
	global knownBuildStyles, autoCompilerFlags, compilerTranslations
	knownBuildStyles = []
	autoCompilerFlags = additionalCompilerFlags or []
	compilerTranslations = compilerTranslateDict or {}
	nestingLevel = -1 #file starts with { which will take us to nestingLevel 0
	tokens = [i for i in tokenizeFile(xcodeFileName)]
	tokenIter = iter(tokens)
	domTree = readDOMObject(tokenIter)
	domTree.resolveReferences()
	domTree.resolvePaths('')
	for bs in knownBuildStyles: bs.standardizeNames()
	domTree.rootObject.expandTargets(knownBuildStyles)
	allFiles = SourceFileTree()
	for t in domTree.rootObject.targets:
		t.consolidateFiles(allFiles)
	sources = allFiles.catalogSourceTrees()
	tttMap = typeToTypeMap or []
	if len(tttMap) > 0:
		translateFromXCodeType(allFiles, tttMap)
		translateFromXCodeType(domTree.rootObject, tttMap)
	return (domTree.rootObject, sources)
			
if __name__ == '__main__':
	if len(sys.argv) < 2:	print 'Usage: %s xcode_project' % sys.argv[0]; sys.exit(1)
	xcodeFileName = sys.argv[1]
	project, files = readXCodeProject(xcodeFileName)
	print project
