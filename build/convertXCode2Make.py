#!/usr/bin/python2.3
from __future__ import generators
import sys, re, copy, cStringIO, os
from sets import Set

versionString = '0.01'

autoCompilerFlags = [] # filled based on platform (with no connection to what is in the xcode project)
compilerTranslations = {} # dictionary used to convert options recognized by mac's g++ (used by xcode) to the appropriate flag for the targetted compiler/platform
headersMatch = re.compile(r'.+\.h(pp)?$')
sourceMatch = re.compile(r'.+\.c(pp)?$')
			
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
#	
#	xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
# 	xsi:noNamespaceSchemaLocation="file:/Users/mholder/Documents/projects/Phycas/phycasdev/build/autoVC/phycas/vcproj.xsd"

vcTemplateStrings = [
'''<?xml version="1.0" encoding="Windows-1252"?>
<VisualStudioProject
 ProjectType="Visual C++"
 Version="7.10"
 Name="phycas"
 ProjectGUID="{A3BCC73B-EBCD-49E2-80C1-B4797A427BF5}"
 RootNamespace="phycas"
 Keyword="Win32Proj">
	<Platforms>
		<Platform
			Name="Win32"/>
	</Platforms>
	<Configurations>
''',
'''	</Configurations>
	<References>
	</References>
	<Files>
''',
'''	</Files>
	<Globals>
	</Globals>
</VisualStudioProject>
''']
vcCommonConfigEnd = '''<Tool
				Name="VCMIDLTool"/>
			<Tool
				Name="VCPostBuildEventTool"/>
			<Tool
				Name="VCPreBuildEventTool"/>
			<Tool
				Name="VCPreLinkEventTool"/>
			<Tool
				Name="VCResourceCompilerTool"/>
			<Tool
				Name="VCWebServiceProxyGeneratorTool"/>
			<Tool
				Name="VCXMLDataGeneratorTool"/>
			<Tool
				Name="VCWebDeploymentTool"/>
			<Tool
				Name="VCManagedWrapperGeneratorTool"/>
			<Tool
				Name="VCAuxiliaryManagedWrapperGeneratorTool"/>
		</Configuration>'''
vcResourceFileFilter = '''<Filter
			Name="Resource Files"
			Filter="rc;ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe;resx"
			UniqueIdentifier="{67DA6AB6-F800-4c08-8B7A-83BB121AAD01}">
		</Filter>'''
		
class PBXBuildFile(XCodeDOMObject): pass
class PBXBuildStyle(XCodeDOMObject):
	def standardizeNames(self):
		if self.name == 'Development': self.name = 'Debug'
		elif self.name == 'Deployment': self.name = 'Release'
class PBXProject(XCodeDOMObject):
	'''the root object of an xcode project file.'''
	def writeVCProj(self, out, allFiles):
		global vcPhycasRoot, vcBoostRoot, vcTemplateStrings, vcCommonConfigEnd, vcResourceFileFilter
		out.write(vcTemplateStrings[0])
		for t in self.targets:
			t.buildSettings.writeVCConfiguration(out)
		out.write(vcTemplateStrings[1])
		self.writeFileFilters(out, allFiles)
		out.write(vcTemplateStrings[2])
	def expandTargets(self, buildStyleList):
		et = []
		for t in self.targets:
			et.extend(t.createResolvedCopies(buildStyleList))
		self.targets = et
	def writeGroupFileFilters(self, out, allFiles, s, isHeader = None):
		if isHeader == None:
			self.writeGroupFileFilters(out, allFiles, s, False)
			self.writeGroupFileFilters(out, allFiles, s, True)
			return
		global sourceMatch, headersMatch
		if isHeader:
			m = headersMatch
			filterStr = 'h;hpp;hxx;hm;inl;inc;xsd'
			groupName = s.capitalize() + ' Header'
		else:
			m = sourceMatch
			filterStr = 'cpp;c;cxx;def;odl;idl;hpj;bat;asm;asmx'
			groupName = s.capitalize() + ' Source'
		out.write('\t\t<Filter\n\t\t\tName="%s Files"\n\t\t\tFilter="%s">\n' % (groupName, filterStr))
		allFiles.writeVCFiles(out, s, m, self.targets)
		out.write('\t\t</Filter>\n')
		
	def writeFileFilters(self, out, allFiles):
		self.writeGroupFileFilters(out, allFiles, 'phycas')
		out.write('\t\t<Filter\n\t\tName="Resource Files"\n\t\tFilter="rc;ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe;resx">\n\t\t</Filter>\n')
		self.writeGroupFileFilters(out, allFiles, 'ncl')
		self.writeGroupFileFilters(out, allFiles, 'gui')
		
	
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
def makePathAbsolute(p): #searches the global sourceTreeAlias dictionary and replaces a path prefix that matches a source tree name
	global sourceTreeAlias
	for k, v in sourceTreeAlias.iteritems():
		if p.startswith(k):
			return v + p[len(k):]
		if p.startswith(k[1:]):			#in some context the file omits the $
			return v + p[len(k) - 1:]
	return p

dirSepPattern = re.compile(r'/')
def makeVCPath(p):
	global dirSepPattern
	vcReplacements = { 	'$PHYCAS_ROOT' : r'$(SolutionDir)..\..',		# these shouldn't be hard coded
						'PHYCAS_ROOT' : r'$(SolutionDir)..\..',
						'$FULL_NCL_ROOT' : r'$(SolutionDir)..\..',
						'FULL_NCL_ROOT' : r'$(SolutionDir)..\..',
						'$BOOST_ROOT' : r'$(SolutionDir)..\..\..\boost-1.30.2',
						'BOOST_ROOT' : r'$(SolutionDir)..\..\..\boost-1.30.2'}
	for k,v in vcReplacements.iteritems():
		if p.startswith(k) : 
			p = v + p[len(k):]
			break
	return dirSepPattern.sub(r'\\', p) #note we double escape \, because it will be printed as a string
	
def makeTempDirPath(p, dirName):
	global sourceTreeAlias
	for k in sourceTreeAlias.iterkeys():
		if p.startswith(k):
			return './temp/' + dirName + '/' + '_' + p[1:]
		if p.startswith(k[1:]):			#in some context the file omits the $
			return './temp/' + dirName + '/' + '_' + p
	return './temp/' + dirName + '/' + '_' + p

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
		This class merges the two levels to make a target-specific set of flags'''
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
				
	def writeFlags(self, out):
		out.write('# Compiler Flags\n')
		if self.headerSearchPaths != []:
			out.write('INCLUDE_FLAGS =')
			for h in self.headerSearchPaths:
				nonQuoteIndex = h[0] == '"' and 2 or 1
				p = os.environ.get(h[nonQuoteIndex:])
				if p == None: 
					s = 'Unknown source tree (%s not defined)' % h; 
					raise ValueError, s
				sourceTreeAlias[h] = p
				out.write(' -I%s' % p)
			out.write('\n')
		out.write('INCLUDE_FILE =')
		if self.GCC_PREFIX_HEADER != '':
			out.write(' -include %s' % makePathAbsolute(self.GCC_PREFIX_HEADER))
		out.write('\nCPPFLAGS = %s -0%d\n' % (makePathAbsolute(self.OTHER_CPLUSPLUSFLAGS), self.optimizationLevel))
		out.write('\nWARNING_CFLAGS = %s\n' % makePathAbsolute(self.WARNING_CFLAGS))

	def writeVCConfiguration(self, out):
		global vcCommonConfigEnd
		uniqueName = self.PRODUCT_NAME + self.buildStyleName
		out.write('\t\t\t<Configuration\n\t\t\tName="%s|Win32"\n\t\t\tOutputDirectory="%s"\n\t\t\tIntermediateDirectory="%s"\n\t\t\tConfigurationType="1"\n\t\t\tCharacterSet="2">\n\t\t\t' % (uniqueName, uniqueName, uniqueName))
		out.write('<Tool\n\t\t\t\tName="VCCLCompilerTool"\n\t\t\t\tOptimization="%d"\n\t\t\t\t' % self.optimizationLevel)
		out.write('AdditionalIncludeDirectories="%s"\n\t\t\t\t' % self.getVCIncludeDirectories())
		out.write('PreprocessorDefinitions="%s"\n\t\t\t\t' % self.getVCPreprocessorDefs())
		if self.isDebugTarget:
			out.write('MinimalRebuild="TRUE"\n\t\t\t\tBasicRuntimeChecks="3"\n\t\t\t\t')
		out.write('RuntimeLibrary="%d"\n\t\t\t\tRuntimeTypeInfo="TRUE"\n\t\t\t\tUsePrecompiledHeader="0"\n\t\t\t\tWarningLevel="3"\n\t\t\t\t' % (self.isDebugTarget and 1 or 0))
		out.write('Detect64BitPortabilityProblems="TRUE"\n\t\t\t\tDebugInformationFormat="%d"\n\t\t\t\t' % (self.isDebugTarget and 4 or 3))
		out.write('ForcedIncludeFiles="%s"/>\n\t\t\t' % self.getVCForceInclude())
		out.write('<Tool\n\t\t\t\tName="VCCustomBuildTool"/>\n\t\t\t')
		out.write('<Tool\n\t\t\t\tName="VCLinkerTool"\n\t\t\t\tAdditionalDependencies="ws2_32.lib"\n\t\t\t\tOutputFile="$(OutDir)/%s.exe"\n\t\t\t\t' % uniqueName)
		out.write('LinkIncremental="%d"\n\t\t\t\t' % (self.isDebugTarget and 2 or 1))
		out.write('GenerateDebugInformation="%s"\n\t\t\t\t' % ((True or self.isDebugTarget) and 'TRUE' or 'FALSE')) # note we are always creating debug symbols
		if self.isDebugTarget:
			out.write('ProgramDatabaseFile="$(OutDir)/%s.pdb"\n\t\t\t\t' % uniqueName)
		out.write('SubSystem="1"\n\t\t\t\t')
		if not self.isDebugTarget:
			out.write('OptimizeReferences="2"\n\t\t\t\tEnableCOMDATFolding="2"\n\t\t\t\t')
		out.write('TargetMachine="1"/>\n\t\t\t')
		out.write(vcCommonConfigEnd)
		out.write('\n')

	def getVCIncludeDirectories(self):
		l = []
		for incD in self.headerSearchPaths:
			l.append('&quot;' + makeVCPath(incD) + '&quot;')
		return ';'.join(l)
	def getVCPreprocessorDefs(self):
		return 'WIN32;%s;_CONSOLE' % (self.isDebugTarget and '_DEBUG' or 'NDEBUG')
	def getVCForceInclude(self):
		return makeVCPath(self.GCC_PREFIX_HEADER)

def verifyDirs(p):
	if not os.path.exists(p):
		os.makedirs(p)
def findCommonParent(l):
	cp = os.path.commonprefix(l)
	lastInd = cp.rfind('/')
	if lastInd < 0: return ''
	return cp[:lastInd+1]
	
class FileClassifier:
	def __init__(self):
		self.children = {}
	def append(self, path):
		splitpath = path.split('/')
		if len(splitpath) == 1:
			if self.children.get('__hidden_file__') == None:
				self.children['__hidden_file__'] = Set()
			self.children['__hidden_file__'].add(splitpath[0])
		else:
			if self.children.get(splitpath[0]) == None:
				self.children[splitpath[0]] = FileClassifier()
			toAdd = '/'.join(splitpath[1:])
			self.children[splitpath[0]].append(toAdd)
	def getFiles(self, fileMatch, dirSep = '/', prefix = ''):
		ret = []
		if len(prefix) > 0:
			prefix = prefix + dirSep			
		for k, v in self.children.iteritems():
			if k == '__hidden_file__':
				for f in v:
					if fileMatch.match(f):
						ret.append(prefix + f)
			else:
				ret.extend(v.getFiles(fileMatch, dirSep, prefix + k))
		return ret
	def writeVCFiles(self, out, child, fileMatch, targets):
		fList = self.getFiles(fileMatch)
		fileSet = Set()
		for f in fList:
			splitF = f.split('/')
			if splitF[0] == child or (len(splitF) > 1 and splitF[1] == child):
				fileSet.add(f)
		l = [f for f in fileSet]
		l.sort()
		commonPrefix = findCommonParent(l)
		self.writeHierarchichalVCFiles(out, targets, l, commonPrefix)
	def writeHierarchichalVCFiles(self, out, targets, listOfFiles, commonPrefix, indentLevel = '\t\t\t', createAllDirs = False):
		cpi = len(commonPrefix)
		fileChildren = []
		dirChildren = {}
		for f in listOfFiles:
			p = f[cpi:].split('/')
			if len(p) == 1: 
				fileChildren.append(f)
			elif dirChildren.get(p[0]) == None:
				dirChildren[p[0]] = [f]
			else:
				dirChildren[p[0]].append(f)
		for f in fileChildren:
			self.writeFlatVCFile(out, f, targets, indentLevel)
		if createAllDirs or len(dirChildren) > 1 or len(fileChildren) > 1:
			newIndentLevel = indentLevel + '\t'
			keys = [k for k in dirChildren.iterkeys()]
			keys.sort()
			for k in keys:
				out.write('%s<Filter\n%sName="%s"\n%sFilter="">\n' % (indentLevel, newIndentLevel, k, newIndentLevel))
				self.writeHierarchichalVCFiles(out, targets, dirChildren[k], commonPrefix + k + '/', newIndentLevel, True)
				out.write('%s</Filter>\n' % indentLevel)
		
	def writeFlatVCFile(self, out, f, targets, indentLevel):
		out.write('%s<File\n%s\tRelativePath="%s">\n' % (indentLevel, indentLevel, makeVCPath(f)))
		for t in targets:
			t.writeVCFileConfiguration(out, f)
		out.write('%s</File>\n' % indentLevel)

class PBXNativeTarget(XCodeDOMObject):
	def writeVCFileConfiguration(self, out, f):
		if f.endswith('phycas/taxa/taxa_manager.cpp'):
			out.write('\t\t\t\t<FileConfiguration\n\t\t\t\t\tName="%s|Win32">\n\t\t\t\t\t' % (self.buildSettings.PRODUCT_NAME + self.buildSettings.buildStyleName))
			out.write('<Tool\n\t\t\t\t\tName="VCCLCompilerTool"\n\t\t\t\t\tObjectFile="$(IntDir)/$(InputName)1.obj"/>\n\t\t\t\t</FileConfiguration>\n')
					
	def createResolvedCopies(self, listOfBuildStyles):
		resCopies = []
		for bs in listOfBuildStyles:
			c = copy.copy(self)
			c.buildSettings = BuildSettings(bs, c.buildSettings) # note the higher level (non-target specific) build settings have precedence
			resCopies.append(c)
		return resCopies
	def consolidateFiles(self, fileClassifier):
		for bp in self.buildPhases:
			if isinstance(bp, PBXSourcesBuildPhase) or isinstance(bp, PBXHeadersBuildPhase):
				for f in bp.files:
					fileClassifier.append(f.fileRef.path)
		
	def writeMakefile(self, out):
		global versionString, xcodeFileName, domTree, makefileName
		buildName = self.buildSettings.PRODUCT_NAME
		out.write('# product = %s, buildStyle = %s\n' % (buildName, self.buildSettings.buildStyleName))
		out.write('# written by %s version %s, from the XCode project: %s (which has objectVersion = %s)\n' % (sys.argv[0], versionString, xcodeFileName, domTree.objectVersion)) 
		out.write('CXX = g++\n')
		self.buildSettings.writeFlags(out)
		out.write('%s_OBJS = ' % buildName.upper())
		objectsToBuild = {}
		for bp in self.buildPhases:
			if isinstance(bp, PBXSourcesBuildPhase):
				for f in bp.files:
					fp = makeTempDirPath(f.fileRef.path, buildName + self.buildSettings.buildStyleName)
					ind = fp.rfind('.')
					if ind < 0: raise ValueError, ('Expecting . in filename: %s' %fp)
					oFile = fp[:ind] + '.o'
					objectsToBuild[oFile] = makePathAbsolute(f.fileRef.path)
					out.write('\\\n\t%s ' % oFile)
					ind = fp.rfind('/')
					if ind >= 0:
						verifyDirs(fp[:ind + 1])
		out.write('\nclean: \n')
		for o, c in objectsToBuild.iteritems():
			out.write('\trm -f %s\n' % o)
		out.write('\n# now write rules for each file (crude, but it is autogenerated so why not?)\n')
		out.write('%s : %s\n' % (buildName, makefileName))
		for o, c in objectsToBuild.iteritems():
			out.write('\t$(CXX) $(INCLUDE_FILE) $(INCLUDE_FLAGS) $(CPPFLAGS) $(WARNING_CFLAGS) -c %s -o %s \n' % (c, o))
		out.write('\t$(CXX) -o %s $(%s_OBJS) -lm -lstdc++' % (buildName, buildName.upper()))
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

def setPlatformSpecificInfo(targetOS): #currently assumes g++ (Mac or Linux is the only question)
	#print 'targetting ', targetOS
	global autoCompilerFlags, compilerTranslations
	if targetOS == 'LINUX':
		autoCompilerFlags.append('-DLINUX_PHOREST')
		compilerTranslations['-Wno-long-double'] = ''
		compilerTranslations['-Wmost'] = ''
		compilerTranslations['-Wno-four-char-constants'] = ''

	
if __name__ == '__main__':
	if len(sys.argv) < 2:	print 'Usage: %s xcode_project {Makefile_name}' % sys.argv[0]; sys.exit(1)
	targettedOS  = len(sys.argv) > 2 and sys.argv[2] or 'LINUX'
	targettedOS = targettedOS.upper()
	setPlatformSpecificInfo(targettedOS)
	xcodeFileName = sys.argv[1]
	nestingLevel = -1 #file starts with { which will take us to nestingLevel 0
	tokens = [i for i in tokenizeFile(xcodeFileName)]
	tokenIter = iter(tokens)
	domTree = readDOMObject(tokenIter)
	domTree.resolveReferences()
	domTree.resolvePaths('')
	if False:
		print domTree
	else:
		if False:
			print domTree.rootObject
		else:
			for bs in knownBuildStyles: bs.standardizeNames()
			domTree.rootObject.expandTargets(knownBuildStyles)
			allFiles = FileClassifier()
			for t in domTree.rootObject.targets:
				t.consolidateFiles(allFiles)
				makefileName = 'make_' + t.buildSettings.PRODUCT_NAME + '_' + t.buildSettings.buildStyleName
				out = open(makefileName, 'w')
				t.writeMakefile(out)
			#print allFiles.getFiles(headersMatch)
			#print allFiles.getFiles(sourceMatch)
			#stringOut = cStringIO.StringIO()
			out = open('autoVC/phycas/phycas.vcproj', 'w')
			domTree.rootObject.writeVCProj(out, allFiles)
			#print stringOut.getvalue()
	sys.exit(0)
