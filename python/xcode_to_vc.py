#!/usr/bin/python
from read_xcode import *
import cStringIO
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
	
class VCProject(PBXProject):
	def writeVCProj(self, out, allFiles):
		global vcPhycasRoot, vcBoostRoot, vcTemplateStrings, vcCommonConfigEnd, vcResourceFileFilter
		out.write(vcTemplateStrings[0])
		for t in self.targets:
			t.buildSettings.writeVCConfiguration(out)
		out.write(vcTemplateStrings[1])
		self.writeFileFilters(out, allFiles)
		out.write(vcTemplateStrings[2])
	
	def writeGroupFileFilters(self, out, fileTree, isHeader = None):
		if isHeader == None:
			self.writeGroupFileFilters(out, fileTree, False)
			self.writeGroupFileFilters(out, fileTree, True)
			return
		global sourceMatch, headersMatch
		if isHeader:
			m = headersMatch
			filterStr = 'h;hpp;hxx;hm;inl;inc;xsd'
			groupName = fileTree.name.capitalize() + ' Header'
		else:
			m = sourceMatch
			filterStr = 'cpp;c;cxx;def;odl;idl;hpj;bat;asm;asmx'
			groupName = fileTree.name.capitalize() + ' Source'
		out.write('\t\t<Filter\n\t\t\tName="%s Files"\n\t\t\tFilter="%s">\n' % (groupName, filterStr))
		fileTree.writeHierarchichalVCFiles(out, self.targets, m)
		out.write('\t\t</Filter>\n')
		
	def writeFileFilters(self, out, allFiles):
		out.write('\t\t<Filter\n\t\tName="Resource Files"\n\t\tFilter="rc;ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe;resx">\n\t\t</Filter>\n')
		for f in allFiles:
			self.writeGroupFileFilters(out, f)
		
class VCBuildSettings(BuildSettings):
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

class VCFiles(SourceFileTree):
	
	def writeHierarchichalVCFiles(self, out, targets, fileMatch, indentLevel = '\t\t\t', createAllDirs = False):
		fileChildren = self.getFileList()
		fileChildren.sort()
		for f in fileChildren:
			if fileMatch.match(f):
				self.writeFlatVCFile(out, self.getPathWithSourceTreeRoot() + '/' + f, targets, indentLevel)
		if createAllDirs or len(self.subdirs) > 1 or len(fileChildren) > 1:
			newIndentLevel = indentLevel + '\t'
			keys = [k for k in self.subdirs.iterkeys()]
			keys.sort()
			for k in keys:
				out.write('%s<Filter\n%sName="%s"\n%sFilter="">\n' % (indentLevel, newIndentLevel, k, newIndentLevel))
				self.subdirs[k].writeHierarchichalVCFiles(out, targets, fileMatch, newIndentLevel, True)
				out.write('%s</Filter>\n' % indentLevel)
		
	def writeFlatVCFile(self, out, f, targets, indentLevel):
		out.write('%s<File\n%s\tRelativePath="%s">\n' % (indentLevel, indentLevel, makeVCPath(f)))
		for t in targets:
			t.writeVCFileConfiguration(out, f)
		out.write('%s</File>\n' % indentLevel)
class VCTarget(PBXNativeTarget):
	def writeVCFileConfiguration(self, out, f):
		if f.endswith('phycas/taxa/taxa_manager.cpp'):
			out.write('\t\t\t\t<FileConfiguration\n\t\t\t\t\tName="%s|Win32">\n\t\t\t\t\t' % (self.buildSettings.PRODUCT_NAME + self.buildSettings.buildStyleName))
			out.write('<Tool\n\t\t\t\t\tName="VCCLCompilerTool"\n\t\t\t\t\tObjectFile="$(IntDir)/$(InputName)1.obj"/>\n\t\t\t\t</FileConfiguration>\n')
					
if __name__ == '__main__':
	if len(sys.argv) < 2:	print 'Usage: %s xcode_project' % sys.argv[0]; sys.exit(1)
	xcodeFileName = sys.argv[1]
	project, allFiles = readXCodeProject(xcodeFileName, [], {}, {	SourceFileTree : VCFiles, 
																	PBXNativeTarget : VCTarget,
																	BuildSettings : VCBuildSettings,
																	PBXProject : VCProject
																	})
	out = cStringIO.StringIO()
	project.writeVCProj(out, allFiles)
	print out.getvalue()
