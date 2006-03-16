#!/usr/bin/python
from read_xcode import *
def makeRelDir(p):
	global sourceTreeAlias
	for k in sourceTreeAlias.iterkeys():
		if p.startswith(k):
			return p[len(k) + 1:] # + 1 removes the dir separator between alias and subdir
		if p.startswith(k[1:]):
			return p[len(k):]
	return p

class MakeBuildSettings(BuildSettings):
	def writeFlags(self, out):
		global sourceTreeAlias
		for p in sourceTreeAlias.itervalues():
			out.write('\t<include>%s\n' % p)
		flags = self.OTHER_CPLUSPLUSFLAGS.split()
		flags.append('-O' + str(self.optimizationLevel))
		flags.extend(self.WARNING_CFLAGS.split())
		for f in flags:
			if f != None and f != '' and f != []:
				out.write('\t<cxxflags>%s\n' % f)
	def makeSourceTreeDict(self):
		global sourceTreeAlias
		if self.headerSearchPaths != []:
			for h in self.headerSearchPaths:
				nonQuoteIndex = h[0] == '"' and 2 or 1
				p = os.environ.get(h[nonQuoteIndex:])
				if p == None: 
					s = 'Unknown source tree (%s not defined)' % h; 
					raise ValueError, s
				sourceTreeAlias[h] = p
		
class MakeTarget(PBXNativeTarget):
	def writeJamfile(self, out, project):
		global versionString, xcodeFileName, makefileName
		buildName = self.buildSettings.PRODUCT_NAME
		out.write('# written by %s version %s, from the XCode project: %s\n' % (sys.argv[0], versionString, xcodeFileName)) 
		out.write('exe jam_%s\n\t:\n' % buildName)
		self.buildSettings.makeSourceTreeDict()
		for f in self.getSources():
			print >>out, '\t', makeRelDir(f)
		out.write('\t:\n\t<define>FORCE_INCLUDE_SOCKET_CONFIG\n')
		self.buildSettings.writeFlags(out)
		out.write('\t;\n')
		
def getCompilerPlatformSettings(targettedOS):
	targettedOS = targettedOS.upper()
	if targettedOS == 'LINUX':
		return (['-DLINUX_PHOREST'],
				{'-Wno-long-double' : '',
				'-Wmost' : '',
				'-Wno-four-char-constants' : ''})
	return ([], {})
if __name__ == '__main__':
	if len(sys.argv) < 2:	print 'Usage: %s xcode_project {Makefile_name}' % sys.argv[0]; sys.exit(1)
	xcodeFileName = sys.argv[1]
	addCompFlags, compilerTrans = getCompilerPlatformSettings(len(sys.argv) > 2 and sys.argv[2] or 'LINUX')
	typeToTypeMap = {PBXNativeTarget : MakeTarget, BuildSettings : MakeBuildSettings}
	project, files = readXCodeProject(xcodeFileName, addCompFlags, compilerTrans, typeToTypeMap)
	for t in project.targets:
		if t.buildSettings.PRODUCT_NAME == 'phycas':
			out = open('Jamfile', 'w')
			t.writeJamfile(out, project)
			sys.exit(0)
			
