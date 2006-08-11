#!/usr/bin/python2.3
'''Script that must be invoked from the base of the Phycas tree to perform
both the GUI and NCL transformations from the description of commands.
Creates a tmp directory with a directory structure that reflects the auto-generated portion of the
phycas source tree.  
Transformations are applied to populate the tmp directory structure.
Directory and file comparisons are used to propogate the files to the real source tree iff the files have been changed (to minimize recompilation).'''

import os,sys
from os.path import join as toPath
if os.environ.get('PHYCAS_ROOT') is not None:
	sys.path.append(os.path.join(os.environ['PHYCAS_ROOT'], 'python'))
from phycas_util import quoteDirIfNeeded
def removeFiles(d, excludeList = []):
	'''Removes all files in d except f[-6:]!='readme' and  (but does not affect subdirectories in d).
	Specifically designed for cleaning out oldfiles from the gui/phycasGUI/swixml handler and model directories'''
	lsDir = os.listdir(d)
	for f in lsDir:
		if os.path.isfile(f) and excludeList.count(f) == 0:
			os.remove(f)
def clearOrCreateDir(d):
		if os.path.exists(d):
			print 'clearing ', d
			removeFiles(d)
		else:
			try:
				print 'creating ', d
				os.makedirs(d)
			except OSError:
				sys.exit('Could not create the directory ', d, '\nExiting.')
		
def syncSharedFilesInDirs(sourceDir, destDir):
	print 'comparing %s and %s' %(sourceDir, destDir)
	import filecmp, shutil
	changedFiles = filecmp.dircmp(sourceDir, destDir).diff_files
	changedFiles.extend(filecmp.dircmp(sourceDir, destDir).left_only)
	for f in changedFiles:
		print '\tupdating ', f
		shutil.copy(toPath(sourceDir, f), toPath(destDir, f))
	dirsDiffer = len(changedFiles) > 0
	return dirsDiffer

def syncFile(sourceDir, destDir, filePath):
	dp, fn = os.path.split(filePath)
	fullSourceDir = toPath(sourceDir, dp)
	fullDestDir= toPath(destDir, dp)
	import filecmp, shutil
	dc  = filecmp.dircmp(fullSourceDir, fullDestDir)
	changedFiles = dc.diff_files
	changedFiles.extend(dc.left_only)
	if changedFiles.count(fn) > 0:
		fullSourcePath = toPath(fullSourceDir, fn)
		fullDestPath = toPath(fullDestDir, fn)
		print '\tcopying ', fullSourcePath, ' to ', fullDestPath
		shutil.copy(fullSourcePath, fullDestPath)
		return True
	return False
if __name__ == '__main__':

	doGUITransforms = True
	doNCLTransforms =  True
	badArgMsg = sys.argv[0] + ' accepts only one command line argument {CPP|GUI}. Pass no argument to produce both transformations'
	if len(sys.argv) > 2:
		sys.exit(badArgMsg)
	if len(sys.argv) > 1:
		if sys.argv[1].upper() == 'CPP':
			doGUITransforms = False
		elif sys.argv[1].upper() == 'GUI':
			doNCLTransforms = False
		else:
			sys.exit(sys.argv[1] + ' unrecognized.  ' + badArgMsg)
	try:
		from sax_constructible import SAXConstructible
		necessaryVersion = '0.1'
		if SAXConstructible.version() != necessaryVersion:
			sys.exit('You do not have version %s of the phycasdev/python/sax_constructible.py file.\nYou might need to update from cvs.\nIf you are up-to-date you might need to verify that your PYTHON_PATH variable points to the correct phycasdev/python directory.' % necessaryVersion)
	
	except ImportError:
		sys.exit('you do not have the newest phycasdev/python directory in your PYTHONPATH setting')
	
	xsltTransformCmd = 'java org.apache.xalan.xslt.Process' # executable for xslt transformations
	cwd = os.getcwd()
	guiXMLPath = 	toPath('gui', 'xml')
	cmdXMLFilePath = 	toPath(guiXMLPath, 'all_commands.xml')
	guiSwixmlPath = 	toPath('gui', 'phycasGUI', 'swixml')
	guiTranslatePath = 	toPath(guiSwixmlPath, 'translate') # directory that holds the gui translation scripts
	guiHelpPath = 	toPath('gui', 'help')
	guiHelpHTMLPath = 	toPath(guiHelpPath, 'master')
	guiHelpTranslatePath = 	toPath(guiHelpPath, 'translate') # directory that holds the JavaHelp system translation scripts
	guiPythonScript = 	toPath(guiTranslatePath, 'commandsJavaTransform.py' )	# path to the python script for GUI transformations
	nclTranslatePath = toPath('phycas', 'command') # directory that holds the ncl translation scripts
	nclPythonScript = 	toPath(nclTranslatePath, 'commands_ncl_transform.py') # path to the python script for NCL transformations
	
	guiPythonExecuteCmd = 'python ' + quoteDirIfNeeded(toPath(cwd, guiPythonScript))	# invocation of the python script for GUI transformations
	nclPythonExecuteCmd = 'python ' + quoteDirIfNeeded(toPath(cwd, nclPythonScript)) # invocation of the python script for NCL transformations
	relTempDir= 'tmp'
	tempDirPrefix = toPath(cwd, relTempDir)
	fullPathToXML = quoteDirIfNeeded(toPath(cwd, cmdXMLFilePath))
	
		# xsltTransformations is a list of (output_path, stylesheet_path) tuples
	xsltTransformations = [	(toPath(guiSwixmlPath, 'PhycasMain.java'), 	toPath(guiTranslatePath, 'mainJavaTransform.xsl')),
							(toPath(guiXMLPath, 'PhycasMain.xml'), 		toPath(guiTranslatePath, 'mainXMLTransform.xsl')),
							(toPath(guiXMLPath, 'Commands.xml'), 		toPath(guiTranslatePath, 'commandsXMLTransform.xsl')),
							(toPath(guiHelpPath, 'Master.jhm'), 		toPath(guiHelpTranslatePath, 'mapXMLTransform.xsl')),
							(toPath(guiHelpPath, 'MasterTOC.xml'),		toPath(guiHelpTranslatePath, 'tocXMLTransform.xsl')),
							(toPath(guiHelpPath, 'MasterIndex.xml'), 	toPath(guiHelpTranslatePath, 'indexXMLTransform.xsl')),
							(toPath(guiHelpHTMLPath, 'commands.html'), 	toPath(guiHelpTranslatePath, 'commandsHTMLTransform.xsl')),
						  ]
							
	guiFilesToSync = [i[0] for i in xsltTransformations]				
	nclFilesToSync = []				
	guiDirsToSync= [toPath(guiSwixmlPath, 'handler')]
	nclDirsToSync= [ nclTranslatePath ]
	
	dirsWithSyncedIndivFiles = [os.path.dirname(i) for i in guiFilesToSync]
	mirroredDirs = nclDirsToSync + guiDirsToSync
	dirSkeleton = dirsWithSyncedIndivFiles + mirroredDirs
	for d in dirSkeleton: 
		clearOrCreateDir(toPath(tempDirPrefix , d))
	
	my_cmds = []
	if doGUITransforms:
		for	x in xsltTransformations:
			s = xsltTransformCmd + ' -IN ' + fullPathToXML + ' -XSL ' +  quoteDirIfNeeded(toPath(cwd, x[1])) + ' -OUT ' +  quoteDirIfNeeded(toPath(tempDirPrefix, x[0]))
			my_cmds.append(s)
		my_cmds.append(guiPythonExecuteCmd + ' ' + quoteDirIfNeeded(cwd) + ' ' + fullPathToXML + ' ' + quoteDirIfNeeded(toPath(tempDirPrefix, guiSwixmlPath)))
	if doNCLTransforms:
		my_cmds.append(nclPythonExecuteCmd + ' ' + quoteDirIfNeeded(cwd) + ' ' + quoteDirIfNeeded(toPath(relTempDir, nclTranslatePath)))
	for c in my_cmds:
		print c
		if os.system(c):
			sys.exit('Failure: ' + c)
	print 'New files created in ', tempDirPrefix
	guiFilesTouched = False
	nclFilesTouched = False
	if doGUITransforms:
		for d in guiDirsToSync:
			guiFilesTouched = syncSharedFilesInDirs(toPath(tempDirPrefix, d), toPath(cwd, d)) or guiFilesTouched
		for f in guiFilesToSync:
			guiFilesTouched = syncFile(tempDirPrefix, cwd, f) or guiFilesTouched
	if doNCLTransforms:
		for d in nclDirsToSync:
			nclFilesTouched = syncSharedFilesInDirs(toPath(tempDirPrefix, d), toPath(cwd, d)) or nclFilesTouched
		for f in nclFilesToSync:
			nclFilesTouched = syncFile(tempDirPrefix, cwd, f) or nclFilesTouched
	
	print '\n\n'
	if guiFilesTouched or nclFilesTouched:
		if nclFilesTouched: print 'C++ program MUST be recompiled !!!'
		else: 'C++ program was not touched but...'
		if guiFilesTouched: print 'Java Gui MUST be recompiled !!!'
		else: print ', but the JAVA GUI was not touched.'
	else:
		print 'No recompilation necessary.'
