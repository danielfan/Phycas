#!/usr/bin/python
import os,sys, re
from os import path
from phycas_util import quoteDirIfNeeded
import time
import threading
usingSubProcess = False
if int(sys.version[0]) > 2 or (int(sys.version[0]) == 2 and int(sys.version[2]) >= 4):
	from subprocess import *
	usingSubProcess = True

if len(sys.argv) < 4 :
	print 'Usage:', sys.argv[0], '<path to executable to test> <path to directory containing input files> <path to directory containing files with reference output> <UseFTS>?'; 
	sys.exit(1) 
executablePath = sys.argv[1]
if not os.path.exists(executablePath):  
	print 'The executable', executablePath, 'does not exist.'
	sys.exit(2)
if not os.access(executablePath, os.X_OK) or os.path.isdir(executablePath):
	print executablePath, 'is not an executable'
	sys.exit(7)
inFilesDir = sys.argv[2]
if not os.path.exists(inFilesDir):
	print 'The directory of input commands', inFilesDir, 'does not exist.'
	sys.exit(5)
if not os.path.isdir(inFilesDir):
	print inFilesDir, 'is not a directory (of input files to test).'
	sys.exit(6)
correctOutFilesDir = sys.argv[3]
if not os.path.exists(correctOutFilesDir):
	print 'The directory of expected output', correctOutFilesDir, 'does not exist.'
	sys.exit(3)
if not os.path.isdir(correctOutFilesDir):
	print correctOutFilesDir, 'is not a directory (of expected files with the expected output)'
	sys.exit(4)
testingPhycasThroughSocket = False
if len(sys.argv) > 4 and sys.argv[4].upper() == 'USEFTS':
	testingPhycasThroughSocket = True


if not usingSubProcess:
	assert(sys.argv[0] + 'currently requires python2.4' == '' and False)
verbose = False
def launchServer(launchedEvent, killedEvent):
	#launchedEvent.set()
	#killedEvent.set()
	#return
	global origDir, launchCmd, usingSubProcess
	silentServer = True
	shelcmd = ' '.join([quoteDirIfNeeded(w) for w in launchCmd])
	if not usingSubProcess:
		assert(sys.argv[0] + 'currently requires python2.4' == '' and False)
	if verbose: print '  launching server'
	if silentServer:
		p = Popen(shelcmd, shell=True, stdout = PIPE, stderr = PIPE)
		child_stdout, child_stderr = p.stdout, p.stderr
	else:
		p = Popen(shelcmd, shell=True)
	if verbose: print '  sending serverlaunched event'
	launchedEvent.set()
	if verbose: print '  waiting for server to die'
	if silentServer:
		c  = filter(lambda x: False, child_stdout)
		c  = filter(lambda x: False, child_stderr)
	p.wait()
	if verbose: print '  server done'
	killedEvent.set()
	if verbose: print '  server thread done'
	
def runFTSCmd(p, f, readyToStartEvent, doneEvent):
	global origDir
	clientToFile = True
	forwardFileCmdPrefix = ['python ', os.path.join(origDir, '..', 'python', 'file_to_socket.py'), 'localhost', '4444', p, '-n'] # hard coding path to file to socket
	print 'running ', os.path.basename(p), '...'
	erredirect = open(os.path.join('.', f + '_err'), 'w')
	redirect = open(os.path.join('.', f), 'w')
	if verbose: print '		client waiting for server to start'
	readyToStartEvent.wait()
	shelcmd = ' '.join([quoteDirIfNeeded(w) for w in forwardFileCmdPrefix])
	if verbose: print '		launching client'
	if clientToFile:
		cliProc = Popen(shelcmd, shell=True, stdout = PIPE, stderr=PIPE)
		fo, fe= cliProc.stdout, cliProc.stderr
		if clientToFile:
			redirect.writelines(fo)
			redirect.write('\n')
			redirect.close()
			erredirect.writelines(fe)
			erredirect.write('\n')
			erredirect.close()
	else:
		cliProc = Popen(shelcmd, shell=True)
	cliProc.wait()
	if verbose: print '		client process done'
	if verbose: print '		client wrapper done'
	doneEvent.set()
	if verbose: print '		client thread done'
	
	
def runTestsInDir(inFilesDir):
	readyToGoEvent = threading.Event()
	clientDoneEvent = threading.Event()
	if not os.path.exists(inFilesDir): 
		sys.exit('The directory of input files ', inFilesDir, ' does not exist') 
	inDirListing = os.listdir(inFilesDir)
	for f in inDirListing:
		p = os.path.join(inFilesDir, f)
		if not f.startswith('.') and not os.path.isdir(p):
			if testingPhycasThroughSocket:
				serverLaunchedEvent = threading.Event()
				serverThread = threading.Thread(name = 'server_thread', target = launchServer, args = ([serverLaunchedEvent, readyToGoEvent]))
				clientThread = threading.Thread(name = 'client_thread', target = runFTSCmd, args = ([p, f, serverLaunchedEvent, clientDoneEvent]))
				if verbose: print 'starting client thread'
				clientThread.start()
				if verbose: print 'starting server thread'
				serverThread.start()
			else:
				command = ' '.join(launchCmd_)+ ' ' + p
			if verbose: print 'waiting for client thread to die'
			clientDoneEvent.wait()
			clientDoneEvent.clear()
			if verbose: print 'waiting for server thread to die'
			readyToGoEvent.wait()
			readyToGoEvent.clear()
	os.chdir(origDir)
	os.remove(tempExecutablePath)
	compareTestOutput()

def compareTestOutput():
		#	write the diffs between expected and obtained output to a directory based on the executable name
	global correctOutFilesDir, outFilesDir
	import filecmp
	dc = filecmp.dircmp(correctOutFilesDir, outFilesDir)
	returnCode = 0
	knownDiffsDict = makeKnownDiffDict()
	for f in dc.common_files:
		diffs = [d for d in os.popen('diff ' + os.path.join(correctOutFilesDir, f) + ' ' + os.path.join(outFilesDir, f), 'r')]
		filteredDiffs = filterKnownDiffs(f.strip(), diffs, knownDiffsDict)
		if filteredDiffs != None and len(filteredDiffs.edits):
			returnCode += 1
			print >>sys.stderr, f
			sys.stderr.write(filteredDiffs.getAsDiffFormat())
	if len(dc.left_only) > 0:
		returnCode += 2
		print >>sys.stderr, 'MISSING_FILES'
		for f in dc.left_only:
			print >>sys.stderr, f,
	if len(dc.right_only) > 0:
		returnCode += 4
		print >>sys.stderr, 'EXTRA_FILES'
		for f in dc.right_only:
			print >>sys.stderr, f,
	sys.exit(returnCode)
				
class DiffElement:
	validEditCodeRE= re.compile(r'^\d+(,\d+)?[acd]\d+(,\d+)?$')
	def getStartAndLength(editRangeStr):
		startEnd = editRangeStr.split(',')
		st = int(startEnd[0])
		if len(startEnd) == 1:
			return (st, 1)
		elif len(startEnd) == 2:
			return (st, 1 + int(startEnd[1]) - st)
		m = 'Expecting only one comma in a diff element range string, got ' +  editRangeStr
		raise ValueError, m
	getStartAndLength = staticmethod(getStartAndLength)
	def isValidChangeTag(editRangeStr):
		return DiffElement.validEditCodeRE.match(editRangeStr)
	isValidChangeTag = staticmethod(isValidChangeTag)
	def __init__(self, editPrefix, diffOutputIterable):
		splitEditPrefix = editPrefix.split('c')
		if len(splitEditPrefix) == 2:
			self.editType = 'change'
			self.firstStart , firstLen =  DiffElement.getStartAndLength(splitEditPrefix[0])
			self.firstContent = []
			for i in range(firstLen):
				self.firstContent.append(diffOutputIterable.next()[2:].strip('\r\n')) # slicing of the string skips the '< ' prefix and we also strip away the newline
			divider = diffOutputIterable.next().strip()
			if divider != '---': 
				e = 'Expecting --- in diff format found %s' %divider
				raise ValueError, e
			self.secondStart , secondLen =  DiffElement.getStartAndLength(splitEditPrefix[1])
			self.secondContent = []
			for i in range(secondLen):
				self.secondContent.append(diffOutputIterable.next()[2:].strip('\r\n')) # slicing of the string skips the '< ' prefix and we also strip away the newline
		else:
			splitEditPrefix = editPrefix.split('d')
			if len(splitEditPrefix) == 2:
				self.editType = 'delete'
			else:
				splitEditPrefix = editPrefix.split('a')
				if len(splitEditPrefix) == 2:
					self.editType = 'add'
				else:
					e = 'Expecting edit code with c, d, or a in it found:' % editPrefix
					raise ValueError, e
			self.firstStart , firstLen =  DiffElement.getStartAndLength(splitEditPrefix[0])
			self.secondStart , secondLen =  DiffElement.getStartAndLength(splitEditPrefix[1])
			self.content = []
			lenToUser = self.editType == 'delete' and firstLen or secondLen
			for i in range(lenToUser):
				self.content.append(diffOutputIterable.next()[2:].strip('\r\n')) # slicing of the string skips the '< ' prefix and we also strip away the newline
	def __eq__(self, other):
		if self.editType != other.editType or self.firstStart != other.firstStart or  self.secondStart != other.secondStart:
			return False
		if self.editType == 'change':
			return self.firstContent == other.firstContent and  self.secondContent == other.secondContent
		return self.content == other.content
	def __str__(self):
		if self.editType == 'change':
			return 'Change\n\t' + '\n\t'.join(self.firstContent) + '\nto\n\t' + '\n\t'.join(self.secondContent)
		elif self.editType == 'add':
			return 'Add\n\t' + '\n\t'.join(self.content) + '\nat line %d' % self.firstStart
		return 'Delete\n\t' + '\n\t'.join(self.content) + '\nfrom lines %d-%d' % (self.firstStart, self.first + len(self.content))
	def toDiffFormatRange(s, content):
		if len(content) < 2:
			return str(s)
		return '%d,%d' % (s, len(content) + s - 1)
	toDiffFormatRange = staticmethod(toDiffFormatRange)
	def toDiffContent(content, pref):
		return '\n'.join([pref + c for c in content])
	toDiffContent = staticmethod(toDiffContent)
	def getAsDiffFormat(self):
		if self.editType == 'change':
			s = DiffElement.toDiffFormatRange(self.firstStart, self.firstContent) + self.editType[0] + DiffElement.toDiffFormatRange(self.secondStart, self.secondContent)
			return '%s\n%s\n---\n%s\n' % (s, DiffElement.toDiffContent(self.firstContent, '< '), DiffElement.toDiffContent(self.secondContent, '> '))
		elif self.editType == 'add':
			s = str(self.firstStart) + self.editType[0] + DiffElement.toDiffFormatRange(self.secondStart, self.content)
			return '%s\n%s\n' % (s, DiffElement.toDiffContent(self.content, '> '))
		elif self.editType == 'delete':
			s = DiffElement.toDiffFormatRange(self.firstStart, self.content) + self.editType[0] + str(self.secondStart)
			return '%s\n%s\n' % (s, DiffElement.toDiffContent(self.content, '< '))
		e = 'Invalid DiffElement.editType= %s' % self.editType
		raise ValueError, e

class Diff:
	def __init__(self, diffEls):
		self.edits = diffEls
	def __eq__(self, other):
		return other != None and other.__dict__.has_key('edits') and  self.edits == other.edits
	def __isub__(self, other): 
		if other.__class__ != Diff:
			raise TypeError
		for o in other.edits:
			for i in range(len(self.edits)):
				ed = self.edits[i]
				if ed == o:
					print ed, '\nEQUALS\n', o
					self.edits.pop(i)
					break
		return self
	def __str__(self):return '\n++++++++++++++++++++++++\n'.join([str(e) for e in self.edits])	
	def getAsDiffFormat(self):
		return '\n'.join([e.getAsDiffFormat() for e in self.edits])
			
def	parseNamedConcatenatedDiffOutput(f, diffOutputIterable):
	diffEls = []
	for s in diffOutputIterable:
		if DiffElement.isValidChangeTag(s):
			diffEls.append(DiffElement(s, diffOutputIterable))
		else:
			yield f, Diff(diffEls)
			diffEls = []
			f = s
	yield f, Diff(diffEls)

def makeKnownDiffDict():
	global origDir
	knownDiffsPath = os.path.join(origDir, 'knownDiffs')
	if not os.path.exists(knownDiffsPath):
		return {}
	knownDiffsFile = open(knownDiffsPath)
	knownDiffsIter = iter(knownDiffsFile)
	knownDiffsDict = {}
	try:
		f = knownDiffsIter.next()
		for f, d in parseNamedConcatenatedDiffOutput(f, knownDiffsIter):
			knownDiffsDict[f.strip()] = d
	except StopIteration: pass
	return knownDiffsDict

def filterKnownDiffs(f, d, knownDiffsDict):
	diffsInF = knownDiffsDict.get(f)
	dlist = [readDiff for f, readDiff in parseNamedConcatenatedDiffOutput(f, iter(d))]
	observedDiffs = dlist[0]
	if diffsInF != None:
		observedDiffs -= diffsInF
	return observedDiffs
	
	
	
if __name__ == '__main__':
	executableTested = os.path.basename(executablePath)
	outFilesDir = 'output_' + executableTested
	if not os.path.exists(outFilesDir):
		print 'making outFilesDir'
		os.makedirs(outFilesDir)
	
		# move the executable to test to a standard temp location in the executor directory
	tempExecutablePath = os.path.join(outFilesDir, 'temp.exe')
	import shutil
	shutil.copy2(executablePath, tempExecutablePath)
	if not os.path.exists(tempExecutablePath) or not os.access(tempExecutablePath, os.X_OK):
		print 'Copying', executablePath, 'to', tempExecutablePath, 'failed.'
		sys.exit(8)
	launchCmd = [os.path.join('.','temp.exe')]
	
		# cd to the directory where output should be stored
	origDir = os.getcwd()
	if not os.path.isabs(inFilesDir):
		inFilesDir = os.path.abspath(inFilesDir)
	os.chdir(outFilesDir)
	
	
	testThread = threading.Thread(name = 'testloop', target = runTestsInDir, args = ([inFilesDir]))
	testThread.start()
			#if inFilesDir:
			#	import signal
			#	serverExpectWrapper.kill(signal.SIGTERM)
	
	
