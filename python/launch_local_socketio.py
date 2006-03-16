#!/usr/bin/python
import os,sys, re, time, socket, pexpect

from os import path
if len(sys.argv) < 2 or len(sys.argv) > 3:
	print 'Usage:', sys.argv[0], '<path to executable to test> <path to input file>?'; 
	sys.exit(1) 
executablePath = sys.argv[1]
if not os.path.exists(executablePath):  
	print 'The executable', executablePath, 'does not exist.'
	sys.exit(2)
if not os.access(executablePath, os.X_OK) or os.path.isdir(executablePath):
	print executablePath, 'is not an executable'
	sys.exit(7)
inFilesPath = ''
if len(sys.argv) == 3:
	inFilesPath = sys.argv[2]
	if not os.path.exists(inFilesPath):
		print 'The input file,', inFilesPath, ', does not exist.'
		sys.exit(5)
	if os.path.isdir(inFilesPath):
		print inFilesDir, 'is a directory, (an input file was expected)'
		sys.exit(6)


for socketNumber in range(4444, 4446):
	try:
		print 'going to try',socketNumber
		free_sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		free_sock.connect(('localhost', socketNumber))
		print socketNumber,'is in use.'
		free_sock.close()
		
	except socket.error, x:
		free_sock.close()
		clientCommand = 'file_to_socket.py localhost ' +str(socketNumber) + inFilesPath
		#	serverCommand = executablePath + ' ' + str(socketNumber) +  ' > /dev/null 2>/dev/null'
		#	print serverCommand
		#	os.spawnv(os.P_NOWAIT, executablePath, [socketNumber])
		#	time.sleep(1)
		serverExpectWrapper = pexpect.spawn(executablePath + ' ' + str(socketNumber))
		ind = serverExpectWrapper.expect([pexpect.EOF, pexpect.TIMEOUT], 2)
		if ind == 1:
			clientExpectWrapper = pexpect.spawn(clientCommand)
			while ind == 1:
				ind = clientExpectWrapper.expect([pexpect.EOF, pexpect.TIMEOUT], 1)
				ind = serverExpectWrapper.expect([pexpect.EOF, pexpect.TIMEOUT], 1)
	del free_sock
print 'could not find an open port'
sys.exit(1)
