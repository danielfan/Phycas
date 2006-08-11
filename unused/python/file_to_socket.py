#! /usr/bin/python
''' Takes 2 or 3 arguments:
		file_to_socket.py IP(=localhost) portNumber(=4444) filename(= None - go straight to prompt)
	passes each line of the file "filename" through the socket,
	printing the returned characters.
	Waits until '<idle/>' is received before sending the next line.
	After the file is read the program enters interactive mode, passing each line typed.
	'quit' or 'exit' close the socket (after being passed to through the socket).
	'QUITFILETOSOCKET' closes the socket (without passing any final command through)
	takes stdinput to respond to user_query tags'''
		
import socket
import sys

def sendLine(line):
	global sock
	quitFileToSocketCmd = 'QUITFILETOSOCKET'
	fileToSocketOutputPrefix = 'FILE_TO_SOCKET_OUTPUT'
	if line.startswith(quitFileToSocketCmd): 
		sock.close()
		sys.exit(0)
	if line.startswith(fileToSocketOutputPrefix):
		sys.stdout.write(line[len(fileToSocketOutputPrefix):])
		return True
	if len(line) == 0 or line[0] == '#': return True
	sock.sendall(line)
	print >> sys.stderr,  'Sent: ', line
	msgTail = ''
	idleBreaker = '<idle/>'
	lenIdle  = len(idleBreaker)
	userQueryBreaker = '</user_query>'
	lenUserQuery = len(userQueryBreaker)
	lenLongestBreaker = max(lenIdle, lenUserQuery)
	while True:
		response = sock.recv(16384)
		if len(response) == 0: 
			print >> sys.stderr,  'Socket disconnected'
			return False
		sys.stdout.write(response)
		if len(response) < lenLongestBreaker:
			msgTail = msgTail + response
		else:
			msgTail = response[-lenLongestBreaker:]
		if len(msgTail) >= lenIdle and msgTail[-lenIdle:] == idleBreaker:
			print >> sys.stderr,  'End Of Transmission'
			if line.upper().startswith('QUIT') or line.upper().startswith('EXIT'):
				sock.close()
				sys.exit(0)
			return True
		if len(msgTail) >= lenUserQuery and msgTail[-lenUserQuery:] == userQueryBreaker:
			r = raw_input('response to query')
			sock.sendall(r)
	

sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
portN =  len(sys.argv) > 2 and sys.argv[2] or 4444
ipN = len(sys.argv) > 1 and sys.argv[1] or 'localhost'
ipHostTuple = (ipN, int(portN))
if False:
	import time
	maxConnectAttempts = 5 # temp
	for i in range(maxConnectAttempts+1):
		try:
			if i != 0:
				print >> sys.stderr, 'Trying again.'
			sock.connect(ipHostTuple)
			break
		except socket.error, x:
			if i == maxConnectAttempts:
				raise socket.error, x
			print >> sys.stderr,  'Connection failed',x
			time.sleep(1)
else:
	sock.connect(ipHostTuple)
			
print >> sys.stderr,  'connected to ', ipHostTuple
if len(sys.argv) > 3:
	inF = open(sys.argv[3], 'r')
	for line in inF:
		if not sendLine(line): 
			print >> sys.stderr, 'Exiting from file'
if len(sys.argv) < 5 or sys.argv[4] != '-n':
	prompt = str(ipHostTuple) + '>'
	while True:
		line = raw_input(prompt) + '\n'
		if not sendLine(line): 
			print >> sys.stderr, 'Exiting'
			sock.close()
			sys.exit()
	sock.close()
