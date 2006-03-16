#! /usr/bin/python
''' Takes 2 or 3 arguments:
		server_file_to_socket.py IP portNumber filename
	opens a socket for listening,
	passes each line of the file "filename" through the socket to the Gui program,
	printing the sent and returned characters.
	After the file is read the program enters interactive mode, passing each line typed.'''
		
import socket
import sys

def sendLine(sock,line):
	if len(line) == 0 or line[0] == '#': return True
	sock.sendall(line)
	print 'Sent: ', line
	return True

sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
ipHostTuple = (sys.argv[1], int(sys.argv[2]))
sock.bind(ipHostTuple)
sock.listen(1)
print 'Listening on ', ipHostTuple

try:
	while True:
		newSocket, address = sock.accept()
		print 'Connected from ', address
		receivedData = newSocket.recv(16384)
		print 'Initially Received ', receivedData
		if not receivedData: break
	
		if len(sys.argv) > 3:
			for f in sys.argv[3:]:
				inF = open(f, 'r')
				for line in inF:
					if not sendLine(newSocket,line): 
						print 'Exiting from file'
						newSocket.close()
						sock.close()
						sys.exit()
		print 'Waiting ...'
		receivedData = newSocket.recv(16384)
		print 'Received after file ', receivedData
		if not receivedData: break

		prompt = str(ipHostTuple) + '>'
		while True:
			line = raw_input(prompt) + '\n'
			if not sendLine(newSocket,line):  
				print 'Exiting'
				newSocket.close()
				sock.close()
				sys.exit()
			else:
				print 'Waiting ...'
				receivedData = newSocket.recv(16384)
				print 'Received after input ', receivedData
				if not receivedData: break

finally:
	sock.close()


		
