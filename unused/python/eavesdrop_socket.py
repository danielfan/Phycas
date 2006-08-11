#! /usr/bin/python
''' Debugging tool for intercepting all information flow between GUI and C++ version of Phycas
	Takes 4 arguments:
		eavesdrop_socket.py listen-IP listen-portNumber forward-IP forward-port 
	Forwards info from the listen socket (to which the Java GUI of phycas will connect) to
	the forward-socket.
	'''
		
import socket
import sys

def sendLine(sock, line, listeningSock):
	if len(line) == 0 or line[0] == '#': return True
	sock.sendall(line)
	print 'Sent: ', line
	msgTail = ''
	while True:
		response = sock.recv(16384)
		if len(response) == 0: 
			print 'Socket disconnected'
			return False
		print 'Received and returning: ', response
		if len(response) < 13:
			msgTail = msgTail + response
		else:
			msgTail = response[-13:]
		if (len(msgTail) > 6 and response[-7:] == '<idle/>') or (len(msgTail) > 12 and response[-13:] == '</user_query>'):
			print 'End Of Transmission'
			listeningSock.sendall(response)
			return line.upper() != 'QUIT' and line.upper() != 'EXIT'	
		listeningSock.sendall(response)
		
if len(sys.argv) < 5:
	print 'excpeting ', sys.argv[0], ' listen-IP listen-portNumber forward-IP forward-port'
	sys.exit(1)
forward_sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
forward_ipHostTuple = (sys.argv[3], int(sys.argv[4]))
forward_sock.connect(forward_ipHostTuple)
listen_sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
listen_ipHostTuple = (sys.argv[1], int(sys.argv[2]))
listen_sock.bind(listen_ipHostTuple)
listen_sock.listen(5)
print 'connected to ', listen_ipHostTuple
print 'forwarding to ', forward_ipHostTuple
new_listen_sock, address = listen_sock.accept()
print 'new connection from ', address
while True:
	request = new_listen_sock.recv(16384)
	if not request: break
	print 'Heard: ', request
	print 'Forwarding...'
	if not sendLine(forward_sock, request, new_listen_sock): 
		print 'Exiting'
		sock.close()
		sys.exit()
sock.close()
