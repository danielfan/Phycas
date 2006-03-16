#! /usr/bin/python
''' '''
		
import socket, sys, os, re, pexpect
import pexpect

someWordsRE = re.compile(r'.*\S+.*')

expected_prompt = 'converge>'
cmd_dictionary = {'cumulative' : ['treefile', 'burnin', 'increment', 'outfile'] , 'slide' : ['treefile', 'burnin', 'window', 'outfile'] , 'gnuplot' : ['from', 'to', 'all', 'outfile', 'gnucmds', 'plottype', 'plot', 'plotfile'] , 'compare' : ['outfile', 'burnin'] , 'var' : ['outfile', 'window', 'reps', 'metric', 'plot', 'stats', 'amongchain'] , 'showsplits' : ['outfile', 'burnin'] , 'treefileinfo' : ['treefile'] , 'splitpresence' : ['treefile', 'splittree', 'splitpattern', 'burnin', 'outfile'] , 'tracksplit' : ['splittree', 'splitpattern', 'outfile'] , 'factory' : [''] , 'run' : [''] , 'help' : [''] , 'quit' : ['']}

def handleAvailableQuery(s):
	global cmd_dictionary
	if len(s.split(' ')) == 1:
		return '\n'.join([i for i in cmd_dictionary.iterkeys()])
	c = s.split(' ')[1].lower()
	l = cmd_dictionary.get(c)
	if l != None: return '\n'.join(l)
	return ''

# following code from "Python in a Nutshell" 2003 page 435
sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
ipHostTuple = (sys.argv[1], int(sys.argv[2]))
sock.bind(ipHostTuple)
sock.listen(5)
child = pexpect.spawn(sys.argv[3])
child.expect(expected_prompt)
try:
	if True:
		newSocket, address = sock.accept()
		print 'connection from  ', address
		while True:
			print 'waiting...'
			receivedData = newSocket.recv(8192)
			print 'received:  ', receivedData
			if not receivedData: break
			dlist  = receivedData.split('\n')
			for l in dlist:
				if someWordsRE.match(l):
					if l.split(' ')[0].lower() != 'available;':
						print 'sending to converge: ', l
						child.sendline(l)
						child.expect(expected_prompt)
						toSend = '<out>%s</out><idle/>' % child.before
					else:
						print 'intercepting: ', l
						toSend = '<hidden_query>%s</hidden_query><out></out><idle/>' % handleAvailableQuery(l)
					print 'About to return:  ', toSend
					newSocket.sendall(toSend)
					
		newSocket.close()
		print 'Disconnected from', address
finally:
	sock.close()
