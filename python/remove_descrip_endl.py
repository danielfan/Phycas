#!/usr/bin/python2.3
import re, sys
if __name__ == '__main__':
	nonWrappedRe = re.compile(r'.+[>]$')
	#for filename in sys.argv[1:]:
	filename = '../gui/xml/all_commands.xml'
	f = open(filename, 'r')
	storedStr = ''
	for l in f:
		noEndl = l.rstrip()
		if nonWrappedRe.match(noEndl):
			if storedStr == '':
				print noEndl
			else:
				print '%s %s' % (storedStr, noEndl.lstrip())
			storedStr = ''
		else:
			if storedStr == '':
				storedStr = noEndl
			else:
				storedStr = storedStr + ' ' + noEndl.lstrip()
			
	sys.exit(0)
