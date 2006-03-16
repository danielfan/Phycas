#!/usr/bin/python2.3
import os,sys

inFile = '../../xml/all_commands.xml'
transformCmd='java org.apache.xalan.xslt.Process'
executableTransformer = 'python translate/commandsJavaTransform.py ./'
	# xsltTransformations is a list of (ouput_path, stylesheet_path) tuples
xsltTransformations = [	('./PhycasMain.java', 'translate/mainJavaTransform.xsl'),
						('../../xml/PhycasMain.xml', 'translate/mainXMLTransform.xsl'),
						('../../xml/Commands.xml', 'translate/commandsXMLTransform.xsl')]
						
						
for	x in xsltTransformations:
        if os.path.exists(x[0]): os.remove(x[0])
	#s = 'rm -f %s' % x[0]
	#print s
	#if os.system(s): print 'Failed'; sys.exit(1)
	s = '%s -IN %s -XSL %s -OUT %s' %(transformCmd, inFile, x[1], x[0])
	print s
	if os.system(s): print 'Failed'; sys.exit(1)
s = executableTransformer + ' ' + inFile
print s
if os.system(s): print 'Failed'; sys.exit(1)
print "OK"
