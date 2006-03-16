#!/usr/bin/python
import sys
from corba_util import *


nameContext = initOrbAndGetNameService(sys.argv)
tdFacObjRef = nameContext.getService(cipresCORBA.TreeDrawFactory)
if tdFacObjRef: 
	tdObjRef = tdFacObjRef.build()
else:	
	print 'trying state-less treedrawer'
	tdObjRef = nameContext.getService(cipresCORBA.TreeDrawer)
if not tdObjRef: sys.exit('could not get a tree drawer from a factory or the NameService')

tdObjRef.findOrSetTaxa(['aa', 'bb', 'cc'])
	# draw a Tree
trees = ['(aa,bb,cc)', '(aa,(bb,cc))', '((aa,bb),cc)', '((aa,cc),bb)']
for t in trees:
	tdObjRef.drawTree(t)
	print "The Tree %s should have been drawn." % t
