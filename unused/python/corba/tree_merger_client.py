#!/usr/bin/python
import sys
from omniORB import CORBA
import CosNaming, CipresIDL
import re
from corba_util import *
from cipres_tree import CipresTree
 
  # draw a Tree
cipresTrees = [CipresTree('(7:0.089964,(8:0.050264,9:0.034803):0.009528,16:0.055375)'), 
  CipresTree('(1:0.232837,(((2:0.073950,4:0.049896):0.006671,3:0.084989):0.018558,(7:0.062061,((8:0.042580,9:0.042487):0.014108,16:0.050796):0.026376):0.013668):0.019356,(5:0.098338,6:0.081727):0.021046)'),
  CipresTree('(7:0.092413,(((((8:0.032800,(19:0.039050,20:0.031417):0.004826):0.002089,(12:0.044456,13:0.033808):0.000695):0.000952,(17:0.019621,18:0.018780):0.038687):0.007391,11:0.043775):0.005896,9:0.039206):0.005125,16:0.052926)'),
  CipresTree('(7:0.091431,((8:0.053804,(9:0.030934,10:0.035387):0.003327):0.002629,(14:0.047441,15:0.025483):0.011044):0.004281,16:0.053908)')]

nameContext = initOrbAndGetNameService(sys.argv)
tmObjRef = nameContext.getServiceOrExit(CipresIDL.TreeMerge)

print cipresTrees
corbaTrees = map(CipresTree.makeCorbaTree, cipresTrees)

print tmObjRef.mergeTrees(corbaTrees).m_newick


if False:
	  # initialize the ORB and get the NameService
	orb = CORBA.ORB_init(sys.argv, CORBA.ORB_ID)
	obj         = orb.resolve_initial_references("NameService")
	rootContext = obj._narrow(CosNaming.NamingContext)
	if not rootContext: sys.exit('Failed to narrow the root naming context')
	
	  # get a treeDrawer reference
	  # previously had "Object" as the second arg, but this seemed to inhibit finding of the NameService
	name = [CosNaming.NameComponent("TreeMerge", "")]
	try:
		obj = rootContext.resolve(name)
	except CosNaming.NamingContext.NotFound, ex:
		sys.exit('Name not found')
	tmObjRef= obj._narrow(CipresIDL.TreeMerge)
	if not tmObjRef: sys.exit("Object reference is not an Example::Echo")
