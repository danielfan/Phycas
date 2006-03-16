#!/usr/bin/env python

import sys
from omniORB import CORBA, PortableServer

# Import the stubs for the Naming service
import CosNaming

# Import the stubs and skeletons for the cipresCORBA module
import cipresCORBA, cipresCORBA__POA
from tree import *

# Define an implementation of the Echo interface
class TreeDrawer_i (cipresCORBA__POA.TreeDrawer):
	def __init__(self):
		self.treeN = 0
		self.taxaMgr = None
	def execute(self, msg):
		return msg
	def findOrSetTaxa(self, taxList):
		self.taxaMgr = NexusTaxaManager(taxList)
		return 'bogusTaxContextName'
	def drawTree(self, newickStr):
		self.treeN += 1
		t = Tree(newick = newickStr, taxManager = self.taxaMgr)
		nTempFiles = 1
		ind = self.treeN % nTempFiles
		f = 'temp%d.dot' % ind
		t.show(f)

  # Initialise the ORB
orb = CORBA.ORB_init(sys.argv, CORBA.ORB_ID)
  # Find the root POA
poa = orb.resolve_initial_references("RootPOA")
  # Create an instance of Echo_i
tdInstance = TreeDrawer_i()
  # Create an object reference, and implicitly activate the object
tdRef = tdInstance._this()

  # Obtain a reference to the root naming context
obj         = orb.resolve_initial_references("NameService")
rootContext = obj._narrow(CosNaming.NamingContext)
if not rootContext: sys.exit('Failed to narrow the root naming context')

if False:
	# Bind a context named "test.my_context" to the root context
	name = [CosNaming.NameComponent("test", "my_context")]
	
	try:
		testContext = rootContext.bind_new_context(name)
		print "New test context bound"
		
	except CosNaming.NamingContext.AlreadyBound, ex:
		print "Test context already exists"
		obj = rootContext.resolve(name)
		testContext = obj._narrow(CosNaming.NamingContext)
		if testContext is None:
			print "test.mycontext exists but is not a NamingContext"
			sys.exit(1)

# Bind the Echo object to the test context
# previously had "Object" as the second arg, but this seemed to inhibit finding of the NameService
name = [CosNaming.NameComponent("TreeDrawer", "")]

try:
    rootContext.bind(name, tdRef)
    print "New TreeDrawer object bound"

except CosNaming.NamingContext.AlreadyBound:
    rootContext.rebind(name, tdRef)
    print "TreeDrawer binding already existed -- rebound"

    # Note that is should be sufficient to just call rebind() without
    # calling bind() first. Some Naming service implementations
    # incorrectly raise NotFound if rebind() is called for an unknown
    # name, so we use the two-stage approach above

 # Activate the POA
poaManager = poa._get_the_POAManager()
poaManager.activate()

 # Everything is running now, but if this thread drops out of the end
 # of the file, the process will exit. orb.run() just blocks until the
 # ORB is shut down
orb.run()
