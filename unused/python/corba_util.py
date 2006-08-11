#!/usr/bin/python
import sys
from omniORB import CORBA
import CosNaming

class CorbaError(ValueError):
	def __init__(self, msg = ''):
		self.msg = msg

class NameContextNode:
	def __init__(self, namingContextRef, nameOfContext = 'NameService'):
		self.namingContext = namingContextRef
		self.name = nameOfContext
	def getService(self, narrowToClass, name = ''):
		  # get a treeDrawer reference
		  # previously had "Object" as the second arg, but this seemed to inhibit finding of the NameService
		if name == '':
		 	  # hacking omniorb py seems to use _NP_RepositoryId = 'name/context/path/InterfaceName:Version'
		 	  # name = narrowToClass._NP_RepositoryId.split('/')[-1].split(':')[0]
			  # print name #name = narrowToClass.__class__.__name__
			name = narrowToClass.__name__
		nameComponent = [CosNaming.NameComponent(name, "")]
		try:
			obj = self.namingContext.resolve(nameComponent)
		except CosNaming.NamingContext.NotFound, ex:
		 	return None
		objRef= obj._narrow(narrowToClass)
		if not objRef: 
			raise CorbaError, ("Object reference is not an %s (%s based on __class__.__name__)" % (name, narrowToClass.__class__.__name__))
		return objRef
	def getServiceOrExit(self, narrowToClass, name = ''):
		objRef = self.getService(narrowToClass, name)
		if not objRef: 
			sys.exit('could not get a %s object from the NameService using %s as the Name' % (narrowToClass.__name__, name or narrowToClass.__name__))
		return objRef
def initOrbAndGetNameService(argv):
	 # initialize the ORB and get the NameService
	orb = CORBA.ORB_init(sys.argv, CORBA.ORB_ID)
	obj         = orb.resolve_initial_references("NameService")
	rootContext = obj._narrow(CosNaming.NamingContext)
	if not rootContext: 
		raise CorbaError, 'Failed to narrow the root naming context'
	return NameContextNode(rootContext, 'NameService')

