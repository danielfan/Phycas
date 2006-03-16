#include "phycas/force_include.h" 
#if defined(CORBA_PHYCAS)

#include "phycas/corba/tree_merge_wrapper.hpp"
#if defined (CORBA_SERVER_PHYCAS)
#	include "phycas/corba/interfaces/tree_merger/TreeMergeI.h"
#	include "CipresFacilitator.h"
#endif
using std::cerr;
using namespace cipresCORBA;

/*
	Todo: need to know whether the tree_merger_impl is thread safe.
	If it is, we can probably put it in the cipres registry as a SHARED_OBJECT
	and change its Lifecycle.remove() impl to do nothing.  If it isn't 
	thread-safe we should create a single-thread poa.
*/
#if defined (CORBA_SERVER_PHYCAS)
	void CorbaTreeMergeWrapper::RunServer()
		{
		try
			{
			// Get the RootPOA and manager
			CORBA::ORB_var orb = CORBA::ORB::_duplicate(FACILITATOR->getORB());
			CORBA::Object_var poa_object =  orb->resolve_initial_references ("RootPOA");
			PortableServer::POA_var poa = PortableServer::POA::_narrow (poa_object.in ());
			PortableServer::POAManager_var poa_manager = poa->the_POAManager ();
				
			poa_manager->activate();
			
			  // Create the servant
			CipresIDL_TreeMerge_i tree_merger_impl;

			 // Activate it to obtain the object reference
			::CipresIDL::TreeMerge_var treeMerger = tree_merger_impl._this ();

			// Facilitator decides whether to use nameservice or write ior to a file 
			// based on the command line args.

			FACILITATOR->publishIOR(treeMerger.in(), "TreeMerge");
			FACILITATOR->run();
			}
		catch (CORBA::Exception &) 
			{
			std::cerr << "CORBA exception raised!" << std::endl;
			}
			FACILITATOR->cleanup();
		}
#endif //CORBA_SERVER_PHYCAS



#endif // defined(CORBA_PHYCAS)
