#include "phycas/force_include.h" 
#if defined(CORBA_PHYCAS)

#include "phycas/corba/corba_wrapper.hpp"
#include "orbsvcs/CosNamingC.h"
#include "tao/PortableServer/ORB_Manager.h"
#include "ncl/misc/string_extensions.hpp"
#if defined (CORBA_SERVER_PHYCAS)
#	include "phycas/corba/interfaces/tree_merger/TreeMergeI.h"
#endif
using std::cerr;
using namespace cipresCORBA;

#if defined(USING_POA_POLICIES)
	void destroyPolicies(CORBA::PolicyList & policies);

	void destroyPolicies(CORBA::PolicyList & policies)
		{
		for (CORBA::ULong i = 0; i != policies.length(); ++i)
			policies[i]->destroy ();
		}
#endif// defined(USING_POA_POLICIES)

OrbWrapper::OrbWrapper(CORBA::ORB_var  o) // CORBA::ORB_var appears to  be pointer
	:orb(o)
	{ 
	orb = o;
	}

OrbWrapper::~OrbWrapper()
	{
	orb->destroy ();
	}

CorbaWrapper * CorbaWrapper::GetInstance(int argc, char * argv[])
	{
	STATIC_SINGLETON CorbaWrapper * instance = NULL;
	if (instance == NULL)
		{
		try {
			instance = new CorbaWrapper(argc, argv);
			}
		catch (CORBA::Exception & ex)
			{
    		std::cerr << "Could not create a CorbaWrapper instance" << std::endl;
    		}
		}
	return instance;
	}
	
CorbaWrapper::CorbaWrapper(int argc, char * argv[])
	{
	try {
		std::cerr << "ORB_init(";
		for (int i = 0 ; i < argc; ++i)
			std::cerr << argv[i] << ' ';
		std::cerr << ")\n";
		orbWrapperPtr = OrbWrapperShPtr(new OrbWrapper(CORBA::ORB_init (argc, argv, "cipres_TAO_orb"))); //arbitrary orb name
		}
	catch (CORBA::Exception & ex)
		{
		std::cerr << "Could not initialize ORB" << std::endl;
		throw;
		}
	try {
		// contact the NameService
		CORBA::Object_var naming_context_object = orbWrapperPtr->orb->resolve_initial_references ("NameService");
		naming_context = CosNaming::NamingContext::_narrow (naming_context_object.in ());
		std::cerr << "Got NameService" << std::endl;
		}
	catch(CORBA::Exception & ex)
		{
		std::cerr << "Could not get the NameService" << std::endl;
		throw;
		}
#	if defined (CORBA_CLIENT_PHYCAS)
#		if defined(USING_FACTORY)
			CosNaming::Name name (1);
			name.length (1); //@ 2 seems redundant with ctor
			name[0].id = CORBA::string_dup ("TreeDrawFactory");
				// get a reference to the factory service
			CORBA::Object_var factory_object = naming_context->resolve (name);
			treeDrawerFactory = TreeDrawFactory::_narrow (factory_object.in ());
			
			name[0].id = CORBA::string_dup ("TreeDBFactory");
				// get a reference to the factory service
			factory_object = naming_context->resolve (name);
			treeDBFactory = TreeDBFactory::_narrow (factory_object.in ());
#		else
			CosNaming::Name name (1);
			name.length (1);
			name[0].id = CORBA::string_dup ("TreeDrawer");
				// get a reference to the TreeDrawer service
			CORBA::Object_var td_object = naming_context->resolve (name);
			treeDrawer = TreeDrawer::_narrow (td_object.in ());
#		endif
			// old system (beginning of the Quoter tutorial uses arg
			// to find service 
		//CORBA::Object_var factory_object = orbWrapper.orb->string_to_object (argv[1]);
		//TreeDrawerFactory_var factory = TreeDrawerFactory::_narrow(factory_object.in());
#	endif
	}

#if defined (CORBA_SERVER_PHYCAS)
	void CorbaWrapper::RunServer(PhoFloor &)
		{
		try
			{
			// Get the RootPOA and manager
			CORBA::Object_var poa_object =  orbWrapperPtr->orb->resolve_initial_references ("RootPOA");
			PortableServer::POA_var poa = PortableServer::POA::_narrow (poa_object.in ());
			PortableServer::POAManager_var poa_manager = poa->the_POAManager ();
				
#			if defined (USING_POA_POLICIES)
					//	implementing object creation policies
				CORBA::PolicyList policies(2);
				policies.length(2);
				policies[0] = poa->create_id_assignment_policy(PortableServer::USER_ID);
				policies[1] = poa->create_implicit_activation_policy (PortableServer::NO_IMPLICIT_ACTIVATION);
					// Create a child POA (using the same manager)
				PortableServer::POA_var tree_merger_poa = poa->create_POA ("TreeMerger_POA", poa_manager.in(), policies);
				destroyPolicies(policies); //destroy our copies (the poa has made its own copies)
#			endif //defined (USING_POA_POLICIES)
			
			poa_manager->activate();
			
			      // Create the servant
#			if defined (USING_POA_POLICIES)
				CipresIDL_TreeMerge_i tree_merger_impl(tree_merger_poa.in());
#			else
				CipresIDL_TreeMerge_i tree_merger_impl;
#			endif

			  // Activate it to obtain the object reference
			::CipresIDL::TreeMerge_var treeMerger = tree_merger_impl._this ();
			CosNaming::Name name (1);
			name.length (1);
			name[0].id = CORBA::string_dup ("TreeMerge");
			naming_context->rebind (name, treeMerger.in ());
			std::cerr << "Running as CORBA Server"<< std::endl;
			orbWrapperPtr->orb->run ();
			poa->destroy (1, 1);
			orbWrapperPtr->orb->destroy ();
			}
		catch (CORBA::Exception &) 
			{
			std::cerr << "CORBA exception raised!" << std::endl;
			}
		}
#endif //CORBA_SERVER_PHYCAS


#   if 0 // old constructor code
#		if defined (USING_FACTORY)
			int argc = 1;
			char * argv[] = {"corba_client"};
#		else
			// First initialize the ORB, that will remove some arguments...
			int argc = 3; 
			//@ hard-coding port 1050 and localhost
			std::string nsArg;
			nsArg << "NameService=corbaloc:iiop:localhost:" << nameServicePort << "/NameService";
			boost::scoped_ptr<char> c(strdup(nsArg.c_str()));
			char * argv[] = {"corba_client", "-ORBInitRef", c.get()};
			//POL char * argv[] = {"corba_client", "-ORBInitRef", "NameService=corbaloc:iiop:192.168.1.102:1050/NameService"};
#		endif
#   endif


#endif // defined(CORBA_PHYCAS)
