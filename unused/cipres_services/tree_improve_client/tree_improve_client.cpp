#include "phycas/force_include.h"
#include <iostream>
#include <string>
#include "CipresCommlib/ConfigDependentHeaders.h"

#include GENERATED_STUB_HEADER(ReadNexus)
#include GENERATED_STUB_HEADER(TreeImprove)
#include NAMING_SERVICE_STUB_HEADER
#include "boost/shared_ptr.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::string;

class OrbWrapper
	{
	public:
		OrbWrapper(CORBA::ORB_var  o);
		~OrbWrapper();
	
		CORBA::ORB_var orb;  ///
	};
typedef boost::shared_ptr<OrbWrapper> OrbWrapperShPtr;

OrbWrapper::OrbWrapper(CORBA::ORB_var  o) // CORBA::ORB_var appears to  be pointer
	:orb(o)
	{ 
	orb = o;
	}

OrbWrapper::~OrbWrapper()
	{
	orb->destroy ();
	}

int main(int argc, char *argv[])
	{
#	if defined(MONITORING_ALLOCATION) && !defined(NDEBUG)
		CREATE_MEMCHK
		{
#	endif
	try
		{
		string filename;
		OrbWrapperShPtr orbWrapperPtr;
		try {
			cerr << "ORB_init(";
			for (int i = 0 ; i < argc; ++i)
				cerr << argv[i] << ' ';
			cerr << ")\n";
			orbWrapperPtr = OrbWrapperShPtr(new OrbWrapper(CORBA::ORB_init (argc, argv, "cipres_TAO_orb"))); //arbitrary orb name
			if (argc < 2)
				{
				cerr << "Expecting the filename of a NEXUS file with characters and a tree as an argument";
				abort();
				}
			filename = argv[1];
			cout << "Filename = "<< filename << endl;
			}
		catch (const CORBA::Exception & ex)
			{
			cerr << "Could not initialize ORB" << endl;
			abort();
			}
		CosNaming::NamingContext_var naming_context;
		try {
			// contact the NameService
			CORBA::Object_var naming_context_object = orbWrapperPtr->orb->resolve_initial_references ("NameService");
			naming_context = CosNaming::NamingContext::_narrow (naming_context_object.in ());
			cerr << "Got NameService" << endl;
			}
		catch(const CORBA::Exception & ex)
			{
			cerr << "Could not get the NameService" << endl;
			abort();
			}
		CosNaming::Name name (1);
		name.length (1);
		CipresIDL::ReadNexus_var nexusReader;
		try	{
			name[0].id = CORBA::string_dup ("ReadNexus");
			CORBA::Object_var rn_object = naming_context->resolve (name);
			nexusReader = CipresIDL::ReadNexus::_narrow (rn_object.in ());
			}
		catch (...)
			{
			cerr << "Could not get ReadNexus object reference"<<endl;
			abort();
			}
		CipresIDL::TreeImprove_var treeImprover;
		
		try	{
			name[0].id = CORBA::string_dup ("TreeImprove");
			// get a reference to the TreeDrawer service
			CORBA::Object_var ti_object = naming_context->resolve (name);
			treeImprover = CipresIDL::TreeImprove::_narrow (ti_object.in ());
			}
		catch (...)
			{
			cerr << "Could not get TreeImprove object reference"<<endl;
			abort();
			}
		
		::CipresIDL::ReadNexus::NumBlockReadSequence * readBlocks = NULL;
		try
			{
			readBlocks = nexusReader->readNexusFile(filename.c_str(), CipresIDL::ReadNexus::NEXUS_TAXA_BLOCK_BIT | CipresIDL::ReadNexus::NEXUS_TREES_BLOCK_BIT | CipresIDL::ReadNexus::NEXUS_CHARACTERS_BLOCK_BIT);
			}
		catch (const CipresIDL::NexusException &e)
			{
			cerr <<"Nexus Error on line "<< e.lineNumber <<":\n"<< e.errorMsg <<endl;
			}
		if (readBlocks->length() < 3 || (*readBlocks)[CipresIDL::ReadNexus::NEXUS_TREES_BLOCK_INDEX] < 1)
			{
			cerr << "Expecting to find a trees block in " << filename<<endl;
			abort();
			}
		if ((*readBlocks)[CipresIDL::ReadNexus::NEXUS_CHARACTERS_BLOCK_INDEX] < 1)
			{
			cerr << "Expecting to find a characters block in " << filename << endl;
			abort();
			}
		cout << "Getting character matrix from NexusReader" << endl;
		::CipresIDL::DataMatrix * chars = nexusReader->getCharacters(0);
		if (!chars)
			{
			cerr << "No matrix obtained." <<endl;
			abort();
			}
		cout << "Sending characters to TreeImprover" << endl;
		treeImprover->setMatrix(*chars);
		
		cout << "Getting trees from NexusReader" << endl;
		::CipresIDL::TreeSeq * trees = nexusReader->getTrees(0);
		if (!trees || trees->length() < 1)
			{
			cerr << "No trees obtained." <<endl;
			abort();
			}
		cout << "Sending the first  tree to TreeImprover" << endl;
		treeImprover->setTree((*trees)[0]);
		
		cout << "Calling improveTree" << endl;
		::CipresIDL::Tree * returnedTree = treeImprover->improveTree(CosEventChannelAdmin::ProxyPushConsumer::_nil());
		if (!trees)
			{
			cerr << "No trees returned." <<endl;
			abort();
			}
		cout << "Tree returned = "<< returnedTree->m_newick << endl;
		
		}
	catch (...)
		{
		cerr << "Terminating as the result of an uncaught exception" << endl;
		exit(2);
		}
#	if defined(MONITORING_ALLOCATION) && !defined(NDEBUG)
		}
		ofstream memf("memcheck.txt");
		MEMCHK_REPORT(memf)
		memf.close();
#	endif
	return 0;
	}
