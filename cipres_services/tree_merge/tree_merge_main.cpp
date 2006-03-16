#include <string>
#include "phycas/force_include.h"
#include "cipres_services/tree_merge/tree_merge_i.hpp"
#include "CipresCommlib/CipresSimpleServer.h"

using std::string;

CipresSimpleServer::NamedServantPtr createTreeMergeServant(const char *);

CipresSimpleServer::NamedServantPtr createTreeMergeServant(const char *)
	{
	const bool verboseMode = true;
	return CipresSimpleServer::NamedServantPtr(string("PhycasTreeMerge"), new CipresIDL_TreeMerge_i(verboseMode));
	}


int main(int argc, char *argv[])
	{
	CipresSimpleServer::ServantFactoryFunc factoryFunc(&createTreeMergeServant);
	CipresSimpleServer server( 	"TreeMerge", // interface name 
								factoryFunc, // servant factory callback function
								false, 		 // multithreading ?
								true		 // persistent ?
								);
	return server.initializeAndRun(argc, argv);
	}
	
