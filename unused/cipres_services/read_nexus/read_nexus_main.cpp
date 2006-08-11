#include <string>
#include "phycas/force_include.h"
#include "cipres_services/read_nexus/read_nexus_i.hpp"
#include "CipresCommlib/CipresSimpleServer.h"

using std::string;

CipresSimpleServer::NamedServantPtr createReadNexusServant(const char *);

CipresSimpleServer::NamedServantPtr createReadNexusServant(const char *)
	{
	const bool verboseMode = true;
	return CipresSimpleServer::NamedServantPtr(string("NCLReadNexus"), new CipresIDL_ReadNexus_i(verboseMode));
	}


int main(int argc, char *argv[])
	{
	CipresSimpleServer::ServantFactoryFunc factoryFunc(&createReadNexusServant);
	CipresSimpleServer server( 	"ReadNexus", // interface name 
								factoryFunc, // servant factory callback function
								false, 		 // multithreading ?
								true		 // persistent ?
								);
	return server.initializeAndRun(argc, argv);
	}
	
