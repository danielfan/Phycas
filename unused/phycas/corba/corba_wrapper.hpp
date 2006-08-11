#if ! defined(PHYCAS_CORBA_WRAPPER_HPP)
#define PHYCAS_CORBA_WRAPPER_HPP

#include<vector>
#if defined (CORBA_CLIENT_PHYCAS)
#	include "phycas/corba/interfaces/TreeDrawerC.h"
#	include "phycas/corba/interfaces/TreeDatabaseC.h"
#endif
#undef USING_POA_POLICIES
#include "orbsvcs/orbsvcs/CosNamingC.h"
#include <boost/shared_ptr.hpp>
	// Mesquite is not currently using factories to generate interface objects
#define USING_FACTORY
class PhoFloor;
	/// wrapper around CORBA::ORB_var which calls orb->destroy on destruction (for exception-safety)
namespace cipresCORBA 
{
class CorbaWrapper;
	// wrapper that calls orb->destroy in its destructor for exception safety
class OrbWrapper
	{
	public:
		OrbWrapper(CORBA::ORB_var  o);
		~OrbWrapper();
	protected:
		CORBA::ORB_var orb;  ///
		friend class CorbaWrapper;
	};
typedef boost::shared_ptr<OrbWrapper> OrbWrapperShPtr;

class CorbaWrapper
	{
	public:
#		if defined (CORBA_CLIENT_PHYCAS)
#			if defined(USING_FACTORY)
				typedef TreeDrawFactory_var TreeDrawerRef;
				typedef TreeDBFactory_var TreeDBRef;
#			else
#				error deprecated
#			endif
#		endif
		STATIC_SINGLETON CorbaWrapper * GetInstance(int argc = 0, char * argv[] = NULL);
		CorbaWrapper(int argc, char * argv[]);

		~CorbaWrapper()
			{
			}
#		if defined (CORBA_SERVER_PHYCAS)
			void CorbaWrapper::RunServer(PhoFloor &floor);
#		endif
	protected:
		OrbWrapperShPtr orbWrapperPtr;
		CosNaming::NamingContext_var naming_context;
	public:
#		if defined (CORBA_CLIENT_PHYCAS)
#			if defined(USING_FACTORY)
				TreeDrawerRef treeDrawerFactory;
				TreeDBFactory_var treeDBFactory;
#			else
				StoredRemoteRef treeDrawer;
#			endif
#		endif  //if defined (CORBA_CLIENT_PHYCAS)
	};
} // namespace cipresCORBA

#endif //! defined(PHYCAS_CORBA_WRAPPER_HPP)
