#if !defined(LOKI_TOP_LEVEL_SINGLETON_H)
# define LOKI_TOP_LEVEL_SINGLETON_H
/////////////////////////////////
// Generated header: Singleton.h
// Forwards to the appropriate code
// that works on the detected compiler
// Generated on Sun Sep 15 15:31:17 2002
/////////////////////////////////

#ifdef LOKI_USE_REFERENCE
#	include "Reference/Singleton.h"
#else
#	if (__INTEL_COMPILER)
#		include "Reference/Singleton.h"
#	elif (__MWERKS__)
#		include "Reference/Singleton.h"
#	elif (__BORLANDC__ >= 0x560)
#		include "Borland/Singleton.h"
#	elif (_MSC_VER >= 1300)
#		include "MSVC/1300/Singleton.h"
#	elif (_MSC_VER >= 1200)
#		include "MSVC/1200/Singleton.h"
#	else
#		include "Reference/Singleton.h"
#	endif
#endif

namespace Loki
{
////////////////////////////////////////////////////////////////////////////////
// mth calls Create and Destroy (static functions) in client class's code
// class template CreateUsingCreate
// Implementation of the CreationPolicy used by SingletonHolder
////////////////////////////////////////////////////////////////////////////////

    template <class T> struct CreateUsingCreate
    {
        static T* Create()
        {
        return T::Create();
        }
        
        static void Destroy(T* p)
        {
        T::Destroy(p);
        }
    };

} //namespace Loki
#endif

