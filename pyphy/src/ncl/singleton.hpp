#if ! defined (NCL_SINGLETON_HPP)
#define NCL_SINGLETON_HPP

#if defined(USE_LOKI_SINGLETON)
#	include "loki/realLoki/Singleton.h"
#else

namespace ncl
{
// very simple singleton implementation
template<typename T>
struct SingletonHolder
	{
	public:
		STATIC_SINGLETON T& Instance()
			{
			STATIC_SINGLETON T instance;
			return instance;
			}
	};

} //namespace phycas    

#endif //USE_LOKI_SINGLETON

#endif // NCL_SINGLETON_HPP
