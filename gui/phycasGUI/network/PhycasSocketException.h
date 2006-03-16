// SocketException class


#ifndef SocketException_class
#define SocketException_class

#include <string>

class SocketException
	{
	public:
		SocketException(std::string s)
			:m_s(s) 
			{}
		std::string description() const
			{
			return m_s;
			}
	private:
		std::string m_s;
	};

#endif
