//	http://www.linuxgazette.com/issue74/tougher.html#3.3
// Definition of the ServerSocket class

#ifndef ServerSocket_class
#define ServerSocket_class

#include "./PhycasSocket.h"


class ServerSocket : public Socket
	{
	public:

# if defined(WIN_PHOREST)
		ServerSocket ( int port, unsigned nConnections = 0, TypeSocket type=BlockingSocket);
# else
		ServerSocket (int port, unsigned nConnections = 0);
		ServerSocket (){}
# endif
		virtual ~ServerSocket() {}

		Socket * Accept () const ;
	private:
		unsigned _port;
	};


#endif
