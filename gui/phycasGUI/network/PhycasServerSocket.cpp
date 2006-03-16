#include "./PhycasServerSocket.h"
#include "./PhycasSocketException.h"
#if defined(WIN_PHOREST)
ServerSocket::ServerSocket(int port, unsigned connections, TypeSocket type)
	{
	sockaddr_in sa;
	memset(&sa, 0, sizeof(sa));

	sa.sin_family = PF_INET;             
	sa.sin_port = htons(static_cast<u_short>(port));          
	s_ = socket(AF_INET, SOCK_STREAM, 0);
	if (s_ == INVALID_SOCKET) 
		{
		throw SocketException("INVALID_SOCKET");
		}

	if(type==NonBlockingSocket) 
		{
		u_long arg = 1;
		ioctlsocket(s_, FIONBIO, &arg);
		}

	/* bind the socket to the internet address */
	if (bind(s_, (sockaddr *)&sa, sizeof(sockaddr_in)) == SOCKET_ERROR) 
		{
		closesocket(s_);
		throw SocketException("INVALID_SOCKET");
		}

	listen(s_, static_cast<int>(connections));                               
	}

Socket* ServerSocket::Accept() const {
  SOCKET new_sock = accept(s_, 0, 0);
  if (new_sock == INVALID_SOCKET) {
          int rc = WSAGetLastError();
          if(rc==WSAEWOULDBLOCK) {
                  return 0; // non-blocking call, no request pending
          }
          else {
            throw SocketException("Invalid Socket");
      }
  }



  Socket* r = new Socket(new_sock);
  return r;
}
#else

// http://www.linuxgazette.com/issue74/tougher.html#3.3
// Implementation of the ServerSocket class



ServerSocket::ServerSocket (int port, unsigned /*nConnections*/)
	:_port(port)
	{
	if (!Socket::create())
		throw SocketException ( "Could not create server socket." );
	if (!Socket::bind(port))
		throw SocketException ( "Could not bind to port." );
	if (!Socket::listen())
		throw SocketException ( "Could not listen to socket." );
	}

Socket * ServerSocket::Accept() const
	{
	Socket * sock = new Socket();
	if (!Socket::accept(*sock))
		{
		delete sock;
		throw SocketException ( "Could not accept socket." );
		}
	return sock;
	}
#endif