

#ifndef __SOCKET_H__
#define __SOCKET_H__

#if defined (WIN_PHOREST) 
/* 
   Socket.h

   Copyright (C) 2002-2004 René Nyffenegger

   This source code is provided 'as-is', without any express or implied
   warranty. In no event will the author be held liable for any damages
   arising from the use of this software.

   Permission is granted to anyone to use this software for any purpose,
   including commercial applications, and to alter it and redistribute it
   freely, subject to the following restrictions:

   1. The origin of this source code must not be misrepresented; you must not
      claim that you wrote the original source code. If you use this source code
      in a product, an acknowledgment in the product documentation would be
      appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
      misrepresented as being the original source code.

   3. This notice may not be removed or altered from any source distribution.

   René Nyffenegger rene.nyffenegger@adp-gmbh.ch
*/

#include <WinSock2.h>

#include <string>


	enum TypeSocket {BlockingSocket, NonBlockingSocket};

#if 0
	class Socket {
	public:

	  virtual ~Socket();
	  Socket(const Socket&);
	  Socket& operator=(Socket&);

	  std::string receiveLine() const;
	  std::string receiveBytes();

	  void   close();

	  // The parameter of SendBytes is a const reference
	  // because SendBytes does not modify the std::string passed 
	  // (in contrast to SendLine).
	 void   send(const std::string&);

	protected:
	  friend class SocketServer;
	  friend class SocketSelect;

	  Socket(SOCKET s);
	  Socket();


	  SOCKET s_;

	  int* refCounter_;

	private:
	  static void Start();
	  static void End();
	  static int  nofSockets_;
	};

	class SocketClient : public Socket {
	public:
	  SocketClient(const std::string& host, int port);
	};


	class SocketServer : public Socket {
	public:
	  SocketServer(int port, int connections, TypeSocket type=BlockingSocket);

	  Socket* Accept();

	};

	// http://msdn.microsoft.com/library/default.asp?url=/library/en-us/winsock/wsapiref_2tiq.asp
	class SocketSelect {
	  public:
	    SocketSelect(Socket const * const s1, Socket const * const s2=NULL, TypeSocket type=BlockingSocket);

	    bool Readable(Socket const * const s);

	  private:
	    fd_set fds_;
	}; 

#endif //0

#else  // LINUX



// http://www.linuxgazette.com/issue74/tougher.html#3.3
// Definition of the Socket class


#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <string>
#include <arpa/inet.h>


const int MAXHOSTNAME = 200;
const int MAXCONNECTIONS = 5;
const int MAXRECV = 1024;
#endif //WIN_PHOREST

class Socket
	{
	public:
		Socket();
#		if defined(WIN_PHOREST)
	  		Socket(const Socket&);
			Socket(SOCKET s);
	  		Socket& operator=(Socket&);
		  	void   close();
#		else
			bool is_valid() const 
				{
				return m_sock != -1; 
				}
			bool bind(const int port);
			bool listen() const;
			bool accept(Socket&) const;
		private: 
	  		Socket& operator=(Socket&);
	  	public:
#		endif
			
		virtual ~Socket();

		// Server initialization
		bool create();

		// Client initialization
		bool connect(const std::string host, const int port);

		// Data Transimission
		bool send(const std::string &) const;
		int recv(std::string&) const ;


		void set_non_blocking(const bool);
		
		const Socket& operator<<(const std::string& ) const;
		const Socket& operator >> (std::string & s);
		
		std::string ReceiveLine()
			{
			std::string s;
			*this >> s;
			return s;
			}
			
		void SendLine(const std::string &s) const
			{
			*this << s;
			}
		
#		if defined(WIN_PHOREST)
			protected:
			  SOCKET s_;

			  int* refCounter_;

			private:
			  static void Start();
			  static void End();
			  static int  nofSockets_;
#		else
			private:
				int			m_sock;
				sockaddr_in m_addr;
#		endif
	};


#endif 