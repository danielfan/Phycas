#include <iostream>
#include "./PhycasSocket.h"
#include "./PhycasSocketException.h"
#if defined(WIN_PHOREST)
/* 
   Socket.cpp

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




int Socket::nofSockets_= 0;

void Socket::Start() {
  if (!nofSockets_) {
    WSADATA info;
    if (WSAStartup(MAKEWORD(2,0), &info)) {
      throw SocketException("Could not start WSA");
    }
  }
  ++nofSockets_;
}

void Socket::End() {
  WSACleanup();
}

Socket::Socket() : s_(0) {
  Start();
  // UDP: use SOCK_DGRAM instead of SOCK_STREAM
  s_ = socket(AF_INET,SOCK_STREAM,0);

  if (s_ == INVALID_SOCKET) {
    throw SocketException("INVALID_SOCKET");
  }

  refCounter_ = new int(1);
}

Socket::Socket(SOCKET s) : s_(s) {
  Start();
  refCounter_ = new int(1);
};

Socket::~Socket() 
	{
	if (! --(*refCounter_)) 
		{
		close();
		delete refCounter_;
		}

	--nofSockets_;
	if (!nofSockets_) 
		End();
	}

Socket::Socket(const Socket& o) {
  refCounter_=o.refCounter_;
  (*refCounter_)++;
  s_         =o.s_;


  nofSockets_++;
}

Socket& Socket::operator =(Socket& o) {
  (*o.refCounter_)++;


  refCounter_=o.refCounter_;
  s_         =o.s_;


  nofSockets_++;


  return *this;
}

void Socket::close() {
  closesocket(s_);
}

/*
std::string Socket::recv() {
  std::string ret;
  char buf[1024];
 
 for ( ; ; ) {
  u_long arg = 0;
  if (ioctlsocket(s_, FIONREAD, &arg) != 0)
   break;
  if (arg == 0)
   break;
  if (arg > 1024)
   arg = 1024;
  int rv = recv (s_, buf, arg, 0);
  if (rv <= 0)
   break;
  std::string t;
  t.assign (buf, rv);
  ret += t;
 }
 
  return ret;
}

*/
int Socket::recv(std::string & ret) const
	{
	for (;;)
		{
		char r;
		switch(::recv(s_, &r, 1, 0)) 
			{
			case -1:
				if (errno == EAGAIN) 
					return 1;
				// fall through to next case
			case 0: // not connected anymore;
				ret.clear();
				return 0;
			}
		ret += r;
		if (r == '\n')  
			return 1;
		}
	}

/*void Socket::send(std::string s) {
  s << '\n';
  send(s_,s.c_str(),s.length(),0);
}
*/
bool Socket::send(const std::string& s) const
	{
	::send(s_, s.c_str(), (int) s.length(), 0);
	return true;
	}


/*SocketClient::SocketClient(const std::string& host, int port) : Socket() {
  std::string error;


  hostent *he;
  if ((he = gethostbyname(host.c_str())) == 0) {
    error = strerror(errno);
    throw error;
  }


  sockaddr_in addr;
  addr.sin_family = AF_INET;
  addr.sin_port = htons(port);
  addr.sin_addr = *((in_addr *)he->h_addr);
  memset(&(addr.sin_zero), 0, 8); 


  if (::connect(s_, (sockaddr *) &addr, sizeof(sockaddr))) {
    error = strerror(WSAGetLastError());
    throw error;
  }
}

SocketSelect::SocketSelect(Socket const * const s1, Socket const * const s2, TypeSocket type) {
  FD_ZERO(&fds_);
  FD_SET(const_cast<Socket*>(s1)->s_,&fds_);
  if(s2) {
    FD_SET(const_cast<Socket*>(s2)->s_,&fds_);
  }     


  TIMEVAL tval;
  tval.tv_sec  = 0;
  tval.tv_usec = 1;


  TIMEVAL *ptval;
  if(type==NonBlockingSocket) {
    ptval = &tval;
  }
  else { 
    ptval = 0;
  }


  if (select (0, &fds_, (fd_set*) 0, (fd_set*) 0, ptval) 
      == SOCKET_ERROR) throw "Error in select";
}

bool SocketSelect::Readable(Socket const * const s) {
  if (FD_ISSET(s->s_,&fds_)) return true;
  return false;
}
*/
#else //WIN_PHOREST
// http://www.linuxgazette.com/issue74/tougher.html#3.3
// Implementation of the Socket class.


#include "./PhycasSocket.h"
#include "string.h"
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <iostream>
#include "./PhycasSocketException.h"



Socket::Socket() :
m_sock (-1)
	{
	memset (&m_addr, 0, sizeof(m_addr));
	}

Socket::~Socket()
	{
	if (is_valid())
		::close(m_sock);
	}

bool Socket::create()
	{
	m_sock = socket(AF_INET, SOCK_STREAM, 0);
	if (!is_valid())
		return false;
	// TIME_WAIT - argh
	int on = 1;
	return setsockopt(m_sock, SOL_SOCKET, SO_REUSEADDR, (const char*) &on, sizeof(on)) != -1;
	}

bool Socket::bind(const int port)
	{
	if (!is_valid())
		return false;
	m_addr.sin_family = AF_INET;
	m_addr.sin_addr.s_addr = INADDR_ANY;
	m_addr.sin_port = htons (port);
	const int bindReturn = ::bind(m_sock, (struct sockaddr *) &m_addr, sizeof(m_addr));
	return bindReturn != -1;
	}

bool Socket::listen() const
	{
	if (!is_valid())
		return false;
	const int listenReturn = ::listen (m_sock, MAXCONNECTIONS);
	return listenReturn != -1;
	}


bool Socket::accept (Socket& new_socket) const
	{
	int addr_length = sizeof (m_addr);
	new_socket.m_sock = ::accept(m_sock, (sockaddr *) &m_addr, (socklen_t *) &addr_length);

	return (new_socket.m_sock > 0);
	}

#if defined(__APPLE__) //@should be test for bsd sockets not just __apple__
#   define	PHYCAS_NOSIG_SOCKET	SO_NOSIGPIPE
#else
#   define	PHYCAS_NOSIG_SOCKET	MSG_NOSIGNAL
#endif

bool Socket::send (const std::string & s) const
	{
	return -1 != ::send(m_sock, s.c_str(), s.size(), PHYCAS_NOSIG_SOCKET);
	}


int Socket::recv (std::string& s) const
	{
	char buf [ MAXRECV + 1 ];
	s.clear();
	memset(buf, 0, MAXRECV + 1);
	int status = ::recv(m_sock, buf, MAXRECV, 0);
	if (status == -1)
		{
		std::cout << "status == -1   errno == " << errno << "  in Socket::recv\n";
		return 0;
		}
	if (status != 0)
		s = buf;
	return status;
	}


bool Socket::connect (const std::string host, const int port)
	{
	if (!is_valid())
		return false;
	m_addr.sin_family = AF_INET;
	m_addr.sin_port = htons(port);
	int status = inet_pton(AF_INET, host.c_str(), &m_addr.sin_addr);
	return (status != -1 || errno != EAFNOSUPPORT) && (0 == ::connect(m_sock, (sockaddr *) &m_addr, sizeof(m_addr)));
	}

void Socket::set_non_blocking(const bool b)
	{
	int opts = fcntl(m_sock, F_GETFL);
	if (opts >= 0)
		{
		opts = (b ? (opts | O_NONBLOCK): (opts & ~O_NONBLOCK));
		fcntl(m_sock, F_SETFL, opts);
		}
	}
#endif // WIN_PHOREST

const Socket& Socket::operator<<(const std::string & s ) const
	{
	if (!send(s))
		throw SocketException ( "Could not write to socket." );
	return *this;
	}


const Socket& Socket::operator >> ( std::string& s )
	{
	if (!recv(s))
		throw SocketException ( "Could not read from socket." );
	return *this;
	}

