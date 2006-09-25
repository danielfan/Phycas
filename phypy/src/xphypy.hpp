#ifndef PHYCAS_XPYPHY_H 
#define PHYCAS_XPYPHY_H

#include <string>
#include <exception>

// General exception class for transmitting exceptions from C++ to Python.
// The attached message should be meaningful if caught and displayed in
// a Python environment.
//
class XPyPhy : public std::exception
	{
	public:
		XPyPhy() throw() {}
		XPyPhy(const std::string s) throw()
			:msg() 
			{
			try {
				msg = s;
				}
			catch (...)
				{
				}
			}
		virtual ~XPyPhy() throw() {}
		std::string	msg;	/**< the error to display to the user */
		const char * what () const throw ()
			{
			return msg.empty() ? "Unknown exception" : msg.c_str();
			}
	}; 

#endif
