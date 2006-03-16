#ifndef PYPHY_XPROBDIST_HPP 
#define PYPHY_XPROBDIST_HPP

#include <string>
#include <exception>

/*----------------------------------------------------------------------------------------------------------------------
|	Exception class for transmitting exceptions originating within the ProbDist module from C++ to Python. The attached
|	message should be meaningful if caught and displayed in a Python environment.
*/
class XProbDist : public std::exception
	{
	public:
							XProbDist()						throw();
							XProbDist(const std::string s)	throw();
		virtual				~XProbDist()					throw();

		const char		*	what() const					throw ();

		std::string			msg;	/**< the error to display to the user */
	}; 

/*----------------------------------------------------------------------------------------------------------------------
|	Default constructor.
*/
inline XProbDist::XProbDist() throw ()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor that sets string assigned to the message delivered by this exception (data member `msg').
*/
inline XProbDist::XProbDist(const std::string s) throw()
  :msg() //@POL Mark, did you add this section? If so, why not just assign msg here?
	{
	try {
		msg = s;
		}
	catch (...)
		{
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Virtual destructor.
*/
inline XProbDist::~XProbDist() throw ()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `msg' is empty, returns string "ProbDist module"; otherwise returns string stored by `msg'. Overrides the virtual
|	function in the base class (std::exception).
*/
inline const char *XProbDist::what () const throw ()
	{
	return msg.empty() ? "ProbDist module" : msg.c_str();
	}

#endif
