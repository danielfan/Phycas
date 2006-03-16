#if !defined (COMPILE_ASSERT_H)
#define COMPILE_ASSERT_H

////////////////////////////////////////////////////////////////////////////////
/// \file compile_assert.hpp
///	Macro for compile-time assert.
//////////

namespace cipres
{
namespace hidden
{
////////////////////////////////////////////////////////////////////////////////
/// used by #COMPILE_TIME_ASSERT
//////////
template<bool> struct CompileTimeChecker
	{
	CompileTimeChecker(...); //default 
	};
////////////////////////////////////////////////////////////////////////////////
/// used by #COMPILE_TIME_ASSERT
//////////
template<> struct CompileTimeChecker<false>{};
}
}

////////////////////////////////////////////////////////////////////////////////
///	\def COMPILE_TIME_ASSERT(condition, msg)
///	\brief An error-detection macro for asserting that conditions are true at compile time.
///
///	Usage:
///
///		COMPILE_TIME_ASSERT(condition to test, error_clue)
///
///	\note	The error_clue must be one alphanumeric word (no spaces of punctuation).
///	The condition to test must be known at compile time (if not a cryptic message such as "illegal non-type template argument"
///	error will be generated.
///	If the compile time assertion evaluates to false, a message such as "Illegal conversion 
///		from ERROR_error_clue to CompileTimeChecker<false> ..." will be generated.
///
///	Implementation Details:
///
///		cipres::hidden::CompileTimeChecker is a boolean-templated class.  
///		The constructor of cipres::hidden::CompileTimeChecker<true> accepts any type
///		The the only constructor for cipres::hidden::CompileTimeChecker<false> is the default constructor
///		The macro COMPILE_TIME_ASSERT(condition, msg):
///			1	Declares a dummy class ERROR_msg
///			2	checks if it can instantiate CompileTimeChecker<condition> from an ERROR_msg type object.
///		If the condition is true, the construction will succeed, if not the error message will be generated
///		Note that while the CompileTimeChecker exists in cipres::hidden:: namespace, the macro is visible
///		by any file that includes compile_assert.h, and he dummy class ERROR_msg will be added to the global 
///		namespace.
///		The resulting code is not affected by the insertion of COMPILE_TIME_ASSERT because the entire construction
///		is inside a sizeof() so no objects are really instantiated.
///	
///	\author	Andrei Alexandrescu "Modern C++ Design"
///	
//////////
#define COMPILE_TIME_ASSERT(condition, msg) {class ERROR_##msg {}; (void)sizeof(cipres::hidden::CompileTimeChecker<(condition)>(ERROR_##msg()));}

#endif
