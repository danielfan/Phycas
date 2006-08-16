#if ! defined (CIPRES_ASSERT_HPP)
#define CIPRES_ASSERT_HPP

#ifdef __cplusplus
	/* C++ version */
	/* CIPRES_RELEASE_ASSERT is an assert that is not sensitive to NDEBUG */
#	include "pyphy/src/cipres/compile_assert.hpp"
#	include <iostream>
#	define CIPRES_RELEASE_ASSERT(expression)  \
  	((void) ((expression) ? 0 : __TRIP_CIPRES_RELEASE_ASSERT (expression, __FILE__, __LINE__)))

#	define __TRIP_CIPRES_RELEASE_ASSERT(expression, file, lineno)  \
  	(std::cerr << file << ':' << lineno << ": failed assertion\n",	\
  	 abort (), 0)
  	 
#else //__cplusplus

	/* C version */
#	include <stdio.h>
#	define CIPRES_RELEASE_ASSERT(expression)  \
 	 ((void) ((expression) ? 0 : __TRIP_CIPRES_RELEASE_ASSERT (expression, __FILE__, __LINE__)))

#	define __TRIP_CIPRES_RELEASE_ASSERT(expression, file, lineno)  \
 	 (printf ("%s:%u: failed assertion\n", file, lineno),	\
 	  abort (), 0)

#endif //__cplusplus

#endif //#if ! defined (CIPRES_ASSERT_HPP)
