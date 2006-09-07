/*	dls.h
|
|	Dave Swofford's standard include file.  It should be included before any any other non-system headers.
|
|	Copyright (c) 2003 by David L. Swofford, Florida State University.
|	All rights reserved.
*/

#ifndef DLS_H_INCLUDED
#define DLS_H_INCLUDED

#define LOCAL static	/* allows grep'ing for names of functions that are declared static */
#define ENTRY			/* allows grepping for external function names */
#define FOREVER ;;		/* makes "for (;;)" more readable */

/*	Declare a "Boolean" type.
|
|	IMPORTANT: In PAUP* and some of my other programs, sizeof(Bool) MUST be same size as sizeof(int), e.g., for
|	"*(int *)p = x" constructs where the code doesn't know the precise type of the pointer p.  Thus, we don't use
|	enum or C99/nonstandard types such as 'bool' (C++) or Boolean (Macintosh).
*/
typedef unsigned int Bool;

#if !defined(FALSE)
#	define FALSE	0
#	define TRUE		1
#endif
#if !defined(OFF)
#	define OFF		0
#	define ON		1
#endif
#if !defined(NO)
#	define NO		0
#	define YES		1
#endif

/*	Macros for global variable definition/declaration:
*/
#ifdef MAIN
#	define EXTERN
#	define INIT(x) = x
#else
#	define EXTERN extern
#	define INIT(x)
#endif

/*	Optionally define a program flag here.  This is the ONLY program-specific modification that should be made to
|	this file.
*/
#if !defined(PAUP)
#	define PAUP
#endif

#endif	/* DLS_H_INCLUDED */
