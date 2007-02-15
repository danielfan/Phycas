/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#if !defined (GENERIC_TYPE_MAPPING_HPP)
#define GENERIC_TYPE_MAPPING_HPP

////////////////////////////////////////////////////////////////////////////////
///	Int2Type<compile time constant integer> defines a unique (and stateless) 
///		class associated with a given integer.  Used for compile time dispatching 
///		of function calls or creation of appropriate templated classes.
///
///	defines an unnamed enum "value" that is equal to the integer used to 
///		define the class. 
///	\author	Andrei Alexandrescu "Modern C++ Design"
//////////

template<int v>
class Int2Type
	{
		public:
			enum {value = v};
	};

typedef Int2Type<true>  TrueAsAType;
typedef Int2Type<false> FalseAsAType;

////////////////////////////////////////////////////////////////////////////////
///	Type2Type<typename> defines a unique (and stateless) class for each type
///		that is specified as the template argument.
///	This is useful in controlling the return type of templated functions in
///		lieu of partial template specialization of templated functions (which is 
///		not allowed by the C++ standard)
///	Defines the typedef OriginalType which corresponds to the template argument
///	\author	Andrei Alexandrescu "Modern C++ Design"
//////////

template<typename T>
class Type2Type
	{
		public:
			typedef T OriginalType;
	};

#endif
