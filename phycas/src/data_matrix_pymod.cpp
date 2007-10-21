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

#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

#include <boost/python.hpp>

//#include "phycas/force_include.h"
#include "phycas/src/cipres/CipresDataMatrixHelper.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(_DataMatrixBase)
{
	class_<CipresNative::DiscreteMatrix, boost::noncopyable>("DataMatrixBase", no_init)
		.def("getNChar", &CipresNative::DiscreteMatrix::getNChar)
		.def("getNTax", &CipresNative::DiscreteMatrix::getNTax)
		.def("getNStates", &CipresNative::DiscreteMatrix::getNStates)
		.def("getSymbolsList", &CipresNative::DiscreteMatrix::getSymbolsList)
		.def("getStateList", &CipresNative::DiscreteMatrix::getStateList, return_value_policy<copy_const_reference>())
		.def("getStateListPos", &CipresNative::DiscreteMatrix::getStateListPos, return_value_policy<copy_const_reference>())
		.def("getCodedDataMatrix", &CipresNative::DiscreteMatrix::getConstNativeC, return_value_policy<copy_const_reference>())
		.def("getRow", &CipresNative::DiscreteMatrix::getRowAsVector)
		;
}
