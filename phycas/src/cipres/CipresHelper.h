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

#ifndef CIPRES_HELPER_H__
#define CIPRES_HELPER_H__
#include <iostream>
#include "ConfigDependentHeaders.h"
#include <CipresIDL/api1/CipresS.h>

std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::DataMatrix &matrix);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::RawMatrix &matrix);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::Characters &characters);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::Tree  &tree);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::TaxonSeqSeq  &taxonSubsets);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::TaxonSeq  &taxonSeq);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::CharStateLookup &matrix);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::StateSet &stateSet);




#endif // CIPRES_HELPER_H__
