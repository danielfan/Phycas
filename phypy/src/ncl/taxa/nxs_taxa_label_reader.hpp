/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
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

#ifndef NCL_NXSTAXALABELREADER_H
#define NCL_NXSTAXALABELREADER_H

#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/command/nxs_command_decl.hpp"
#include "phypy/src/ncl/nxs_basic_manager.hpp"
class NxsToken;
/*----------------------------------------------------------------------------------------------------------------------
|	Parses the taxlabels command (which can occur in taxa, characters, data, unaligned, distances, etc blocks)
|	NOTE:	No tests are added to the TaxLabels command by CreateTaxLabelsCommand.  Derived classes should add tests to
|		assure that the dimension command has been read.  When the dimensions are read taxLabels should be filled with
|		numbers for the taxa names.  These numbers are replaced in ParseTaxLabels
*/
class NxsTaxaLabelReader
	{
	protected:
		NxsCommandShPtr CreateTaxLabelsCommand();	
		bool 				ParseTaxLabels(NxsToken &);	
		
		OrderedCaseInsensitiveLabels	taxLabels;		/* known taxon labels */
		bool				labelsRead;		/*should be set to false in derived class's reset function.  Set to true at the end of ParseTaxLabels*/
	};

#endif
