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

#if !defined(NCL_NXSCHARACTERSBLOCK_CMD_SETTINGS_H)
#define NCL_NXSCHARACTERSBLOCK_CMD_SETTINGS_H

#include "phycas/src/ncl/misc/nxs_index_set.hpp"

class NxsFormatCmdSettings
	{
	public:
		unsigned 	dataTypeIndex;
		std::string	dataTypeName;
		bool		respectingCase;
		char		missingSymbol;
		char		gapSymbol;
		std::string 	rawSymbols;
		std::string 	rawEquate;
		char		matchSymbol;
		bool		labels;
		bool		transposing;
		bool		interleaving;
		bool		tokens;
		unsigned	itemsIndex;
		std::string 	itemTypeName;
		unsigned	stateFormIndex;
		std::string	stateFormName;
		
		NxsFormatCmdSettings()
			:dataTypeIndex(0),
			dataTypeName("STANDARD"),
			respectingCase(false),
			missingSymbol('?'),
			gapSymbol('\0'),
			matchSymbol('\0'),
			labels(true),
			transposing(false),
			interleaving(false),
			tokens(false),
			itemsIndex(0),
			itemTypeName("STATES"),
			stateFormIndex(0),
			stateFormName("STATESPRESENT")
			{}
	};

class NxsEliminateCmdSettings
	{
	public:
		NxsIndexSet 	toEliminate;
	};
	
#endif
