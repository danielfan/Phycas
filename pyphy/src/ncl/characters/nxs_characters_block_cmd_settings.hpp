#if !defined(NCL_NXSCHARACTERSBLOCK_CMD_SETTINGS_H)
#define NCL_NXSCHARACTERSBLOCK_CMD_SETTINGS_H

#include "phypy/src/ncl/misc/nxs_index_set.hpp"

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
