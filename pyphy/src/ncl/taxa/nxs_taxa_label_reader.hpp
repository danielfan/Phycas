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
