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

#ifndef NCL_NEWICK_VERIFIER_H
#define NCL_NEWICK_VERIFIER_H

#include "phycas/src/ncl/trees/full_tree_description.hpp"

class BaseTaxaManager;
class NxsToken;
/*----------------------------------------------------------------------------------------------------------------------
|	Class that reads trees blocks and stores descriptions of the trees
*/
class  NewickVerifier
	{
	public :
		NewickVerifier(BaseTaxaManager & taxaMgr)
		  	:baseTaxaMgr(taxaMgr),
			allowImplicitNames(false)
		  	{
		  	}
   	
	
		void				AllowImplicitNames(bool doIt = true) {allowImplicitNames = doIt;}
		FullTreeDescription ReadNewickTree(NxsToken &token, VecString *newNames = NULL, std::string treeName = std::string() ) const;
		void 				ReadTreeAndAlertTaxaManager(NxsToken &token, FullTreeDescription *ftd);
		unsigned			FindIndexForTaxon(std::string, bool allowNewNames) const;
		
	protected:
		
		typedef std::map<std::string, std::pair<std::string, unsigned> > TaxNameTransTable;
		
		BaseTaxaManager		  & baseTaxaMgr;
		TaxNameTransTable 		translationTable; 	/*	User specified tranlations - mapping of arbitrary string to real taxon information */
		bool					allowImplicitNames;  /* if true, an empty taxon block can be used, and taxa will be added as they are read */
	
		friend class SimpleNode;
	};

#endif
