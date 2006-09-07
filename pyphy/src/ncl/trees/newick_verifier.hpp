#ifndef NCL_NEWICK_VERIFIER_H
#define NCL_NEWICK_VERIFIER_H

#include "phypy/src/ncl/trees/full_tree_description.hpp"

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
