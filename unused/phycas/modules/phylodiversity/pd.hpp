#ifndef PHO_PD_H
#define PHO_PD_H

#include <fstream>
#include "phycas/trees/split.hpp"
#include "phycas/command/phylodiversity_settings.hpp"
#include "phycas/modules/phylodiversity/phylodiversity_summary.hpp"
class NxsCommandManager;
class PhoTreesManager;
class PhoTaxaManager;
class Tree;

class Phylodiversity
	{
	public:
		void				SetupPhylodiversity(NxsCommandManager *cmdMgr);
		CmdResult			HandlePhylodiversity(PhylodiversitySettings *s);
		STATELESS_FUNC PhylodiversitySummary ProcessTree(const Split & selSplit, const Tree & treee);
	private:
		Phylodiversity()
		 :taxaMgr(NULL)
		 {}
	public:
		void SetTaxaManager(PhoTaxaManager * inTaxaMgr)
			{
			taxaMgr = inTaxaMgr;
			}
		void SetTreesManager(PhoTreesManager * inTreesMgr)
			{
			treesMgr = inTreesMgr;
			}
		friend PhylodiversityCreator;
		PhoTaxaManager * taxaMgr;
		PhoTreesManager * treesMgr;
	};
	
#endif

