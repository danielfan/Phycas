#ifndef PHO_SCM_H
#define PHO_SCM_H

#include <fstream>
#include "phycas/trees/split.hpp"
#include "phycas/command/scm_settings.hpp"
#include "phycas/modules/scm/scm_summary.hpp"
class NxsCommandManager;
class PhoTreesManager;
class PhoTaxaManager;
class Tree;
class PhoTreesManager;
class Tree;
typedef boost::shared_ptr<Tree> TreeShPtr;

std::string NewickSCM(const StrVec & toMerge);
TreeShPtr StrictConsensusMerger(const PhoTreesManager & treesMgr);
std::string NewickSCM(const StrVec & toMerge);

class SCM
	{
	public:
							SCM(NxsOutputManager & o, const PhoTreesManager & trm);
		void				SetupSCM(NxsCommandManager *cmdMgr);
		CmdResult			HandleSCM(SCMSettings *s);

	private:
		TreeShPtr			InvokeMergeTrees(const SCMSettings & s);
		STATIC_DATA_FUNC TreeShPtr	MergeTrees(const PhoTreesManager & treesMgr, NxsOutputStream * const outStream = NULL, NxsOutputStreamWrapper * const outFilePtr = NULL);
		NxsOutputManager	  & outputMgr;
		//const PhoTaxaManager  & taxaMgr;
		const PhoTreesManager & treesMgr;
		
		friend TreeShPtr StrictConsensusMerger(const PhoTreesManager & treesMgr);

	};
class XBadTreeInSCM{};
class XNotATaxon: public XBadTreeInSCM {};
class XBadRootTaxon: public XBadTreeInSCM {};

#endif

