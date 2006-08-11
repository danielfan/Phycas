#ifndef NCL_NXSTREESBLOCK_H
#define NCL_NXSTREESBLOCK_H

#include "pyphy/src/ncl/nxs_defs.hpp"
#include "pyphy/src/ncl/nxs_basic_manager.hpp"
#include "pyphy/src/ncl/nxs_block.hpp"
#include "pyphy/src/ncl/taxa/nxs_taxa_manager.hpp"
#include "pyphy/src/ncl/trees/newick_verifier.hpp"

class NxsTreesManager;
/*----------------------------------------------------------------------------------------------------------------------
|	Class that reads trees blocks and stores descriptions of the trees
*/
class NxsTreesBlock : public NxsCommandManagerBlock , public NxsIndexInfo, public NewickVerifier
	{
	public :
		NxsTreesBlock(NxsTaxaManager &	taxaMgr, NxsTreesManager &	treesMgr);
	
		//	Accessors
		//
		unsigned			FindIndex(const std::string &label ) const;	// returns UINT_MAX if not found
		unsigned   			FindIndexFromNamesOnly(const std::string &label ) const;
		const std::string & GetUserSuppliedLabel(unsigned i) const;
		unsigned 			GetMaxLabelLength() const;
		unsigned			GetNumTrees() const;
		unsigned   			GetSize() const	{return GetNumTrees();}
		void				AddTree(const FullTreeDescription &d);
		void				Clear();
		CmdResult 			EndEncountered();
		void				ResetCmdMgrNxsBlock();

		const VecOfTreeDescriptions & GetTreeDescriptions() const;
		unsigned 			GetDefaultTree() const;
		void				ReadNewickRepresentationAndAddToTaxa(NxsToken & token);
	private:
		
		void				AddRecognizedCommands();
		bool				ParseTree(NxsToken &);
		bool				ParseTranslate(NxsToken &);
		
		VecOfTreeDescriptions	trees;				/*	full descriptions of all known trees */
		unsigned 				defTree;			/*	index of the default tree */
		friend class SimpleNode;
		friend class NxsBasicListManager;
		
		NxsTaxaManager	  &	taxaMgr;
		NxsTreesManager	  &	treesMgr;
	};

#if defined (NCL_SUPPORT_OLD_NAMES)
	typedef NxsTreesBlock TreesBlock;
#endif

inline const VecOfTreeDescriptions &NxsTreesBlock::GetTreeDescriptions() const
	{
	return trees;
	}
	
inline unsigned NxsTreesBlock::GetDefaultTree() const
	{
	return defTree;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns the number of trees stored
*/
inline unsigned NxsTreesBlock::GetNumTrees() const 
	{
	return (unsigned) trees.size();
	}
	
inline void NxsTreesBlock::ResetCmdMgrNxsBlock()
	{
	trees.clear();
	defTree = UINT_MAX;
	translationTable.clear();
	}

inline void NxsTreesBlock::Clear()
	{
	Reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	returns the index of the tree with the name "label", or -1 if that name is unknown.  (not case sensitive)
*/
inline unsigned NxsTreesBlock::FindIndex(
  const std::string &label) const
	{
	return GetIndexFromNexusNameGeneric(trees, label);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	returns the index of the tree with the name "label", or -1 if that name is unknown.  (not case sensitive)
*/
inline unsigned NxsTreesBlock::FindIndexFromNamesOnly(
  const std::string &label) const
	{
	return GetIndexFromNexusNameOnlyGeneric(trees, label);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the label (name) of tree i
|	Throws IndexInfo::XNotFound() if the index ind is not valid.
*/
inline const std::string & NxsTreesBlock::GetUserSuppliedLabel(unsigned i) const
	{
	return GetLabelFromVecGeneric(trees, i);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the length of the longest tree name.
*/
inline unsigned NxsTreesBlock::GetMaxLabelLength() const
	{
	return (unsigned) GetMaxNexusLabelGeneric(trees);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Adds a new tree description to the vector of trees.
*/
inline void NxsTreesBlock::AddTree(const FullTreeDescription &d)
	{
	trees.push_back(d);
	}

#endif
