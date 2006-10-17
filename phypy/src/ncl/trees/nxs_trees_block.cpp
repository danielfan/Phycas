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

//#include "phycas/force_include.h"
#include "phypy/src/ncl/nxs_defs.hpp"
#include "phypy/src/ncl/trees/nxs_trees_block.hpp"
#include "phypy/src/ncl/trees/nxs_trees_manager.hpp"
#include "phypy/src/ncl/nxs_exception.hpp"
#include "phypy/src/ncl/command/nxs_auto_command.hpp"
#include "phypy/src/ncl/output/nxs_output.hpp"
#include "phypy/src/ncl/taxa/nxs_taxa_manager.hpp" 
using std::string;

#define TREAT_TAXA_NAMES_AS_TAXA_NUMBERS_ONLY
/*----------------------------------------------------------------------------------------------------------------------
|	Relies on the the branch lengths and tree description being legal! 
*/
void FullTreeDescription::RemoveEdgeLengths()
	{
	if (newick.empty() || hasEdgeLens == kNoEdgeLengthsFound)
		return;
	string newNewick;
	newNewick.reserve(newick.size());	//too long, but a good first approx.
		{
		NxsToken token(newick);
		for (;;)
			{
			++token;
			if (token.AtEOF())
				break;
			if (token.GetTokenReference() == ':')
				{
				++token;
				NXS_ASSERT(!token.AtEOF());
				}
			else
				newNewick += token.GetTokenReference();
			}
		}
	hasEdgeLens = kNoEdgeLengthsFound;
	newick = newNewick;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Alerts the NxsBasicListManager that another trees block has been read via NewBlockRead.  
*/
CmdResult NxsTreesBlock::EndEncountered()
	{
	return treesMgr.NewBlockRead(this);
	}

void NxsTreesBlock::AddRecognizedCommands()
	{
	typedef NxsExternalCommandParser<NxsTreesBlock>							TreesBlockParser;
	typedef NxsManualCommand<TreesBlockParser, NxsDummyExecuteCommandCallback> TreesBlockNoExecuteParsedCommand;
	typedef boost::shared_ptr<TreesBlockNoExecuteParsedCommand> 				TreesBlockNoExecuteParsedCommandPtr;
	
	TreesBlockParser translateParser = TreesBlockParser(this, &NxsTreesBlock::ParseTranslate);
	TreesBlockNoExecuteParsedCommandPtr translateCommand =  TreesBlockNoExecuteParsedCommandPtr(new TreesBlockNoExecuteParsedCommand("Translate", translateParser, NxsDummyExecuteCommandCallback()));
	AddCommand(translateCommand);
	
	TreesBlockParser treeParser = TreesBlockParser(this, &NxsTreesBlock::ParseTree);
	TreesBlockNoExecuteParsedCommandPtr treeCommand =  TreesBlockNoExecuteParsedCommandPtr(new TreesBlockNoExecuteParsedCommand("Tree", treeParser, NxsDummyExecuteCommandCallback()));
	AddCommand(treeCommand);
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	The Callback function used to parse the TRANSLATE command.  Reads in a new translation table to be used in reading 
|	tree descriptions.  Throws NxsException errors for any illegal command.
*/
bool NxsTreesBlock::ParseTranslate(
  NxsToken &token)
  	{
  	VecString newNames; //used if we are allowing implicit naming of taxa (no/inaccurate taxa block)
  	const unsigned origNtax = taxaMgr.GetSize();
  	++token;
  	translationTable.clear();
	for (;;)
		{
		//	Read the key and make sure it isn't a taxon name
		//
		string key(token.GetTokenReference());
		unsigned n = taxaMgr.FindIndexFromNamesOnly(key);
		if (n != UINT_MAX)
			{
			file_pos pos(token.GetFilePosition()); 
			unsigned line(token.GetFileLine()); 
			unsigned col(token.GetFileColumn());
			++token;
			if (!EqualsCaseInsensitive(key, token.GetTokenReference()))
				{
				string s;
				s << key << " is the name of a taxon, so it cannot be used to stand for " << token.GetTokenReference() << " in the Translate command.";
				throw NxsException(s, pos, line, col);
				}
			}
		else
			{
			//	Read the mapped value and make sure it is valid taxon label
			//
			++token;
			string taxName(token.GetTokenReference());
			n = taxaMgr.FindIndex(taxName);
			if (n == UINT_MAX)
				{
				if (allowImplicitNames)
					{
					n = origNtax + newNames.size();
					newNames.push_back(taxName);
					}
				else
					{
					string ss;
					ss << "Unknown taxon label (" << taxName << ") found in the Translate command.";
					throw NxsException(ss, token);
					}
				}
			NStrCaseInsensitiveEquals compF(key);
			TaxNameTransTable::iterator kIt = translationTable.begin();
			for (; kIt != translationTable.end(); ++kIt)
				{
				if (compF(kIt->first))
					break;
				}
			if (kIt != translationTable.end())
				{
				string sss;
				sss << "The taxon translation key " << key << " occurs more than once in the Translate command.";
				throw NxsException(sss, token);
				}
			//	map the key to the taxon name and the taxon index that it represents
			//
			translationTable[key] =  std::make_pair<string, unsigned>(taxName, n);
			}
		++token;
		if (token.Equals(";" ))
			break;
		token.ThrowIfNot(",");
		++token;
		}
	if (allowImplicitNames && !newNames.empty())
		{
		taxaMgr.AppendTaxa(newNames);
		allowImplicitNames = false; //@disallowing implicit names in trees if there are new names in the translate block.  Seems reasonable (to me [mth], at least)
		}
	return true;
  	}

		
/*----------------------------------------------------------------------------------------------------------------------
|	Creates a NxsTreesBlock reader using the taxa list manager supplied to verify taxon names and the output manager to 
|	display output to the user.  
*/
NxsTreesBlock::NxsTreesBlock(NxsTaxaManager &	inTaxaMgr, NxsTreesManager &	inTreesMgr)
  	:NxsCommandManagerBlock("TREES"),
	NewickVerifier(inTaxaMgr),
  	defTree(UINT_MAX),
	taxaMgr(inTaxaMgr),
	treesMgr(inTreesMgr)
   	{
    InitializeRecognizedCommands();
	Reset();
	}

string FullTreeDescription::BriefReport() const
	{
	string s;
	s << name << (rooted ? " (Rooted" : " (Unrooted");
	if (!hasPolytomies)
		s << ", Bifurcating";
	if (hasEdgeLens == kEdgeLengthsFound)
		s << ", with branch lengths" ;
	s << ')';
	return s;
	}
	

/*----------------------------------------------------------------------------------------------------------------------
|	The Callback function used to parse the TREE.  Throws NxsException errors for any illegal tree description.
*/
bool NxsTreesBlock::ParseTree(
  NxsToken &token)
  	{
 	++token;
	bool isDefault = false;
	if (token.GetTokenReference() == '*')
		{
		isDefault = true;
		++token;
		}
  	//	Read the tree name and verify that it is unique
	//
	string treeName(token.GetTokenReference());
	unsigned n = FindIndex(treeName);
	if (n != UINT_MAX)
		{
		errorMsg << "The name " << treeName << " has already been used in this block of trees.";
		throw NxsException(errorMsg, token);
		}
	++token;
	token.ThrowIfNot("=");
	ReadNewickRepresentationAndAddToTaxa(token);	
	token.ThrowIfNot(";");
	if (isDefault)
		defTree = GetNumTrees() - 1;
	return true;
  	}
  	
void NxsTreesBlock::ReadNewickRepresentationAndAddToTaxa(NxsToken & token)
	{
	FullTreeDescription ftd;
	ReadTreeAndAlertTaxaManager(token, &ftd);
	AddTree(ftd);
	}

