#include "phycas/force_include.h"
#if defined (READ_PHYCAS_BLOCK)
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "phycas/floor.hpp"
#include "phycas/command/set_settings.hpp"
#include "phycas/command/execute_settings.hpp"
#include "phycas/command/alias_settings.hpp"
#include "ncl/output/temp_basic_output_operators.hpp"
#include "ncl/nxs_token.hpp"
#if defined(SUPPORT_GETTREES)
#   include "phycas/command/gettrees_settings.hpp"
#endif
#include "phycas/trees/tree_inl.hpp"
using ncl::endl;
using std::ifstream;
using std::map;
using std::string;

#if defined(SUPPORT_GETTREES)
	/*----------------------------------------------------------------------------------------------------------------------
	|	Reads Trees block from specified file
	*/
	CmdResult PhoFloor::HandleGetTrees(GetTreesSettings *s)
		{
#   if 0
		ifstream inStream;
		if (!s->fileToRead.Open(inStream))
			return kCmdFailedSilent;
		RestoreOnExit<bool> noub(notifyOfUnknownBlocks); 
		notifyOfUnknownBlocks = false;
		NxsOutputStream &outStream = *outputMgrRef.GetOutputStreamPtr();
		NxsString fn = s->fileToRead.GetFileName();
		fileStack.push(s->fileToRead);
		NxsToken tokenStream(inStream);
		StrVec enabledBlocks;
		enabledBlocks.push_back("TREES");
		if (phoTaxaMgr->GetNumTaxa() == 0)
			enabledBlocks.push_back("TAXA");
		const vector<bool> prevEnabledStatus = DisableAllBlocksExcept(enabledBlocks);
		bool successfullyRead = Execute(tokenStream);
		RestoreEnabledStatus(prevEnabledStatus);
		fileStack.pop();
		if (!successfullyRead)
			return kCmdFailedSilent;
		outStream << "Finished reading trees from "<< fn <<endl;
#   endif 
		return kCmdSucceeded;
		}
	
#endif // defined (SUPPORT_GETTREES)

/*----------------------------------------------------------------------------------------------------------------------
|	This file contains handler functions for some of the very simple Phorest commands (those that are unlikely to change 
|	much and are too small to make the creation of a separate cpp file worth the effort).
*/

CmdResult PhoFloor::HandleSet(SetSettings *s)
	{
#   if defined(NCL_USE_NXS_CONSOLE_OUTPUT)
		NxsOutputStream * outStream = outputMgrRef.GetOutputStreamPtr();
		if (outStream && s->outputWidth != outputMgrRef.GetOutputWidth())
			{
			outputMgrRef.SetOutputWidth(s->outputWidth);
			*outStream << "The output width is now " << s->outputWidth << endl;
			}
#	elif defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(s)
#   endif // NCL_USE_NXS_CONSOLE_OUTPUT
	return kCmdSucceeded;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This file contains handler functions for some of the very simple Phorest commands (those that are unlikely to change 
|	much and are too small to make the creation of a separate cpp file worth the effort).
*/

CmdResult PhoFloor::HandleExecute(ExecuteSettings *settings)
	{
	ifstream inStream;
	if (!settings->fileToExecute.Open(inStream))
		return kCmdFailedSilent; // should not happen if
	const std::string fn = settings->fileToExecute.GetFileName();
	string msg;
	msg << "Executing " << fn;
	NxsOutputOperationStatusID fileOpID = outputMgrRef.StartStatusDisplay(msg, false);
		//
		//	after we begin executing the command we don't want to use s->fileToExecute, because it may be changed (by another
		//  execute command in the file).
	fileStack.push(settings->fileToExecute);
	NxsToken tokenStream(inStream);
	const bool successfullyRead = Execute(tokenStream);
	fileStack.pop();
	if (!successfullyRead)
		return kCmdFailedSilent;
	msg.clear();
	msg << "Finished reading " << fn;
	outputMgrRef.EndStatusDisplay(fileOpID, msg);
	return kCmdSucceeded;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the NxsBlock::stopReading and quitNow to true. 
|	NxsBlock::stopReading stops the reading of the the current token stream
|	quitNow causes exit from the  checked in the infinie loop in  PhoFloor::Run
*/	
CmdResult  PhoFloor::HandleQuit()
	{
	NxsBlock::stopReading = true;
	quitNow = true;
	return kCmdSucceeded;
	}
#if 0
#include "phycas/rand/probability_distribution.hpp"
#include "phycas/trees/tree_node.hpp"
#include "phycas/trees/tree.hpp"
//@POL temp
// better to do this with functor 
class ObsoleteEdgeLenPriorCalculator
	{
	ProbabilityDistribution *probDist;
	double ln_prior;

	public:
		ObsoleteEdgeLenPriorCalculator(ProbabilityDistribution *d)
			{
			probDist = d;
			ln_prior = 0.0;
			}

		double GetLnEdgeLenPrior(Tree &t)
			{
			ln_prior = 0.0;
			t.PostorderTraverse(boost::bind(&ObsoleteEdgeLenPriorCalculator::NodeWork, this,_1));
			//t.PostorderTraverse(boost::function1<bool, TreeNode*>(std::bind1st(std::mem_fun(&ObsoleteEdgeLenPriorCalculator::NodeWork), this)));
			return ln_prior;
			}

		bool NodeWork(TreeNode * nd)
			{
			ln_prior += probDist->GetRelativeLnPDF(nd->GetEdgeLen().f);
			return true;
			}
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Wish to find p such that (1-p)/(1 - p^(n-2)) = q, where q is the 
|
*/
class CondGeomDistrMatcher : public PhoFunction
	{
	public:
		double Evaluate(double *X);
	};

double CondGeomDistrMatcher::Evaluate(double *X)
	{
	double pi = *X;
	double p = (15.0/26.0);
	double pp = (1.0 - pi) / (1 - pow(pi, 3.0));
	double retval = (p - pp);
	if (retval < 0.0)
		return -retval;
	else
		return retval;
	}

class ObsoleteKappaDistribution : public ProbabilityDistribution
	{
	public:
					ObsoleteKappaDistribution(HKYAdHocEvaluator &h) : hky(h) {}

		// These are pure virtuals in the ProbabilityDistribution base class
		// that are not used by SliceSampler - defined only to satisfy the compiler
		//
		bool		IsDiscrete() const						{return false;}
		string 	GetDistributionName() const				{return string();}
		string 	GetDistributionDescription() const		{return string();}
		double 		GetMean() const							{return 0.0;}
		double 		GetVar() const							{return 0.0;}
		double 		GetStdDev() const						{return 0.0;}
		double		GetCDF(double ) const					{return 0.0;}
		double		Sample(Lot &) const					{return 0.0;}
		void 		SetMeanAndVariance(double , double )	{}

		// This is the only pure virtual function that is actually used by SliceSampler
		//
		//double		GetRelativeLnPDF(double x) const;

		HKYAdHocEvaluator &hky;
	};
#endif
CmdResult PhoFloor::HandleCopying()
	{
	string tmp, newick;
	NxsOutputStream * outStream = outputMgrRef.GetOutputStreamPtr();
	if (outStream != NULL)
		*outStream << "Congratulations, you made it to the copying command" << endl;
	return kCmdSucceeded;
	}

bool PhoFloor::ParseHelp(NxsToken &token)
	{
	NxsCommandManager::ProcessHelp(token);
	return true;
	}

CmdResult PhoFloor::HandleAlias(AliasSettings *s)
	{
	NxsOutputStream *outStream = outputMgrRef.GetOutputStreamPtr();
	typedef map<string, string, NxsStringNoCaseLess> StrMap;
	StrMap &aliasMap = *NxsCommandManager::GetAliasMap();
	if (s->aliasKey.empty())
		{
		if (outStream != NULL)
			{
			if (aliasMap.empty())
				*outStream << "There are no aliases.";
			else
				{
#				if defined(NCL_HAS_TABULAR_OUTPUT)
					NxsTable &outTable = *outStream->GetTablePtr();
					outTable.AddString("Alias");
					outTable.SetLeftMargin();
					outTable.AddString("Expansion");
					outTable.SetRightMargin();
					outTable.HyphenLine();
					outTable.SetTopMargin();
					for (StrMap::const_iterator mIt = aliasMap.begin(); mIt != aliasMap.end(); ++mIt)
						{
						outTable.AddString(mIt->first);
						outTable.AddString(mIt->second);
						}
					outTable.SetBottomMargin();
					outStream->PrintTable(&outTable);
#				endif
				}
			*outStream << endl;
			}
		return kCmdSucceeded;
		}
	if (s->expansion.empty())
		{
		StrMap::iterator amIt = aliasMap.find(s->aliasKey);
		if (amIt != aliasMap.end())
			{
			if (outStream != NULL)
				*outStream << "The alias " << amIt->first << " is now undefined." << endl;
			aliasMap.erase(amIt);
			}
		else
			{
			if (outputMgrRef.GetWarningStreamPtr() != NULL)
				*outputMgrRef.GetWarningStreamPtr() << "The alias " << amIt->first << " does not exist." << endl;
			}
		return kCmdSucceeded;
		}
	NxsToken expansionToken(s->expansion);
	++expansionToken;
	const string firstExpansionToken = expansionToken.GetToken();
	if (EqualsCaseInsensitive(firstExpansionToken, s->aliasKey))
		{
		if (outputMgrRef.GetErrorStreamPtr() != NULL)
			*outputMgrRef.GetErrorStreamPtr() << "The alias (" << s->aliasKey << "->" << s->expansion << ") is illegal because it expands to itself." << endl;
		return kCmdFailedSilent;
		}
	StrMap::const_iterator camIt = aliasMap.find(firstExpansionToken);
	if (camIt != aliasMap.end())
		{
		//	The new alias expands to a match of another alias.
		//	We need to verify that the new key/expansion is legal.
		//	ie there are no circular alias-expands to alias-expands to alias.. structures
		//	(these would lead to infinite loops in command processing.
		//
		while(camIt != aliasMap.end())
			{
			NxsToken nextPutativeAliasToken(camIt->second);
			++nextPutativeAliasToken;
			string nextPutativeAlias = nextPutativeAliasToken.GetToken();
			if (EqualsCaseInsensitive(nextPutativeAlias, s->aliasKey))
				{
				if (outputMgrRef.GetErrorStreamPtr() != NULL)
					*outputMgrRef.GetErrorStreamPtr()  << "The alias (" << s->aliasKey << "->" << s->expansion << ") conflicts with previously defined aliases" << endl;
				return kCmdFailedSilent;
				}
			camIt = aliasMap.find(nextPutativeAlias);
			}
		}
	aliasMap[s->aliasKey] = s->expansion;
	return kCmdSucceeded;
	}

#endif //defined(READ_PHYCAS_BLOCK)
