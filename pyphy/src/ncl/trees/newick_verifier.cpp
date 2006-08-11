//#include "phycas/force_include.h"
#include "pyphy/src/ncl/trees/newick_verifier.hpp"
#include "pyphy/src/ncl/nxs_token.hpp"
#include "pyphy/src/ncl/taxa/base_taxa_manager.hpp"
#include "pyphy/src/ncl/nxs_exception.hpp"
using std::string;
// undef ALLOW_EDGE_LENGTHS_ON_ROOT to throw exceptions when there is :# after the closing parentheses
#define ALLOW_EDGE_LENGTHS_ON_ROOT 
using std::set;
/*----------------------------------------------------------------------------------------------------------------------
|	returns the index of the tree with the name "label", or -1 if that name is unknown.  (not case sensitive)
*/
unsigned NewickVerifier::FindIndexForTaxon(string s, bool allowNewNames) const
	{
	Capitalize(s);
	TaxNameTransTable::const_iterator trIt = translationTable.find(s);
	if (trIt == translationTable.end())
		return allowNewNames ? baseTaxaMgr.FindNewOrExistingIndex(s) : baseTaxaMgr.FindIndex(s);
	return trIt->second.second;
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Very rudimentary tree node class used to check for errors Newick tree descriptions are being read. 
*/
class SimpleNode
	{
		BaseTaxaManager & baseTaxaMgr;
		SimpleNode *lChild;		/* pointer to left child node */
		SimpleNode *rSib;		/* pointer to right sibling node */
		SimpleNode *par;		/* pointer to parent node */
		double 		edgelen;	/* branchlength (whatever follows the : interpreted as a double ) */
		
		void 		ReadTaxSetIntoTree(const NxsIndexSet *taxSet, const NxsToken &taxSetName, set<unsigned> &taxIndicesFound);
	public :
		STATIC_DATA FullTreeDescription::EdgeLengthMode edgeLenReading;	/* indicates whether there are branch lengths */
		STATIC_DATA const NewickVerifier  *treesBlock;			/*pointer to the current trees block reader*/
		STATIC_DATA string 	translatedTreeDesc;	/*Newick representation with taxon indices according to the taxa manager (not the translate block)*/
		STATIC_DATA bool 		polytomyFound;		/*True if there is at least one polytomy in the tree*/
		STATIC_DATA bool 		nodeDegreeTwoFound; /*True if there is at least one node with only one descendant in the tree*/
		STATIC_DATA bool 		internalNodeLabelsFound; /*True if there is at least one of the internal nodes was given a name */
		
		string 	name;	/*the node's name */
		SimpleNode(BaseTaxaManager &);
		void 	ReadTree(NxsToken &, set<unsigned> &taxIndicesFound, VecString *newTaxaNamesIfFound);
		~SimpleNode();
		friend class NewickVerifier;
	};

/*----------------------------------------------------------------------------------------------------------------------
|	recursively deletes the subtree of simple nodes.
*/
SimpleNode::~SimpleNode()	
	{
	delete lChild; 
	delete rSib;
	}


	
FullTreeDescription::EdgeLengthMode SimpleNode::edgeLenReading = FullTreeDescription::kUnknownEdgeLengthMode;	
const NewickVerifier *SimpleNode::treesBlock = NULL;	
string 	SimpleNode::translatedTreeDesc;
bool 		SimpleNode::polytomyFound = false;
bool 		SimpleNode::nodeDegreeTwoFound = false;
bool 		SimpleNode::internalNodeLabelsFound = false;


/*----------------------------------------------------------------------------------------------------------------------
|	Creates an unattached node (all pointers NULL) with edgelen of DBL_MAX
*/
SimpleNode::SimpleNode(BaseTaxaManager & taxMgr)
	:baseTaxaMgr(taxMgr),
	lChild(NULL),
  	rSib(NULL),
  	par(NULL),
  	edgelen(DBL_MAX)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called when a taxset is found in the tree description.  The tax set taxGroup contains taxA, taxB, and taxC, then
|	a tree description with (taxGroup), should be treated the same as if it said ( taxA, taxB,  taxC ).  The branches 
|	cannot be assigned a length.
|	This function performs the expansion in the STATIC_DATA string translatedTreeDesc (so this variable never has taxset 
|	names) and it check to make sure all taxa in tax haven't been added to the tree already.
*/
void SimpleNode::ReadTaxSetIntoTree(
  const NxsIndexSet *taxSet, /* indices of the taxa that make up the taxset that was found in the tree description */
  const NxsToken &token, /*token that is supplying the tree description (used only if an NxsException exception need to be generated */
  set<unsigned> &taxIndicesFound)	/*set that stores the indices of all of the taxa included in the tree thus far */
	{
	NxsIndexSet::const_iterator numIt = taxSet->begin();
	std::pair< set<unsigned>::iterator, bool> tifRet = taxIndicesFound.insert(*numIt);
   	bool duplicate = false;
	if (!tifRet.second)
		duplicate = true;
	else
		{
		translatedTreeDesc << *numIt;
		name << *numIt;
		}
	SimpleNode *tempNode = this;
	if (numIt != taxSet->end())
		{
		if (edgeLenReading == FullTreeDescription::kEdgeLengthsFound)
			edgeLenReading = FullTreeDescription::kSomeEdgeLengths;
		else if (edgeLenReading == FullTreeDescription::kUnknownEdgeLengthMode)
			edgeLenReading = FullTreeDescription::kNoEdgeLengthsFound;
	  	for (numIt++; !duplicate && numIt != taxSet->end(); ++numIt)
			{
			translatedTreeDesc << ',' << *numIt;
			tempNode->rSib = new SimpleNode(baseTaxaMgr);
			tempNode = tempNode->rSib;
			tempNode->par = par;
			tempNode->name << *numIt;
			tifRet = taxIndicesFound.insert(*numIt);
			if (!tifRet.second)
	  			{
	  			duplicate = true;
	  			break;
	  			}
	  		}
  		}
	if (duplicate)
		{
		string eMsg;
		eMsg << "The taxset " << token.GetTokenReference() << " includes taxon number " << *numIt + 1 << " which was already included in the tree description";
		throw NxsException(eMsg, token);
		}	
	}
	
/*----------------------------------------------------------------------------------------------------------------------
|	Recursive function used to read in a sub-clade.  Checks for legal tree decription and taxa names.  Builds up the
|	tree description (using the taxa indices from the trees block's taxa manager) in the STATIC_DATA translatedTreeDesc.
|	Stores information such as whether or not there are polytomies in SimpleTree STATIC_DATA fields.
*/
void SimpleNode::ReadTree(
  NxsToken &token, 
  set<unsigned> &taxIndicesFound,
  VecString *newTaxaNames)
	{
#	if defined(C_FUNCS_IN_STD_NAMESPACE)
		using std::isdigit;
#	endif
	assert(treesBlock != NULL);
	delete rSib;
	rSib = NULL;
	delete lChild;
	lChild = NULL;
  	name.clear();
  	string eMsg;
  	if (token.GetTokenReference() == '(')
  		{
  		translatedTreeDesc << '(';
  		++token;
  		lChild = new SimpleNode(baseTaxaMgr);
  		lChild->par = this;
  		lChild->ReadTree(token, taxIndicesFound, newTaxaNames);
  		}
  	
  	if (!token.GetTokenReference().empty())
  		{
  		const char &ch(token.GetTokenReference()[0]);
  		if (ch != ',' && ch != ':' && ch != ';' && ch != ')')
	  		{
	  		//	This token should be name (label, index, or translation table key)
	  		//
	  		unsigned n = treesBlock->FindIndexForTaxon(token.GetTokenReference(), false);
			if (n == UINT_MAX)
	  			{
	  			//	see if the token is a tax set name
	  			//
	  			const NxsIndexSet *taxSet = baseTaxaMgr.GetSet(token.GetTokenReference());
	  			if (taxSet != NULL)
	  				 {
	  				 if (lChild != NULL)
		  				{
		  				eMsg << "Illegal use of a TaxSet (" << token.GetTokenReference() << ") as the label for an internal node in the tree.";
		  				throw NxsException(eMsg, token);
		  				}
		  			string taxSetName(token.GetTokenReference());
		  			ReadTaxSetIntoTree(taxSet, token, taxIndicesFound);
		  			token.ReadSingleCharacter();
		  			if (token.GetTokenReference()[0] == ':')
		  				{
		  				eMsg << "Illegal use of a branch length specifier following aTaxSet (" << taxSetName << ").";
		  				throw NxsException(eMsg, token);
		  				}
		  			}
	  			else 
	  				{
	  				if (newTaxaNames != NULL)
		  				{
						n = treesBlock->FindIndexForTaxon(token.GetTokenReference(), true);
						if (n != UINT_MAX)
							{
							const UInt storedNNames = baseTaxaMgr.GetSize();
							const UInt nNewTaxaAdded = (unsigned)newTaxaNames->size();
							for (UInt i = storedNNames + nNewTaxaAdded + 1; i <= (n+1); ++i)
								{
								std::string stringForm;
								stringForm << i;
								newTaxaNames->push_back(stringForm);
								}
							}
						else
							{
							n = baseTaxaMgr.GetSize();
							string newName(token.GetTokenReference());
							NStrCaseInsensitiveEquals elementCompare(newName);
							VecString::iterator nnIt = find_if(newTaxaNames->begin(), newTaxaNames->end(), elementCompare);
							if (nnIt == newTaxaNames->end())
								{
		#						if defined(TREAT_TAXA_NAMES_AS_TAXA_NUMBERS_ONLY)
									pair<bool, unsigned> nPair = IsAnUnsigned(newName);
									if (!nPair.first || nPair.second == 0)
										{
										eMsg << "Expecting a (positive) taxon number instead of \"" << newName << "\" because this version was compiled in TREAT_TAXA_NAMES_AS_TAXA_NUMBERS_ONLY mode.";
										throw NxsException(eMsg, token);
										}
									unsigned currSize = newTaxaNames->size();
									string tmp;
									while (currSize != nPair.second)
										{
										++currSize;
										tmp.clear();
										tmp << currSize;
										newTaxaNames->push_back(tmp);
										}
									n += nPair.second -1;
		#						else
									newTaxaNames->push_back(newName); //@@ should look for the name in the vector to
									n += (unsigned)newTaxaNames->size() - 1;
		#						endif
								}
							else
								n += (unsigned)distance(newTaxaNames->begin(), nnIt);
							}
		  				}
		  			if (n == UINT_MAX)
						{
		  				eMsg << "An unknown taxon name (" << token.GetTokenReference() << ") was found in a tree description.";
		  				throw NxsException(eMsg, token);
		  				}
	  				}
	  			}
	  		if (n != UINT_MAX) 	// if a tax set was already read n == UINT_MAX  (so we can skip this loop)
	  			{
		  		std::pair< set<unsigned>::iterator, bool> tifRet = taxIndicesFound.insert(n);
		  		if (!tifRet.second)
		  			{
		  			eMsg << "The taxon " << token.GetTokenReference() << " appears more than once in the tree description";
		  			throw NxsException(eMsg, token);
		  			}
		  		translatedTreeDesc << n;
		  		name << n;
		  		if (lChild != NULL)
		  			internalNodeLabelsFound = true;
		  		token.ReadSingleCharacter();
	  			}
	  		}
  		}
  	else
  		{
  		if (par == NULL)
  			return;
  		throw NxsException("Missing closing parentheses in tree.", token);
  		}
  	//	name or subtree description must be present.  () is not allowed
  	//
 	if (lChild == NULL && name.empty())
		{
		eMsg << "Expecting either a taxon name of subtree description.  Found " << token.GetTokenReference()[0];
		throw NxsException(eMsg, token);
		}
  	if (token.GetTokenReference()[0] == ':')
		{
		//	Read in edgelength.  Throw an exception if it can't be converted to a double.
  		//
#		if !defined(ALLOW_EDGE_LENGTHS_ON_ROOT)
	 		if (par == NULL)
				throw NxsException("Unexpected colon.  Branch length specifiers are not allowed for the root of the tree (they have no meaning).", token);
#		endif
		if (edgeLenReading == FullTreeDescription::kNoEdgeLengthsFound)
			edgeLenReading = FullTreeDescription::kSomeEdgeLengths;
		else if (edgeLenReading == FullTreeDescription::kUnknownEdgeLengthMode)
			edgeLenReading = FullTreeDescription::kEdgeLengthsFound;
		if (!token.ReadDoubleToken(&edgelen))
			{
  			eMsg << "Expecting a positive number for the branch lengths in the tree, but found " << token.GetTokenReference() << ".\n";
  			throw NxsException(eMsg, token);
  			}
  		
		translatedTreeDesc << ':' << token.GetTokenReference();
  		token.ReadSingleCharacter();
  		}
  	else 
  		{
  		if (edgeLenReading == FullTreeDescription::kEdgeLengthsFound && token.GetTokenReference()[0] != ';')
			edgeLenReading = FullTreeDescription::kSomeEdgeLengths;
		else if (edgeLenReading == FullTreeDescription::kUnknownEdgeLengthMode)
			edgeLenReading = FullTreeDescription::kNoEdgeLengthsFound;
  		}
  		
 	if (token.AtEOF() || token.GetTokenReference()[0] == ';')
  		{
  		if (par == NULL)
  			return;
  		throw NxsException("Missing closing parentheses in tree.", token);
  		}
  	else if (token.GetTokenReference()[0] == ',')
  		{
  		//	read the sister subtree description and return
  		//
  		if (par == NULL)
  			throw NxsException("Found a comma outside any parentheses in a tree", token);
  		SimpleNode *tempNode = this;
  		while (tempNode->rSib != NULL)
  			tempNode = tempNode->rSib;
  		tempNode->rSib = new SimpleNode(baseTaxaMgr);
  		tempNode->rSib->par = par;
  		translatedTreeDesc << ',';
  		if (par->lChild != this)
  			polytomyFound = true;
  		++token;
  		tempNode->rSib->ReadTree(token, taxIndicesFound, newTaxaNames);
  		return;
  		}
 	else if (token.GetTokenReference()[0] == ')')
		{
		//	Finished reading this node, return.
  		//
  		if (par == NULL)
			throw NxsException("Unbalanced tree.  Too many closing parentheses were found.", token);
		if (par->lChild->rSib == NULL)
			{
			nodeDegreeTwoFound = true;
			throw NxsException("A node with only one descendant was found.", token);	
			}
		translatedTreeDesc << ')';
  		++token;
		return;
  		}

  	//	found some token other than , ) or ; after the node name or branch length
  	//
  		
  	eMsg << "Unrecognized " << token.GetTokenReference() << " in tree (perhaps a ";
  	if (token.GetTokenReference()[0] == '.' || isdigit(token.GetTokenReference()[0]))
  		eMsg << "colon, comma, or";
  	else 
  		eMsg << "comma or";
  	eMsg << " parentheses is missing).";	
  	throw NxsException(eMsg, token);	
  	}
  


void NewickVerifier::ReadTreeAndAlertTaxaManager(NxsToken &token, FullTreeDescription *ftd)
  	{
  	if (allowImplicitNames)
		{
		VecString newNames;
		*ftd = ReadNewickTree(token, &newNames);
		if (!newNames.empty())
			baseTaxaMgr.AppendTaxa(newNames);
		}
	else
		*ftd = ReadNewickTree(token);
  	}
  	
FullTreeDescription NewickVerifier::ReadNewickTree(NxsToken &token, VecString *newTaxaNames, const string treeName) const
 	{
  	string errorMsg;
  	//	Next token(s) could be a command comments.  Pay attention to [&U] or [&R] or the beginning of the tree
	//
	bool rootedSpecified = false;
	bool unrootedSpecified = false;
	bool cmdCommentFound = token.ReadCommandCommentOrToken();
	while (cmdCommentFound)
		{
		unsigned commandCommentIndex = token.GetNextCommentIndex('&', 0);
		if (commandCommentIndex != UINT_MAX)
			{
			NxsComment cmdComment = token.GetComment(commandCommentIndex);
			if (cmdComment.GetLocationInToken() < 1)
				{
				if (StrEquals(cmdComment.GetCommentText(), "U", ncl::kStringNoRespectCase))
					unrootedSpecified = true;
				else if (StrEquals(cmdComment.GetCommentText(), "R", ncl::kStringNoRespectCase))
					rootedSpecified = true;
				if (rootedSpecified && unrootedSpecified)
					throw NxsException("Tree was specified as being both rooted [&R] and unrooted [&U]" , token);
				}
			}
		cmdCommentFound = token.ReadCommandCommentOrToken();
		}
	//	Prepare all of the SimpleNode statics for the checking of a new tree
	//
	SimpleNode rootN(baseTaxaMgr);
	SimpleNode::treesBlock = this;
 	SimpleNode::edgeLenReading = FullTreeDescription::kUnknownEdgeLengthMode;
	SimpleNode::translatedTreeDesc.clear();
	SimpleNode::polytomyFound = false;
	SimpleNode::nodeDegreeTwoFound = false;
	SimpleNode::internalNodeLabelsFound = false;
	set<unsigned> leafSet;
	
	// Use the SimpleNode class to verify that the tree description is legal
	//
	if (token.GetTokenReference() != '(')
		{
		errorMsg << "Expecting the tree description to start with \"(\" instead of " << token.GetToken();
		throw NxsException(errorMsg, token);
		}
	rootN.ReadTree(token, leafSet, newTaxaNames);
		
	// make sure there are enough nodes for a minimal tree description
	//
	if (rootN.lChild == NULL)
		{
		errorMsg << "A tree must contain three or more nodes, only one was found";
		throw NxsException(errorMsg, token);
		}
	assert(rootN.lChild != NULL && rootN.lChild->rSib != NULL);
	bool rooted = rootedSpecified;
	if (unrootedSpecified)
		{
		if (rootN.lChild->rSib->rSib == NULL)
			{
			errorMsg << "An unrooted tree was specified, but the lowest node of the tree has only two descendants";
			throw NxsException(errorMsg, token);
			}
		}
	else if(!rootedSpecified)
		{
		if (rootN.lChild->rSib->rSib == NULL)
			rooted = true;
		}
	// 	Store the newick description (with taxon indices in terms of the TaxaListManager) and information 
	//	about the tree (e.g. whether there were polytomies etc.
	//
	return FullTreeDescription(treeName, SimpleNode::translatedTreeDesc, NxsIndexSet(leafSet), rooted, SimpleNode::edgeLenReading, SimpleNode::polytomyFound, SimpleNode::internalNodeLabelsFound);
  	}


