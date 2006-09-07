#ifndef NCL_NXSTREE_DESCRIPTION_H
#define NCL_NXSTREE_DESCRIPTION_H

#include "phypy/src/ncl/misc/nxs_index_set.hpp"

/*----------------------------------------------------------------------------------------------------------------------
|	Class that stores a Newick representation of the tree (with indices for the taxa that correspond to those in the 
|	taxa manager that was used when the block was read) as well as other information about the tree that might be useful
|	in deciding whether the tree can be read by a program (its name, whether or not it is rooted, whether or not there 
|	are edge lengths, whether there are polytomies, and whether or not some of the internal nodes are named).
*/
class FullTreeDescription
	{
	public:
		enum EdgeLengthMode				/* specifies status of edge lengths */
			{
			kUnknownEdgeLengthMode, 	/* do not know whether edge lengths were specified */
			kEdgeLengthsFound, 			/* edge lengths were found when reading tree description */
			kNoEdgeLengthsFound, 		/* no edge lengths were found when reading tree description */
			kSomeEdgeLengths			/* some, but not all, edge lengths were found when reading tree description */
			};
		
							FullTreeDescription();
							FullTreeDescription(
								const std::string &	n, 
								const std::string &	descrip, 
								const NxsIndexSet &	leafSet, 
								bool				isRooted = false, 
								EdgeLengthMode		hasBrLens = kUnknownEdgeLengthMode, 
								bool				polytomous = false, 
								bool				internalNames = false
								);

		bool				HasEdgeLengths() const;
		const std::string & GetName() const {return name;}
		bool				HasEdgeLengths();
		std::string			BriefReport() const;
		void				RemoveEdgeLengths();

		std::string   		name;
		std::string			newick;
		bool				rooted;
		EdgeLengthMode		hasEdgeLens;
		bool				hasPolytomies;
		NxsIndexSet			leafSet;
		bool				hasInternalTaxonNames;
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Initializes `hasEdgeLens' to `kUnknownEdgeLengthMode' and `nTips' to 0.
*/
inline FullTreeDescription::FullTreeDescription()
  : hasEdgeLens(kUnknownEdgeLengthMode)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if edge lengths were found in the tree description, false otherwise.
*/
inline bool FullTreeDescription::HasEdgeLengths() const
	{
	return (hasEdgeLens == kEdgeLengthsFound);
	}

typedef std::vector<FullTreeDescription> 	VecOfTreeDescriptions;

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
inline FullTreeDescription::FullTreeDescription(
  const std::string & n, 			/*	the name of the tree */
  const std::string & description, 	/*	Newick description (with indices to indicate taxa) */
  const NxsIndexSet & leaves,		/*	number of tips in the tree */
  bool 				isRooted, 		/*	true if rooted */
  EdgeLengthMode	hasBrLens,		/*	value from enum indicating whether none, some, or all branches have lengths assigned */
  bool 				polytomous, 	/*	true if there are any polytomies */
  bool 				internalNames) 	/*	true if any internal nodes are named */
	:name(n),
	newick(description),
	rooted(isRooted),
	hasEdgeLens(hasBrLens),
	hasPolytomies(polytomous), 
	leafSet(leaves), 
	hasInternalTaxonNames(internalNames)
	{
	}

#endif
