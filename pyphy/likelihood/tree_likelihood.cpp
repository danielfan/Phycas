#include "phycas/force_include.h"
#include "CipresCommlib/CipresDataMatrixHelper.h"
#include "pyphy/phylogeny/basic_tree.hpp"
#include "pyphy/likelihood/likelihood_models.hpp"
#include "pyphy/likelihood/tree_likelihood.hpp"
#include "pyphy/likelihood/tip_data.hpp"
#include "pyphy/likelihood/internal_data.hpp"
#include "pyphy/likelihood/sim_data.hpp"
#include "pyphy/common/pyphy_string.hpp"
#include "pyphy/prob_dist/basic_lot.hpp"
#include <boost/format.hpp>
#include <numeric>

namespace phycas
{

#if 0
class CalcTransitionMatrixForOneNode : public std::unary_function<TreeNode &, void>
	{
	private:
		TreeLikelihood & treelike;

	public:
		CalcTransitionMatrixForOneNode(TreeLikelihood & tl) : treelike(tl) {}
		void operator()(TreeNode & nd)
			{
#			error do not use unless root node special case is taken into account
			if (!nd.IsRoot())
				{
				double edge_len = nd.GetEdgeLen();
				if (nd.IsTip())
					{
					TipData & ndTD = *(nd.GetTipData());
					treelike.calcTMatForSim(ndTD, edge_len);
					}
				else
					{
					InternalData & ndID	= *(nd.GetInternalData());
					treelike.calcPMat(ndID, edge_len);
					}
				}
			}
	};
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Simulates data for `nchar' characters using the current model and edge lengths. Transition matrices are recomputed
|	if `refresh_probs' is true; if `refresh_probs' is false, assumes transition probability matrices are up-to-date. 
|	Only the `num_states' primary states are generated, and the state created for each node is stored initially in the 
|	`state' data member of the TipData or InternalData structure associated with the tip or internal node, respectively.
|	This function expects these data structures to be in place (the TreeLikelihood::prepareForSimulation and 
|	TreeLikelihood::prepareForLikelihood functions both perform this task). As states are generated for tip nodes, they
|	are copied into a temporary pattern inside `sim_data', and this pattern is then inserted into a pattern map inside 
|	`sim_data' when it is completed. Assumes that `nchar' is greater than zero and that the shared pointers `sim_data' 
|	and `rng' actually point to real objects.
|	
|	The following makes use of T matrices, which are transposed, augmented transition matrices. These are transposed 
|	because the "from" states form the columns rather than the rows. They are augmented because, ordinarily, there are
|	additional rows corresponding to ambiguities observed in some tip nodes. With simulated data, however, there are 
|	never any ambiguities, so in this case T matrices are nothing more than transposed transition matrices. 
|	Important: note that T matrices are only used in the TipData structures; InternalData structures store normal 
|	untransposed transition probability matrices in which the rows form the "from" states. 
|	
|	Here is an example of a T matrix created using the JC model and an edge length equal to 0.1:
|>	
|	            |--------------- from state -------------|
|	                0          1          2          3
|	  	            A          C          G          T
|	t  0   A     0.90638    0.03121    0.03121    0.03121
|	o  1   C     0.03121    0.90638    0.03121    0.03121
|	   2   G     0.03121    0.03121    0.90638    0.03121
|	s  3   T     0.03121    0.03121    0.03121    0.90638
|	t  4   N     1.00000    1.00000    1.00000    1.00000 \
|	a  5 {GT}    0.06241    0.06241    0.93759    0.93759  | These rows not present if prepareForSimulation function
|	t  6 {ACT}   0.96879    0.96879    0.09362    0.96879  | was used to create the TipData structures
|	e  7 {AG}    0.93757    0.06241    0.93759    0.06241 /
|>
|	The `pMatrixTranspose' data member in TipData structures holds the array of T matrices (one T matrix for each 
|	rate category).
*/
void TreeLikelihood::simulateImpl(SimDataShPtr sim_data, TreeShPtr t, LotShPtr rng, unsigned nchar, bool refresh_probs)
	{
	assert(sim_data);
	assert(rng);
	assert(nchar > 0);

	// Recalculate transition probabilities if requested
	//@POL using for_each would simplify this
	if (refresh_probs)
		{
		preorder_iterator nd = t->begin();

		// First preorder node is the root node and represents a special case
		// Its transition matrices must be computed using the "subroot" node's edge length
		// The subroot node's transition matrices need not be calculated 
		TipData & ndTD = *(nd->GetTipData());
		TreeNode * subroot = nd->GetLeftChild();
		assert(subroot);
		assert(!subroot->GetRightSib());	//@POL need to create a IsSubroot() member function for TreeNode
		calcTMatForSim(ndTD, subroot->GetEdgeLen());
		++nd;

		// Skip subroot node as its transition matrices are never used and thus do not need to be computed
		++nd;

		// Process the remaining nodes in the tree
		for (; nd != t->end(); ++nd)
			{
			if (nd->IsTip())
				{
				TipData & ndTD = *(nd->GetTipData());
				calcTMatForSim(ndTD, nd->GetEdgeLen());
				}
			else
				{
				InternalData & ndID	= *(nd->GetInternalData());
				calcPMat(ndID, nd->GetEdgeLen());
				}
			}
		}

	// Create a vector of cumulative state frequencies to use in choosing starting states
	const std::vector<double> & freqs = model->getStateFreqs();
	std::vector<double> cum_freqs(num_states, 0.0);
	std::partial_sum(freqs.begin(), freqs.end(), cum_freqs.begin());

	// Create a vector of cumulative rate probabilities to use in choosing relative rates
	std::vector<double> cum_rate_probs(num_rates, 0.0);
	std::partial_sum(rate_probs.begin(), rate_probs.end(), cum_rate_probs.begin());

	sim_data->resetPatternLength(t->GetNTips());
	sim_data->wipePattern();

	for (unsigned character = 0; character < nchar; ++character)
		{
		// Choose a rate for this character (actually, choose index, the actual rate is rate_means[r])
		unsigned r = 0;
		if (num_rates > 1)
			{
			// warning: removing the if statement will invalidate all examples involving simulated data with rate
			// homogeneity because of the call to rng->Uniform here!
			r = (unsigned)(std::lower_bound(cum_rate_probs.begin(), cum_rate_probs.end(), rng->Uniform()) - cum_rate_probs.begin());
			}

		// Generate the starting state
		int8_t j = (unsigned)(std::lower_bound(cum_freqs.begin(), cum_freqs.end(), rng->Uniform()) - cum_freqs.begin());

		// Assign starting state to the tip node currently serving as the root of the tree
		preorder_iterator nd = t->begin();
		TipData & rootTD = *(nd->GetTipData());
		rootTD.state = j;

		sim_data->setState(nd->GetNodeNumber(), j);

		// Go ahead and generate the state for the (only) descendant of the root node (the "subroot" node)
		// Note that the root node's T matrix is used for this calculation; the P matrix of the subroot node
		// is never computed
		unsigned parent_state = (unsigned)j;

		// Get the T matrix for the tip node serving as the root
		double * * Tmatrix = rootTD.pMatrixTranspose[r];

		// Choose a uniform random deviate
		double u = rng->Uniform();

		// Spin the roulette wheel to choose a state for the subroot node
		double cum = 0.0;
		unsigned i = 0;
		for (; i < num_states; ++i)
			{
			double pr = Tmatrix[i][parent_state];
			//std::cerr << str(boost::format("Tmatrix[%d][%d] = %f") % i % parent_state % pr) << std::endl;
			cum += pr;
			if (u < cum)
				break;
			}

		// Increment iterator so that nd now refers to the subroot (sole descendant of the root)
		++nd;

		// Assign the new state to the subroot node
		InternalData & ndID = *(nd->GetInternalData());
		ndID.state = (int8_t)i;
		//std::cerr << "  Assigning state " << i << " to node " << nd->GetNodeNumber() << std::endl;

		// Walk the remainder of the tree using the preorder sequence, generating data for each node along the way
		for (++nd; nd != t->end(); ++nd)
			{
			// Get state of parent of nd
			TreeNode * parent = nd->GetParent();
			parent_state = UINT_MAX;
			if (parent->IsTip())
				{
				TipData * parentTD = parent->GetTipData();
				parent_state = (unsigned)parentTD->state;
				}
			else
				{
				InternalData * parentID = parent->GetInternalData();
				parent_state = (unsigned)parentID->state;
				}
			assert(parent_state < num_states);

			if (nd->IsTip())
				{
				// Get the T matrix
				TipData & ndTD = *(nd->GetTipData());
				double * * Tmatrix = ndTD.pMatrixTranspose[r];

				// Choose a uniform random deviate
				double u = rng->Uniform();

				// Spin the roulette wheel and assign a state to nd
				double cum = 0.0;
				unsigned i = 0;
				for (; i < num_states; ++i)
					{
					double pr = Tmatrix[i][parent_state];
					//std::cerr << str(boost::format("Tmatrix[%d][%d] = %f") % i % parent_state % pr) << std::endl;
					cum += pr;
					if (u < cum)
						break;
					}
				ndTD.state = (int8_t)i;
				sim_data->setState(nd->GetNodeNumber(), (int8_t)i);

				//std::cerr << "  Assigning state " << i;
				}
			else
				{
				// Get the T matrix
				InternalData & ndID = *(nd->GetInternalData());
				double * * Pmatrix = ndID.pMatrices[r];

				// Choose a uniform random deviate
				double u = rng->Uniform();

				// Spin the roulette wheel and assign a state to nd
				double cum = 0.0;
				unsigned i = 0;
				for (; i < num_states; ++i)
					{
					double pr = Pmatrix[parent_state][i];
					//std::cerr << str(boost::format("Pmatrix[%d][%d] = %f") % parent_state % i % pr) << std::endl;
					cum += pr;
					if (u < cum)
						break;
					}
				ndID.state = (int8_t)i;
				//std::cerr << "  Assigning state " << i;
				}

			//std::cerr << " to node " << nd->GetNodeNumber() << std::endl;
			}
		//std::cerr << std::endl;

		// We are now finished simulating data for one character, so insert the pattern just generated
		// into the pattern map maintained by sim_data; the 1 means that the count for this pattern
		// should be incremented by 1
		sim_data->insertPattern(1);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes and returns the log-likelihood score of the supplied tree using the supplied model. Assumes that the
|	function prepareForLikelihood() has already been called for the same tree. If setNoData function was called more
|	recently than the setHaveData function, returns 0.0 immediately.
*/
double TreeLikelihood::calcLnL(
  TreeShPtr t)		/**< is the tree for which the log-likelihood is to be computed */
	{
	if (no_data)
		{
		++nevals;
		return 0.0;
		}

	TreeNode * rootNd = t->GetFirstPreorder();
	TreeNode * firstNd = rootNd->GetLeftChild();
	assert(firstNd);
	for (TreeNode * nd = t->GetLastPreorder(); nd != 0; nd = nd->GetNextPostorder())
		{
		if (nd->IsInternal())
			{
#			if 0	// if reinstated, also reinstate code in TreeNode::SetEdgeLen
				// If this node needs to be recalculated, so does its parent
				TreeNode * parent = nd->GetParent();
				assert(parent);
				parent->SelectNode();

				// Can go ahead and unselect this internal node now. The only exception is the node
				// that is the only child of the root (tip) node. This node must remain selected until
				// the information from the root is incorporated
				if (nd != firstNd)
					nd->UnselectNode();
#			endif

			InternalData & nd_ID = *(nd->GetInternalData());
			TreeNode * left_child = nd->GetLeftChild();
			TreeNode * second_child	= left_child->GetRightSib();

			// There are four possibilities:
			// 1. both left child and its sibling are tips
			// 2. left child is a tip, but its sibling is an internal node
			// 3. left child is an internal node, but its sibling is a tip
			// 4. both left child and its sibling are internal nodes

			if (left_child->IsTip())
				{
				// left child is a tip
				TipData & leftTD = *(left_child->GetTipData());
				calcPMatTranspose(leftTD, left_child->GetEdgeLen());
				if (second_child->IsTip())
					{
					// 1. both left child and its sibling are tips
					TipData & secondTD = *(second_child->GetTipData());
					calcPMatTranspose(secondTD, second_child->GetEdgeLen());
					calcCLATwoTips(nd_ID, leftTD, secondTD);
					}
				else
					{
					// 2. left child is a tip, but its sibling is an internal node
					InternalData & secondID	= *(second_child->GetInternalData());
					calcPMat(secondID, second_child->GetEdgeLen());
					calcCLAOneTip(nd_ID, leftTD, secondID);
					}
				}
			else
				{
				// left child is an internal node
				InternalData & leftID = *(left_child->GetInternalData());
				calcPMat(leftID, left_child->GetEdgeLen());
				if (second_child->IsTip())
					{
					// 3. left child is an internal node, but its sibling is a tip
					TipData & secondTD = *(second_child->GetTipData());
					calcPMatTranspose(secondTD, second_child->GetEdgeLen());
					calcCLAOneTip(nd_ID, secondTD, leftID);
					}
				else
					{
					// 4. both left child and its sibling are internal nodes
					InternalData & secondID	= *(second_child->GetInternalData());
					calcPMat(secondID, second_child->GetEdgeLen());
					calcCLANoTips(nd_ID, leftID, secondID);
					}
				}

			// Deal with possible polytomy in which second_child has siblings
			for (TreeNode * currNd = second_child->GetRightSib(); currNd != NULL; currNd = currNd->GetRightSib())
				{
				if (currNd->IsTip())
					{
					TipData & currTD = *(currNd->GetTipData());
					calcPMatTranspose(currTD, currNd->GetEdgeLen());
					conditionOnAdditionaTip(nd_ID, currTD);
					}
				else
					{
					InternalData & currID = *(currNd->GetInternalData());
					calcPMat(currID, currNd->GetEdgeLen());
					conditionOnAdditionaInternal(nd_ID, currID);
					}
				}
			} // if (nd->IsInternal())
		} // main postorder loop

	// Deal with special case of internal node that is only child of the tip node serving as the root
	// Currently assuming the tree is unrooted and thus the root node has data like all other tips
	InternalData & firstID = *(firstNd->GetInternalData());
	//if (firstNd->IsSelected())
	//	{
		TipData & rootTD = *(rootNd->GetTipData());
		calcPMatTranspose(rootTD, firstNd->GetEdgeLen());
		conditionOnAdditionaTip(firstID, rootTD);
	//	firstNd->UnselectNode();
	//	}

	// The first internal node (child of root) now holds all the information about the likelihood of each
	// pattern for each site and each relative rate category. The tree has ripened, and now it is harvest time!
	double lnL = harvestLnL(firstID);

#	if 0
		//POL-debug
		std::string s;
		std::vector<std::string> pattern_repr;
		buildPatternReprVector(pattern_repr, t);

		std::cerr << "\ncalcLnL function";
		std::cerr << std::endl;
		std::cerr << str(boost::format("%6s  %12s  pattern") % "i" % "site-like") << std::endl;

		double *siteLike = (double *)(&site_likelihoods[0]);
		for (unsigned i = 0; i < getNPatterns(); ++i)
			{
			std::cerr << str(boost::format("%6d  %12.10f  %s") % i % siteLike[i] % pattern_repr[i]) << std::endl;
			}
#	endif

	++nevals;

	return lnL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates the TipData data structure needed to store the data for one tip (the tip corresponding to the supplied
|	`row' index in the data matrix `mat'. Returns a pointer to the newly-created TipData structure. See documentation 
|	for the TipData structure for more explanation.
*/
TipData * TreeLikelihood::allocateTipData(  //POLBM TreeLikelihood::allocateTipData
  unsigned row) 	/**< is the row of the data matrix corresponding to the data for this tip node */
	{
	std::map<int8_t, int8_t>					globalToLocal;
	std::vector<unsigned int>					stateListVec;
	std::map<int8_t, int8_t>::const_iterator	foundElement;

	int8_t *									tipSpecificStateCode	= new int8_t[(unsigned)pattern_map.size()];
	//@POL 21-Nov-2005 make tipSpecificStateCode a shared_array or a std:Vector - currently I don't think these are being deleted

	const int8_t								ns						= num_states;
	const int8_t								nsPlusOne				= num_states + 1;
	unsigned									nPartialAmbig			= 0;

	// Loop through all patterns for this row of the matrix. For each global state code encountered,
	// determine which local state code it represents, and build up the tipSpecificStateCode array
	// as we go.
	unsigned i = 0;
	for (PatternMapType::const_iterator it = pattern_map.begin(); it != pattern_map.end(); ++it, ++i)
		{
		const int8_t globalStateCode = (it->first)[row];

		if (globalStateCode < nsPlusOne)
			{
			// no partial ambiguity, but may be gap state
			tipSpecificStateCode[i] = (globalStateCode < 0 ? ns : globalStateCode);
			}
		else
			{
			// partial ambiguity
			foundElement = globalToLocal.find(globalStateCode);
			if (foundElement == globalToLocal.end())
				{
				// state code needs to be added to map
				globalToLocal[globalStateCode] = nPartialAmbig + nsPlusOne;
				stateListVec.push_back(state_list_pos[globalStateCode]);
				tipSpecificStateCode[i] = nPartialAmbig + nsPlusOne;
				nPartialAmbig++;
				}
			else
				{
				// state code is already in the map
				tipSpecificStateCode[i] = foundElement->second;
				}
			}
		}

#	if 0
		//POL-debug
		std::map<int8_t, int8_t>::const_iterator mapiter; 
		unsigned z;

		//std::cerr << "\n\nInside allocateTipData function:" << std::endl;

		//std::cerr << "\n  ns: " << (int)ns << std::endl;

		std::cerr << "\n  Adding TipData for node number " << row << ": ";
				for (z = 0; z < num_patterns; ++z)
					std::cerr << ' ' << (int)(myRow[z]);
		std::cerr << std::endl;

		//std::cerr << "\n  Contents of state_list_pos:" << std::endl;
		//for (vector<unsigned int>::const_iterator ziter = state_list_pos.begin(); ziter != state_list_pos.end(); ++ziter)
		//	std::cerr << "\n    " << (*ziter);
		//std::cerr << std::endl;

		//std::cerr << "\n  nPartialAmbig: " << nPartialAmbig << std::endl;

		//std::cerr << "\n  Contents of tipSpecificStateCode:" << std::endl;
		//for (z = 0; z < nPatterns; ++z)
		//	std::cerr << "\n    " << (int)(tipSpecificStateCode[z]);
		//std::cerr << std::endl;

		//if (globalToLocal.empty())
		//	std::cerr << "\n  globalToLocal map is empty" << std::endl;
		//else
		//	{
		//	std::cerr << "\n  Contents of globalToLocal map:" << std::endl;
		//	for (mapiter = globalToLocal.begin(); mapiter != globalToLocal.end(); ++mapiter)
		//		std::cerr << "\n    " << (int)(mapiter->first) << " (global) -> " << (int)(mapiter->second) << " (local)";
		//	std::cerr << std::endl;
		//	}

#	endif

	return new TipData(	stateListVec,												// stateListPosVec
						boost::shared_array<const int8_t>(tipSpecificStateCode),	// stateCodesShPtr
						num_rates,													// number of relative rate categories
						num_states,													// number of states in the model
						NULL,														// pMatTranspose
						true);														// managePMatrices
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates the InternalData data structure needed to store the conditional likelihood arrays and transition matrices
|	needed for likelihood calcuations. Returns a pointer to a newly-constructed InternalData structure. See the 
|	documentation for InternalData for more explanation.
*/
InternalData * TreeLikelihood::allocateInternalData()
	{
	return new InternalData(num_patterns,	// number of site patterns
							num_rates,		// number of relative rate categories
							num_states,		// number of model states
							NULL,		// pMat
							true);		// managePMatrices
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function that deletes the TipData structure allocated by the allocateTipData() function. The intention is for this
|	function to be given to TreeNode objects (in the form of a boost::function object) whenever allocateTipData() is 
|	used to allocate memory for a TipData structure. The TreeNode destructor can then use the boost::function object to
|	delete the memory required by the TipData structure without knowing any details of the TipData structure (i.e.,
|	the file that defines the TreeNode destructor does not need to include a header file that declares the details of
|	the TipData class.
*/
void deallocateTipData(TipData * p)
	{
	delete p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Function that deletes the InternalData structure allocated by the allocateInternalData() function. The intention is 
|	for this function to be given to TreeNode objects (in the form of a boost::function object) whenever 
|	allocateInternalData() is used to allocate memory for an InternalData structure. The TreeNode destructor can then 
|	use the boost::function object to delete the memory required by the InternalData structure without knowing any 
|	details of the InternalData structure (i.e., the file that defines the TreeNode destructor does not need to include 
|	a header file that declares the details of the InternalData class.
*/
void deallocateInternalData(InternalData * p)
	{
	delete p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a basic TipData object for all tip nodes and calls allocateInternalData() for all internal nodes in the 
|	supplied Tree. This function does not require a data matrix and does not add observed data to TipData structures.
|	Silently returns if root node already has a TipData structure. This is taken to indicate that prepareForLikelihood
|	or prepareForSimulation was previously called for this tree, in which case all data structures needed for 
|	simulation are already present.
*/
void TreeLikelihood::prepareForSimulation(
  TreeShPtr t)			/**< is the tree to decorate */
	{
	TreeNode::TipDataDeleter		td_deleter	= &deallocateTipData;
	TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;

	preorder_iterator nd = t->begin();

	// Only proceed if root node does not already have a TipData structure
	if (nd->GetTipData() != NULL)
		return;

	for (; nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			TipData * td = 	new TipData(num_rates,	num_states);	//@POL should be using shared_ptr here?
			nd->SetTipData(td, td_deleter);
			}
		else
			{
			InternalData * cl = allocateInternalData();
			nd->SetInternalData(cl, cl_deleter);
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Builds `pattern_map' and `pattern_counts' using uncompressed data stored in `mat'.
*/
void TreeLikelihood::copyDataFromDiscreteMatrix(
  const CipresNative::DiscreteMatrix & mat)		/**< is the data source */
	{
	nTaxa = mat.getNTax();

	// The compressDataMatrix function first erases, then builds, both pattern_map and 
	// pattern_counts using the uncompressed data contained in mat
	num_patterns = compressDataMatrix(mat);
	site_likelihoods.resize(num_patterns, 0.0);

	state_list = mat.getStateList(); 
	state_list_pos = mat.getStateListPos();

	// size of likelihood_rate_site vector needs to be revisited if the number of rates subsequently changes 
	recalcRelativeRates();
	//num_total_rates = model->getNRatesTotal();
	//likelihood_rate_site.resize(num_total_rates*num_patterns, 0.0);	//POL replaced mat.getNChar() with nPatterns
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Builds `pattern_map' and `pattern_counts' using data stored in `sim_data'.
*/
void TreeLikelihood::copyDataFromSimData(
  SimDataShPtr sim_data)	/**< is the data source */
	{
	// Copy simulated data to pattern_map
	pattern_map = sim_data->getSimPatternMap();

	// Build up counts vector
	pattern_counts.clear();
	for (PatternMapType::iterator it = pattern_map.begin(); it != pattern_map.end(); ++it)
		{
		pattern_counts.push_back(it->second);
		}

	nTaxa = sim_data->getPatternLength();
	num_patterns = (unsigned)pattern_map.size();
	site_likelihoods.resize(num_patterns, 0.0);

	model->buildStateList(state_list, state_list_pos);

	// size of likelihood_rate_site vector needs to be revisited if the number of rates subsequently changes 
	recalcRelativeRates();
	//num_total_rates = model->getNRatesTotal();
	//likelihood_rate_site.resize(num_total_rates*num_patterns, 0.0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds data currently stored in `pattern_map' to the patterns already in `other'. Assumes that `pattern_length' for 
|	this SimData object is identical to the `pattern_length' of `other'.
*/
void TreeLikelihood::addDataTo(SimData & other)
	{
	if (pattern_map.empty())
		return;
	if (other.getTotalCount() == 0)
		{
		assert(nTaxa > 0);
		other.resetPatternLength(nTaxa);
		}
	assert(nTaxa == other.getPatternLength());
	for (PatternMapType::iterator it = pattern_map.begin(); it != pattern_map.end(); ++it)
		{
		PatternCountType count = it->second;
		VecStateList & other_pattern = other.getCurrPattern();
		std::copy(it->first.begin(), it->first.end(), other_pattern.begin());
		other.insertPattern(count);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls allocateInternalData() to add an InternalData structure to `nd' containing the conditional likelihood arrays
|	needed for likelihood calculations.
*/
void TreeLikelihood::prepareInternalNodeForLikelihood(
  TreeNode * nd)	/**< is the node to decorate */
	{
	//InternalData * ndID = nd->GetInternalData();
	if (!nd)
		{
		TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;
		InternalData * cl = allocateInternalData();
		nd->SetInternalData(cl, cl_deleter);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls allocateTipData() for all tip nodes and allocateInternalData() for all internal nodes in the supplied Tree. 
|	Assumes each tip node number in the tree equals the appropriate row in the data matrix. 
*/
void TreeLikelihood::prepareForLikelihood(
  TreeShPtr t)									/**< is the tree to decorate */
	{
	// If no_data is true, it means that calcLnL will always return 0.0 immediately and 
	// will thus never need the TipData or InternalData data structures
	if (no_data)
		return;

	TreeNode::TipDataDeleter		td_deleter	= &deallocateTipData;
	TreeNode::InternalDataDeleter	cl_deleter	= &deallocateInternalData;

	//std::cerr << "In TreeLikelihood::prepareForLikelihood..." << std::endl; //POL temp

	for (preorder_iterator nd = t->begin(); nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			unsigned row = nd->GetNodeNumber();
			TipData * td = allocateTipData(row);
			nd->SetTipData(td, td_deleter);
			}
		else
			{
			InternalData * cl = allocateInternalData();
			nd->SetInternalData(cl, cl_deleter);
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copies data from `mat' to the map `pattern_map'. The resulting map holds pairs whose key is the pattern for one site
|	and whose value is a count of the number of sites having that pattern. The counts from `pattern_map' are transferred  
|	to the `pattern_counts' vector (vectors are more efficient containers for use during likelihood calculations).
*/
unsigned TreeLikelihood::compressDataMatrix(const CipresNative::DiscreteMatrix & mat) //POLBM TreeLikelihood::compressDataMatrix
	{
	pattern_map.clear();
	unsigned ntax = mat.getNTax();
	unsigned nchar = mat.getNChar();

	// Loop across each site in mat
	for (unsigned j = 0; j < nchar; ++j)
		{
		// Build up a vector representing the pattern of state codes at this site
		std::vector<int8_t> pattern;
		for (unsigned i = 0; i < ntax; ++i)
				{
			const int8_t * row  = mat.getRow(i);
			const int8_t   code = row[j];
			pattern.push_back(code);
			}

		// Add the pattern to the map if it has not yet been seen, otherwise increment 
		// the count of this pattern if it is already in the map (see item 24, p. 110, in Meyers' Efficient STL)
		PatternMapType::iterator lowb = pattern_map.lower_bound(pattern);
		if (lowb != pattern_map.end() && !(pattern_map.key_comp()(pattern, lowb->first)))
			{
			// pattern is already in pattern_map, increment count
			lowb->second += 1.0;
			}
		else
			{
			// pattern has not yet been stored in pattern_map
			pattern_map.insert(lowb, PatternMapType::value_type(pattern, 1));
			}
		}

	// Copy counts to pattern_counts before returning
	pattern_counts.clear();
	pattern_counts.reserve(pattern_map.size());
	for (PatternMapType::iterator mapit = pattern_map.begin(); mapit != pattern_map.end(); ++mapit)
		{
		pattern_counts.push_back(mapit->second);
		}

	return (unsigned)pattern_map.size();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a string representation of the supplied `state'. For example, if `state' equals 1 (where model is a 
|	standard DNA model), the string returned would be "C". If, however, `state' was 7 (again, standard DNA model), then
|	there is some ambiguity present and the string returned might look something like "{AC}" (the string returned 
|	depends of course on the actual meaning of the global state code 7).
*/
std::string TreeLikelihood::getStateStr(
  int8_t state)	const /**< is the global state code to be converted to a std::string */
	{
	std::string s;
	const int8_t nsPlusOne = num_states + 1;

	if (state < nsPlusOne)
		{
		// either no ambiguity or complete ambiguity
		s << model->lookupStateRepr((int)state);
		}
	else
		{
		// `state' represents partial ambiguity

		// First, find location of the definition of `state' in the global state list
		unsigned pos = state_list_pos[(unsigned)state];
		VecStateList::const_iterator it = state_list.begin() + pos;

		// Now get the number of basic states composing `state'
		unsigned n = *it++;

		// Walk down global state list converting states into strings
		s << "{";
		for (unsigned i = 0; i < n; ++i)
			{
			s << model->lookupStateRepr((int)*it++);
			}
		s << "}";
		}
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets private data member `nevals' to 0.
*/
void TreeLikelihood::resetNEvals()
	{
	nevals = 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor providing access to the current value of the private data member `nevals'.
*/
unsigned TreeLikelihood::getNEvals()
	{
	return nevals;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Assuming TreeLikelihood::compressDataMatrix has been called, so that `pattern_map' is up-to-date, returns a string 
|	listing all observed patterns and their frequencies.
*/
std::string TreeLikelihood::listPatterns(
  bool show_coded_states)	/**< if true, output global state codes used internally; otherwise, do the translation back to the original state representations */
	{
	std::vector<std::string> state_repr;
	std::string s;
	unsigned i = 0;
	PatternMapType::iterator it = pattern_map.begin();
	for (; it != pattern_map.end(); ++it, ++i)
		{
		const std::vector<int8_t> & p = it->first;
		PatternCountType c = it->second;
		s << str(boost::format("%6d %6.1f ") % i % c);
		unsigned ntax = (unsigned)p.size();
		if (show_coded_states)
			{
			for (unsigned j = 0; j < ntax; ++j)
				{
				s << str(boost::format("%d ") % (int)(p[j]));
				}
			}
		else
			{
			for (unsigned j = 0; j < ntax; ++j)
				{
				s << str(boost::format("%s") % getStateStr(p[j]));
				}
			}
		s << '\n';
		}
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	First erases `pattern_repr', then builds it into a vector of pattern representation strings. When the function 
|	returns, each element of `pattern_repr' is a string containing global state codes separated by whitespace. This 
|	function is intended to be used for debugging purposes, where it is sometimes helpful to be sure of exactly which
|	data pattern is being operated upon. Thus, no particular effort has been made to make this function efficient.
*/
void TreeLikelihood::buildPatternReprVector(std::vector<std::string> & pattern_repr, TreeShPtr t)
	{
	unsigned nTips		= t->GetNTips();
	//unsigned nStates	= getNStates();
	//unsigned nPatterns	= getNPatterns();

	pattern_repr.clear();
		pattern_repr.reserve(num_patterns);

	std::cerr << "\nglobal_pos:" << std::endl;
	for (unsigned z = 0; z < state_list_pos.size(); ++z)
		{
		std::cerr << str(boost::format("%3d") % z) << "  " << (int)state_list_pos[z] << std::endl;
		}

	// Recreate the data matrix in terms of tip-specific state codes
	int8_t * * m = NewTwoDArray<int8_t>(nTips, num_patterns);

	for (preorder_iterator nd = t->begin(); nd != t->end(); ++nd)
		{
		if (nd->IsTip())
			{
			// get node number i
			unsigned i = nd->GetNodeNumber();
			assert(i >= 0 && i < nTips);

			// get tip-specific state code array
			TipData * tipData = nd->GetTipData();
			assert(tipData);
			const int8_t * tipCodes = tipData->getConstStateCodes();

			// get position vector that allows us to translate tip-specific state codes into global state codes
			const std::vector<unsigned int> & global_statelist_pos = tipData->getConstStateListPos();

			std::cerr << "\nglobal_statelist_pos vector for node " << i << ": " << std::endl;
			for (unsigned z = 0; z < global_statelist_pos.size(); ++z)
				{
				std::cerr << "  " << (int)global_statelist_pos[z] << std::endl;
				}

			// fill row i of data matrix
			for (unsigned j = 0; j < num_patterns; ++j)
				{
				int8_t global_code = tipCodes[j];
				int8_t offset = global_code - ((int8_t)num_states + 1);
				if (offset >= 0)
					{
					unsigned pos = global_statelist_pos[offset];
					for (unsigned m = 0; m < (unsigned)state_list_pos.size(); ++m)
						{
						if (state_list_pos[m] == pos)
							global_code = (int8_t)(m);
						}
					}
				m[i][j] = global_code;

				//std::cerr << j << " | " << (int)tipCodes[j] << " | ";
				//if (offset >= 0)
				//	std::cerr << (int)global_statelist_pos[offset];
				//else
				//	std::cerr << "(" << (int)offset << ")";
				//std::cerr << std::endl;
				}
			}
		}

	for (unsigned j = 0; j < num_patterns; ++j)
		{
		std::string s;
		for (unsigned i = 0; i < nTips; ++i)
			{
			s << (int)m[i][j] << " ";
			}
		pattern_repr.push_back(s);
		}

	DeleteTwoDArray<int8_t>(m);
	}

}	// namespace phycas
