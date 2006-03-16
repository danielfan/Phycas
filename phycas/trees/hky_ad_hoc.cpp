#include "phycas/force_include.h"
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include "ncl/output/nxs_output.hpp"
#include "phycas/trees/tree_node.hpp"
#include "phycas/likelihood/tree_node_attributes.hpp"
#include "phycas/trees/tree.hpp"
#include "phycas/taxa/taxa_manager.hpp"
#include "phycas/characters/characters_manager.hpp"
#include "phycas/modules/mcmc/gibbs_param.hpp"
#include "phycas/trees/hky_ad_hoc.hpp"
#include "phycas/misc/compressed_2d_array.hpp"
#include "phycas/trees/underflow.hpp"

using std::string;
#define A		0
#define C		1
#define G		2
#define T		3
#define MISSING	4
#define TOTAL	4

#define IS_MISSING(c) (c == '?' || c == 'N' || c == '-' || c == 'R' || c == 'Y')
#define IS_NUCLEOTIDE(c) (c == 'A' || c == 'a' || c == 'C' || c == 'c' || c == 'G' || c == 'g' || c == 'T' || c == 't')
#define UINT_TO_NUCLEOTIDE(u) ((u)==0 ? 'A' : ((u)==1 ? 'C' : ((u)==2 ? 'C' : 'T')))

#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::fabs;
	using std::exp;
	using std::log;
	using std::exit;
#endif

//	Summary of Mark's code:
//	o AnalysisManager::CreateScoreableTree called with tree description and scoring description
//		- The tree description is actually an NCL FullTreeDescription
//		- The scoring description is a vector<SubsetModelDescription>
//		- A SubsetModelDescription is a vector<ModelDescription> (I assume usually just one element in this vector)
//		- A ModelDescription comprises a UniqueModelTypeID (kJCModelID, kK2PModelID, kHKYModelID, or kGTRModelID), 
//		and an embedded MultiSiteModelDescription that knows whether the model is a site-specific rates model,
//		an invariant sites model, etc.
//	o ScoreableTree constructor called, passed a Tree (created from the FullTreeDescription) and a 
//		ScoreableTreeInitializer (created from the vector<SubsetModelDescription>)
//		- the ScoreableTreeInitializer comprises a vector<SubsetScorerInitializer> and two EdgeLenAccumulatorShPtr
//		- EdgeLenAccumulator changes the edge length for all subsets in the TreeNode attr objects
//	o ScoreableTree::CalculateLnLikelihood function called to compute likelihood of tree

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets all pointers to NULL, initializes `ntaxa', `nchar' and `npatterns' to 0, `kappa' to 1.0 and all 
|	relative base frequencies in `pi' to 0.25.
*/
HKYAdHocEvaluator::HKYAdHocEvaluator(
  GibbsParameterShPtr freqAParam,			/**< is a shared pointer to the dirichlet parameter governing freq. A */
  GibbsParameterShPtr freqCParam,			/**< is a shared pointer to the dirichlet parameter governing freq. C */
  GibbsParameterShPtr freqGParam,			/**< is a shared pointer to the dirichlet parameter governing freq. G */
  GibbsParameterShPtr freqTParam,			/**< is a shared pointer to the dirichlet parameter governing freq. T */
  GibbsParameterShPtr kappaParam,			/**< is a shared pointer to the transition/transversion rate ratio parameter to use */
  GibbsParameterShPtr rateVarParam,			/**< is a shared pointer to the discrete gamma variance parameter to use */
  GibbsParameterShPtr pInvarParam,			/**< is a shared pointer to the proportion invariable sites parameter to use */
  GibbsParameterShPtr meanRatesParam,		/**< is a shared pointer to the mean rates parameter to use */
  unsigned ncateg)							/**< is the number of rate categories */
  : PhoRateManager(pInvarParam, rateVarParam, meanRatesParam, ncateg), modelParams(HKYAdHocEvaluator::kHKYLastParam)
	{
	tree				= NULL;
	ntaxa				= 0;
	nchar				= 0;
	npatterns			= 0;
	kappa				= 1.0;
	rateVar				= 2.0;
	ncat				= ncateg;
	paramOutOfBounds	= false;
	entireTreeDirty		= false;
	pi[0]				= 0.25;
	pi[1]				= 0.25;
	pi[2]				= 0.25;
	pi[3]				= 0.25;

#if defined(POL_UNDERFLOW_CORRECTION)
	underflow.reset(new UnderflowManager());
	//underflow = new UnderflowManager();
#endif

#if !defined(USING_VARIANCE_FOR_RATEHET)
#	error Must define USING_VARIANCE_FOR_RATEHET or modify HKYAdHocEvaluator accordingly
#endif

	relrate = GetRepRates();

	ResetParameterPtr(kHKYFreqAIndex, freqAParam);
	ResetParameterPtr(kHKYFreqCIndex, freqCParam);
	ResetParameterPtr(kHKYFreqGIndex, freqGParam);
	ResetParameterPtr(kHKYFreqTIndex, freqTParam);
	ResetParameterPtr(kHKYKappaIndex, kappaParam);
	ResetParameterPtr(kHKYRateVarIndex, rateVarParam);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Counts the number of elements in the `modelParam' array that point to real objects and thus can be updated during a
|	Gibbs sampling step. Note that some parameters that are updatable are not considered free parameters in the 
|	statistical sense (e.g. all four base frequency dirichlet parameters are updated, but only three of the four
|	base frequencies used in computing likelihoods are free because each of the dirichlet parameters is divided by the
|	sum of all four to normlize them, making it possible to calculate one of the base frequencies by subtraction.
*/
unsigned HKYAdHocEvaluator::GetNumUpdatableParameters()
	{
	unsigned np = 0;
	for (unsigned i = 0; i < HKYAdHocEvaluator::kHKYLastParam; ++i)
		{
		if (modelParams[i])
			++np;
		}
	return np;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Produces a table of updateable parameters similar to the following on the output stream `out':
|>
|	Parameter    Fixed  Current Value
|	---------------------------------
|	Kappa        no     4.0
|	RateVar      no     0.5
|	FreqA        yes    0.25680
|	FreqC        yes    0.21478
|	FreqG        yes    0.29092
|	FreqT        yes    0.23750
|	---------------------------------
|	Average of underlying frequency parameters: 1.103425
|	Number of rate categories: 4
|<
|	In this case, the frequencies are fixed at their empirical values, while Kappa and RateVar will be updated by
|	Gibbs sampling during the MCMC run.
*/
void HKYAdHocEvaluator::ShowUpdatableParameters(
  NxsOutputStream & out)	const /**< is the output stream on which to display the table of parameters */
	{
	//@POL this is admittedly pretty crude at the moment; note especially that no parameters are ever 
	// fixed at this writing because there is currently no mechanism for the user to choose to fix a
	// parameter's value

	bool include_kappa   = modelParams[kHKYKappaIndex];
	bool include_ratevar = modelParams[kHKYRateVarIndex];
	bool include_freqs   = modelParams[kHKYFreqAIndex];
	assert(!include_freqs || modelParams[kHKYFreqCIndex]);
	assert(!include_freqs || modelParams[kHKYFreqGIndex]);
	assert(!include_freqs || modelParams[kHKYFreqTIndex]);

	out << "Parameter    Fixed  Current Value\n";
	out << "---------------------------------\n";
	if (include_kappa)
		out << "Kappa        no     " <<  modelParams[kHKYKappaIndex]->GetValue() << '\n';
	if (include_ratevar)
		out << "RateVar      no     " << modelParams[kHKYRateVarIndex]->GetValue() << '\n';
	
	double cA = 0.0;
	double cC = 0.0;
	double cG = 0.0;
	double cT = 0.0;
	double cTotal = 0.0;

	if (include_freqs)
		{
		cA = modelParams[kHKYFreqAIndex]->GetValue();
		cC = modelParams[kHKYFreqCIndex]->GetValue();
		cG = modelParams[kHKYFreqGIndex]->GetValue();
		cT = modelParams[kHKYFreqTIndex]->GetValue();
		assert(cA >= 0.0);
		assert(cC >= 0.0);
		assert(cG >= 0.0);
		assert(cT >= 0.0);
		cTotal = cA + cC + cG + cT;
		out << "FreqA        no     " << (cA/cTotal) << '\n';
		out << "FreqC        no     " << (cC/cTotal) << '\n';
		out << "FreqG        no     " << (cG/cTotal) << '\n';
		out << "FreqT        no     " << (cT/cTotal) << '\n';
		}

	out << "---------------------------------\n\n" ;
	out << "Average of underlying frequency parameters: " << ((cA + cC + cG + cT)/4.0) << '\n';
	out <<  "Number of rate categories: " << ncat << '\n';
	if (ncat > 1)
		{
		out << "Representative relative rates:\n";
		std::vector<double> ratevec = GetRepRates();
		for (unsigned rr = 0; rr < ratevec.size(); ++rr)
			out << "  rate " << rr << " has mean " << ratevec[rr] << '\n';
		}

	out << ncl::endl;
 	}

/*----------------------------------------------------------------------------------------------------------------------
|	Refreshes local copies of parameter values from their corresponding GibbsParameter objects, then calls 
|	InvalidateTree() because a change in parameter means that the likelihood for the entire tree must be recomputed the 
|	next time a likelihood is needed.
*/
void HKYAdHocEvaluator::ParameterChanged()
	{
	paramOutOfBounds = false;
	kappa		= modelParams[kHKYKappaIndex]->GetValue();
	if (kappa <= 0.0)
		paramOutOfBounds = true;

	rateVar		= modelParams[kHKYRateVarIndex]->GetValue();
	if (rateVar > 0.0)
		{
		// CalculateRates() and GetRepRates() are PhoRateManager member functions that operate with the 
		// GibbsParameter object pointed to by modelParams[kHKYRateVarIndex] (passed to PhoRateManager
		// via the constructor
		//
		relrate.clear();
		CalculateRates();
		relrate = GetRepRates();
		}
	else
		{
		paramOutOfBounds = true;
		}

	double cA	= modelParams[kHKYFreqAIndex]->GetValue();
	if (cA < 0.0)
		paramOutOfBounds = true;

	double cC	= modelParams[kHKYFreqCIndex]->GetValue();
	if (cC < 0.0)
		paramOutOfBounds = true;

	double cG	= modelParams[kHKYFreqGIndex]->GetValue();
	if (cG < 0.0)
		paramOutOfBounds = true;

	double cT	= modelParams[kHKYFreqTIndex]->GetValue();
	if (cT < 0.0)
		paramOutOfBounds = true;

#if 1
	// This overwrites the empirical base frequencies
	//
	if (!paramOutOfBounds)
		{
		double total = cA + cC + cG + cT;
		assert(total > 0.0);
		pi[A] = cA/total;
		pi[C] = cC/total;
		pi[G] = cG/total;
		pi[T] = cT/total;
		}
#endif

	InvalidateTree();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copies data from the supplied character manager `cm' and copies taxon labels from the supplied taxa manager `tm'.
*/
void HKYAdHocEvaluator::CopyData(
  PhoTaxaManager &tm,			/**< is the taxa manager supplying taxon labels */
  PhoCharactersManager &cm)		/**< is the character manager supplying data */
	{
	unsigned i, k;
	string tmp, taxonName;

	unsigned baseCounts[5];
	baseCounts[A] = baseCounts[C] = baseCounts[G] = baseCounts[T] = baseCounts[TOTAL] = 0;

	Clear();

	// Store indices of active taxa
	NxsIndexSet activeTaxa = tm.GetActiveTaxa();

	// Define the data type
	NxsDataType dnaType(NxsDataType::kDNA, false);

	// Get compressed version of data matrix, obviates need to call CrunchData
	matrixShPtr = cm.GetCompressedMatrix(dnaType, activeTaxa);

	// Get dimensions of data matrix
	npatterns = matrixShPtr->GetNumPatterns();
	ntaxa = matrixShPtr->GetNumTaxa();
	unsigned nstates = matrixShPtr->GetDataType().GetNumStates();

	//patfrq = new unsigned[npatterns];
	patfrq.reset(new unsigned[npatterns]);

	// Count number of characters and assign values to patfrq array
	//
	std::vector<double> wts = matrixShPtr->GetPatternWeights();
	nchar = 0;
	for (unsigned  j = 0; j < npatterns; ++j)
		{
		double wt = wts[j];
		patfrq[(std::ptrdiff_t) j] = (unsigned)wt;
		nchar += patfrq[(std::ptrdiff_t) j];
		if (fabs(wt - patfrq[(std::ptrdiff_t) j]) > 1.e-8)
			{
			tmp.clear();
			StrPrintF(tmp, "Warning: pattern frequency %f is not an integer", wt);
			//@@@ std::cerr << tmp << std::endl;
			}
		}

	// Allocate pattern shared array, which will hold the data for the tip nodes
	//
	assert(!pattern);
	pattern = DataMatrix(new StateArray[ntaxa]);

	// Copy and store patterns
	//
	for (i = 0; i < ntaxa; ++i)
		{
		taxonName = tm.GetLabel(i);
		taxonLabels.push_back(taxonName);

		// OrdCodedArrs is derived from Compressed2DArray, differing mainly in having a specialized 
		// const_iterator that can handle ambiguity in state for a particular taxon. Here is how 
		// characters are stored in an OrdCodeArrs object:
		//
		// States for 8 taxa for one pattern: 0 3 1 1 {2,3} 0 {0,1,2,3} 2
		//   (note ambiguity for taxa 5 and 7)
		//
		// As stored in Compressed2DArray:    0 3 1 1 -2 2 3 0 -4 2
		//   -2 indicates ambiguity/polymorphism and next 2 states are for same taxon
		//   -4 indicates completely missing data and there is no need to list states afterwards
		//
		// The OrdCodedArrs::const_iterator would return the following vectors for each of the 
		// eight taxa if the array is walked using DerefThenAdvance:
		//
		//   (0)    0 is the only state for taxon 0
		//   (3)    3 is the only state for taxon 1
		//   (1)    1 is the only state for taxon 2
		//   (1)    1 is the only state for taxon 3
		//   (2,3)  taxon 4 is ambiguous, either state 2 or state 3 is present
		//   (0)    0 is the only state for taxon 5
		//   ()     taxon 6 has completely missing data, so empty vector returned
		//   (2)    2 is the only state for taxon 7
		//
		OrdCodedArrs thisTaxOCA(nstates);
		matrixShPtr->GetTaxonsOrdCodedStates(&thisTaxOCA, i);

		StateArray &states = pattern[(std::ptrdiff_t) i];

		OrdCodedArrs::const_iterator oIt = thisTaxOCA.begin();

		char xval;
		char all_missing = (char)(-4);
		for (unsigned j = 0; j < npatterns; ++j)
			{
			std::vector<unsigned> x = oIt.DerefThenAdvance();
			unsigned sz = (unsigned)x.size();
			if (sz == 0)
				{
				// Vector x will be empty if there is ambiguity for taxon i
				// in pattern j. Push back -4 to indicate this.
				// 
				states.push_back(all_missing);
				}
			else if (sz == 1)
				{
				// Vector x will contain only one element if there is no ambiguity for taxon i
				// in pattern j. Push back the state (0, 1, 2 or 3).
				// 
				xval = (char)x[0];
				assert(xval >= 0);
				assert(xval <= 3);
				states.push_back(xval);
				baseCounts[(int) xval]  += (unsigned)wts[j];
				baseCounts[TOTAL] += (unsigned)wts[j];
				}
			else
				{
				// Vector x will contain more than one element if there is ambiguity for taxon i
				// in pattern j. First push back minus the number of states observed, and then
				// push back each state (0, 1, 2 or 3) observed.
				//
				xval = static_cast<char>(-static_cast<char>(sz)); //POL minus sign should go outside, otherwise VC complains "unary minus operator applied to unsigned type, result still unsigned"
																	// mwerks converts "- char" to a negative int and then complains about "int" to "char" cast
				states.push_back(xval);
				for (k = 0; k < sz; ++k)
					{
					xval = (char)x[k];
					assert(xval >= 0);
					assert(xval <= 3);
					states.push_back(xval);
					baseCounts[(int) xval]  += (unsigned)wts[j];
					baseCounts[TOTAL] += (unsigned)wts[j];
					}
				}
			}	// loop over patterns
		}	// loop over taxa

	assert(baseCounts[TOTAL] > 0);
	pi[A] = (double)baseCounts[A] / (double)baseCounts[TOTAL];
	pi[C] = (double)baseCounts[C] / (double)baseCounts[TOTAL];
	pi[G] = (double)baseCounts[G] / (double)baseCounts[TOTAL];
	pi[T] = (double)baseCounts[T] / (double)baseCounts[TOTAL];

#if 0
	tmp.clear();
	StrPrintF(tmp, "Total number of bases counted = %d\nEmpirical nucleotide relative frequencies:\n  A %12.5f\n  C %12.5f\n  G %12.5f\n  T %12.5f"
		, baseCounts[TOTAL], pi[A], pi[C], pi[G], pi[T]);
	std::ofstream freqf("frequencies.txt");
	freqf << tmp << std::endl;
	freqf.close();
#endif
	}	// CopyData

/*----------------------------------------------------------------------------------------------------------------------
|	This deletes all dynamically allocated memory used by this object and returns the object to the just constructed
|	state. Called by the destructor.
*/
void HKYAdHocEvaluator::Clear()
	{
	taxonLabels.clear();
	patternStr.clear();

	pattern.reset();
	patfrq.reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Debugging function that outputs a table on `outf' showing the conditional likelihood array for the node `nd'.
*/
void HKYAdHocEvaluator::DebugShowCondLikeArrays(
  TreeNode *nd)			/**< is the node to display */
	{
	TreeNodeAttribute *attrib = (nd->attr)[0];

	outf << "\nCond. like. array for " << (nd->IsRoot() ? "root" : (nd->IsShootTip() ? "tip" : "internal")) << " node " << nd->GetNodeNumber() << ", edge length = " << nd->GetFltEdgeLen() << std::endl;
	outf << std::setw(6) << 'i';
	outf << std::setw(6) << "count";
	outf << std::setw(20) << 'A';
	outf << std::setw(20) << 'C';
	outf << std::setw(20) << 'G';
	outf << std::setw(20) << 'T';
	outf << std::endl;

	double *cl = attrib->condLike;
	for (unsigned i = 0; i < npatterns; ++i)
		{
		outf << std::setw(6) << i;
		outf << std::setw(6) << patfrq[(std::ptrdiff_t) i];
		outf << std::setw(20) << std::setprecision(15) << *cl++;
		outf << std::setw(20) << std::setprecision(15) << *cl++;
		outf << std::setw(20) << std::setprecision(15) << *cl++;
		outf << std::setw(20) << std::setprecision(15) << *cl++;
		outf << std::endl;
		}
	outf << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes the conditional likelihood array for node `nd'. Assumes `prMatrices' is up to date for all nodes. 
|	Assumes `nd' is not the root node. The conditional likelihood arrays at interior nodes are laid out as follows:
|>
|	 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+-
|	 | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T |
|	 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+-
|	 |    rate 1     |    rate 2     |    rate 3     |    rate 4     |    rate 1     |    rate 2     |    rate 3     |
|	 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+-
|	 |                          pattern 1                            |                         pattern 2              
|	 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+-
|<
|	Important: this function does not watch for or correct for underflow. Use UFProtectedCondLike member function if
|	underflow has been detected. If changes are made to this function, the corresponding change should be made in
|	UFProtectedCondLike.
*/
void HKYAdHocEvaluator::RecalcCondLike(
  TreeNode *nd)					/**< is the node for which the conditional likelihood array is to be recomputed */
	{
	TreeNodeAttribute *attrib = GetNodeAttr(nd);
	attrib->clDirty = false;

	// For each child, loop through this node's entire condLike array. For first child, 
	//   use operator =, otherwise, use operator *=.
	// Note that nd will not be a tip node, tips are avoided for purposes of
	//   recalculating conditional likelihood arrays.
	// If child is a tip, child's contribution is simply a sum of transition probabilities
	//   pr[i][j], where i is nd's state and j is one of the states present in the child.
	// If child is an internal node, child's contribution is sum over all states of 
	//   transition probability times conditional likelihood at the child node.
	//
	for (TreeNode *child = nd->GetLeftChild(); child != NULL; child = child->GetRightSib())
		{
		TreeNodeAttribute *child_attr = GetNodeAttr(child);
		assert(!child_attr->pMatIsDirty);
		double fromA, fromC, fromG, fromT;

		// cpr is the array of transition probability matrices for child's edge (one 4x4 matrix
		// for each relative rate category)
		//
		double ***cpr = child_attr->prMatrices;

		// cl is nd's conditional likelihood array
		//
		double *cl = attrib->condLike;

		if (child->IsShootTip())
			{
			unsigned cnn = child->GetNodeNumber();

			// states is an array of (ord. coded) states observed in child
			// note: states array is not necessarily npatterns long, so do not be tempted to replace 
			// len with npatterns in loop
			//
			StateArray states = (StateArray) pattern[(std::ptrdiff_t) cnn];	//@#POL lots of L2 cache misses here
			unsigned len = (unsigned)states.size();
			unsigned pattern = 0;

			for (unsigned i = 0; i < len; ++i)
				{
				for (unsigned r = 0; r < ncat; ++r)
					{
					fromA = 0.0;
					fromC = 0.0;
					fromG = 0.0;
					fromT = 0.0;
					char s = states[i]; //@#POL precalculate states[i] before this loop
					if (s == -4)
						{
						// All missing data does not contribute anything to conditional likelihood.
						// Effectively, the contribution of this child to cl[A] would be
						//   (cpr[r][A][A])(1) + (cpr[r][A][C])(1) + (cpr[r][A][G])(1) + (cpr[r][A][T])(1) = 1.0
						// so there is no need to bother. The (1) terms above represent the 
						// conditional likelihood of the child node, which is 1 for every state
						// in this case.
						//
						fromA = 1.0;
						fromC = 1.0;
						fromG = 1.0;
						fromT = 1.0;
						}
					else if (s < 0)
						{
						// Ambiguity present but not complete missing data. For every state j
						// found in child, add transition probability pr[i][j] to fromi.
						// For example, if A and G are both present in child for this pattern, this
						// would be the result (effectively):
						//   fromA = (cpr[r][A][A])(1) + (cpr[r][A][C])(0) + (cpr[r][A][G])(1) + (cpr[r][A][T])(0)
						//   fromC = (cpr[r][C][A])(1) + (cpr[r][C][C])(0) + (cpr[r][C][G])(1) + (cpr[r][C][T])(0)
						//   fromG = (cpr[r][G][A])(1) + (cpr[r][G][C])(0) + (cpr[r][G][G])(1) + (cpr[r][G][T])(0)
						//   fromT = (cpr[r][T][A])(1) + (cpr[r][T][C])(0) + (cpr[r][T][G])(1) + (cpr[r][T][T])(0)
						// Again, the (1) and (0) terms represent the conditional likelihood of the child.
						//
						unsigned ns = (unsigned)(-s);
						assert(i + ns < len);
						for (unsigned z = 0; z < ns; ++z)
							{
							++i;
							char ss = states[i];
							assert(ss == A || ss == C || ss == G || ss == T);
							fromA += cpr[r][A][ss]; //@#POL the cpr array not laid out well, ACGT should be cpr[ss][r][ACGT]
							fromC += cpr[r][C][ss];
							fromG += cpr[r][G][ss];
							fromT += cpr[r][T][ss];
							}
						}
					else
						{
						// No ambiguity present; s should already equal the state of the child for this pattern
						//
						assert(s == A || s == C || s == G || s == T);
						fromA += cpr[r][A][s];
						fromC += cpr[r][C][s];
						fromG += cpr[r][G][s];
						fromT += cpr[r][T][s];
						}

					// Save the child's contribution to this pattern and advance conditional likelihood 
					// array pointer to next pattern
					//
					if (child == nd->GetLeftChild())
						{
						*cl++ = fromA;
						*cl++ = fromC;
						*cl++ = fromG;
						*cl++ = fromT;
						}
					else
						{
						*cl++ *= fromA;
						*cl++ *= fromC;
						*cl++ *= fromG;
						*cl++ *= fromT;
						}
					}	// loop over rate categories
				pattern++;
				}	// loop over patterns
			}	// child->IsShootTip()
		else
			{
			// child is an internal node
			//
			double *ccl = child_attr->condLike;

			for (unsigned pattern = 0; pattern < npatterns; ++pattern)
				{
				for (unsigned r = 0; r < ncat; ++r)
					{
					//@#POL need to unroll the operations below so that not so much array indexing is going on
					fromA = (cpr[r][A][A]*ccl[A] + cpr[r][A][C]*ccl[C] + cpr[r][A][G]*ccl[G] + cpr[r][A][T]*ccl[T]);
					fromC = (cpr[r][C][A]*ccl[A] + cpr[r][C][C]*ccl[C] + cpr[r][C][G]*ccl[G] + cpr[r][C][T]*ccl[T]);
					fromG = (cpr[r][G][A]*ccl[A] + cpr[r][G][C]*ccl[C] + cpr[r][G][G]*ccl[G] + cpr[r][G][T]*ccl[T]);
					fromT = (cpr[r][T][A]*ccl[A] + cpr[r][T][C]*ccl[C] + cpr[r][T][G]*ccl[G] + cpr[r][T][T]*ccl[T]);
					
					// Advance child's conditional likelihood array pointer to set of four
					//
					ccl += 4;

					// Save the child's contribution to this pattern and advance conditional likelihood 
					// array pointer to next set of four
					//
					if (child == nd->GetLeftChild())
						{
						*cl++ = fromA;
						*cl++ = fromC;
						*cl++ = fromG;
						*cl++ = fromT;
						}
					else
						{
						*cl++ *= fromA;
						*cl++ *= fromC;
						*cl++ *= fromG;
						*cl++ *= fromT;
						}
					}	// loop over rate categories
				}	// loop over patterns
			}	// else child is internal node
		}	// loop over children
	}

#if defined(POL_UNDERFLOW_CORRECTION)
/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes the conditional likelihood array for node `nd'. Assumes `prMatrices' is up to date for all nodes. 
|	Assumes `nd' is not the root node. The conditional likelihood arrays at interior nodes are laid out as follows:
|>
|	 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+-
|	 | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T |
|	 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+-
|	 |    rate 1     |    rate 2     |    rate 3     |    rate 4     |    rate 1     |    rate 2     |    rate 3     |
|	 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+-
|	 |                          pattern 1                            |                         pattern 2              
|	 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+-
|<
|	This is a version of RecalcCondLike that checks for underflow. This function should only be used to recalculate 
|	conditional likelihood arrays if underflow has been detected (it is not as fast as RecalcCondLike).
*/
void HKYAdHocEvaluator::UFProtectedCondLike(
  TreeNode *nd)					/**< is the node for which the conditional likelihood array is to be recomputed */
	{
	TreeNodeAttribute *attrib = GetNodeAttr(nd);
	attrib->clDirty = false;

	unsigned z;
	unsigned last = 4*ncat;
	double uf = underflow->GetLikeCutoff();
	
	// For each child, loop through this node's entire condLike array. For first child, 
	//   use operator =, otherwise, use operator *=.
	// Note that nd will not be a tip node, tips are avoided for purposes of
	//   recalculating conditional likelihood arrays.
	// If child is a tip, child's contribution is simply a sum of transition probabilities
	//   pr[i][j], where i is nd's state and j is one of the states present in the child.
	// If child is an internal node, child's contribution is sum over all states of 
	//   transition probability times conditional likelihood at the child node.
	//
	for (TreeNode *child = nd->GetLeftChild(); child != NULL; child = child->GetRightSib())
		{
		TreeNodeAttribute *child_attr = GetNodeAttr(child);
		assert(!child_attr->pMatIsDirty);
		double fromA, fromC, fromG, fromT;

		// cpr is the array of transition probability matrices for child's edge (one 4x4 matrix
		// for each relative rate category)
		//
		double ***cpr = child_attr->prMatrices;

		// cl is nd's conditional likelihood array
		//
		double *cl = attrib->condLike;

		if (child->IsShootTip())
			{
			unsigned cnn = child->GetNodeNumber();

			// states is an array of (ord. coded) states observed in child
			// note: states array is not necessarily npatterns long, so do not be tempted to replace 
			// len with npatterns in loop
			//
			StateArray states = (StateArray) pattern[(std::ptrdiff_t) cnn];
			unsigned len = (unsigned)states.size();
			unsigned pattern = 0;
			unsigned *bounces = underflow->GetBouncesArray();

			for (unsigned i = 0; i < len; ++i)
				{
				unsigned num_uf_corrections = 0;
				//bool has_underflowed = underflow->HasUnderflowed(pattern);
				bool has_underflowed = (*bounces++ < UINT_MAX);
				double *cl_start = cl;

				for (unsigned r = 0; r < ncat; ++r)
					{
					fromA = 0.0;
					fromC = 0.0;
					fromG = 0.0;
					fromT = 0.0;
					char s = states[i];
					if (s == -4)
						{
						// All missing data does not contribute anything to conditional likelihood.
						// Effectively, the contribution of this child to cl[A] would be
						//   (cpr[r][A][A])(1) + (cpr[r][A][C])(1) + (cpr[r][A][G])(1) + (cpr[r][A][T])(1) = 1.0
						// so there is no need to bother. The (1) terms above represent the 
						// conditional likelihood of the child node, which is 1 for every state
						// in this case.
						//
						fromA = 1.0;
						fromC = 1.0;
						fromG = 1.0;
						fromT = 1.0;
						}
					else if (s < 0)
						{
						// Ambiguity present but not complete missing data. For every state j
						// found in child, add transition probability pr[i][j] to fromi.
						// For example, if A and G are both present in child for this pattern, this
						// would be the result (effectively):
						//   fromA = (cpr[r][A][A])(1) + (cpr[r][A][C])(0) + (cpr[r][A][G])(1) + (cpr[r][A][T])(0)
						//   fromC = (cpr[r][C][A])(1) + (cpr[r][C][C])(0) + (cpr[r][C][G])(1) + (cpr[r][C][T])(0)
						//   fromG = (cpr[r][G][A])(1) + (cpr[r][G][C])(0) + (cpr[r][G][G])(1) + (cpr[r][G][T])(0)
						//   fromT = (cpr[r][T][A])(1) + (cpr[r][T][C])(0) + (cpr[r][T][G])(1) + (cpr[r][T][T])(0)
						// Again, the (1) and (0) terms represent the conditional likelihood of the child.
						//
						unsigned ns = (unsigned)(-s);
						assert(i + ns < len);
						for (unsigned z = 0; z < ns; ++z)
							{
							++i;
							char ss = states[i];
							assert(ss == A || ss == C || ss == G || ss == T);
							fromA += cpr[r][A][ss];
							fromC += cpr[r][C][ss];
							fromG += cpr[r][G][ss];
							fromT += cpr[r][T][ss];
							}
						}
					else
						{
						// No ambiguity present; s should already equal the state of the child for this pattern
						//
						assert(s == A || s == C || s == G || s == T);
						fromA += cpr[r][A][s];
						fromC += cpr[r][C][s];
						fromG += cpr[r][G][s];
						fromT += cpr[r][T][s];
						}

					// Save the child's contribution to this pattern and advance conditional likelihood 
					// array pointer to next pattern
					//
					if (has_underflowed)
						{
						if (child == nd->GetLeftChild())
							{
							cl[A] = fromA;
							cl[C] = fromC;
							cl[G] = fromG;
							cl[T] = fromT;
							}
						else
							{
							cl[A] *= fromA;
							cl[C] *= fromC;
							cl[G] *= fromG;
							cl[T] *= fromT;
							}
						if (cl[A] < uf)
							{
							z = underflow->CalcNumCorrectionsNeeded(cl[A]);
							if (z > num_uf_corrections)
								num_uf_corrections = z;
							}
						if (cl[C] < uf)
							{
							z = underflow->CalcNumCorrectionsNeeded(cl[C]);
							if (z > num_uf_corrections)
								num_uf_corrections = z;
							}
						if (cl[G] < uf)
							{
							z = underflow->CalcNumCorrectionsNeeded(cl[G]);
							if (z > num_uf_corrections)
								num_uf_corrections = z;
							}
						if (cl[T] < uf)
							{
							z = underflow->CalcNumCorrectionsNeeded(cl[T]);
							if (z > num_uf_corrections)
								num_uf_corrections = z;
							}
						cl += 4;
						}
					else	// if (has_underflowed)...
						{
						if (child == nd->GetLeftChild())
							{
							*cl++ = fromA;
							*cl++ = fromC;
							*cl++ = fromG;
							*cl++ = fromT;
							}
						else
							{
							*cl++ *= fromA;
							*cl++ *= fromC;
							*cl++ *= fromG;
							*cl++ *= fromT;
							}
						}	// if (has_underflowed) else ...
					}	// loop over rate categories

					if (num_uf_corrections > 0)
						{
						// Need to go back and bounce up each element of the conditional likelihood
						// array that pertains to this pattern
						//
						cl = cl_start;
						double f = underflow->CalcCorrectionFactor(num_uf_corrections);
						for (z = 0; z < last; z++)
							*cl++ *= f;

						// Keep track of the number of times this pattern has been bounced so that 
						// the correction factors can be removed after the log of the site likelihood
						// has been computed
						//
						underflow->BounceSite(pattern, num_uf_corrections);
						}
				pattern++;
				}	// loop over patterns
			}	// child->IsShootTip()
		else
			{
			// child is an internal node
			//
			double *ccl = child_attr->condLike;
			unsigned *bounces = underflow->GetBouncesArray();

			for (unsigned pattern = 0; pattern < npatterns; ++pattern)
				{
				unsigned num_uf_corrections = 0;
				//bool has_underflowed = underflow->HasUnderflowed(pattern);
				bool has_underflowed = (*bounces++ < UINT_MAX);
				double *cl_start = cl;

				for (unsigned r = 0; r < ncat; ++r)
					{
					fromA = (cpr[r][A][A]*ccl[A] + cpr[r][A][C]*ccl[C] + cpr[r][A][G]*ccl[G] + cpr[r][A][T]*ccl[T]);
					fromC = (cpr[r][C][A]*ccl[A] + cpr[r][C][C]*ccl[C] + cpr[r][C][G]*ccl[G] + cpr[r][C][T]*ccl[T]);
					fromG = (cpr[r][G][A]*ccl[A] + cpr[r][G][C]*ccl[C] + cpr[r][G][G]*ccl[G] + cpr[r][G][T]*ccl[T]);
					fromT = (cpr[r][T][A]*ccl[A] + cpr[r][T][C]*ccl[C] + cpr[r][T][G]*ccl[G] + cpr[r][T][T]*ccl[T]);
					
					// Advance child's conditional likelihood array pointer to set of four
					//
					ccl += 4;

					// Save the child's contribution to this pattern and advance conditional likelihood 
					// array pointer to next set of four
					//
					if (has_underflowed)
						{
						if (child == nd->GetLeftChild())
							{
							cl[A] = fromA;
							cl[C] = fromC;
							cl[G] = fromG;
							cl[T] = fromT;
							}
						else
							{
							cl[A] *= fromA;
							cl[C] *= fromC;
							cl[G] *= fromG;
							cl[T] *= fromT;
							}
						if (cl[A] < uf)
							{
							z = underflow->CalcNumCorrectionsNeeded(cl[A]);
							if (z > num_uf_corrections)
								num_uf_corrections = z;
							}
						if (cl[C] < uf)
							{
							z = underflow->CalcNumCorrectionsNeeded(cl[C]);
							if (z > num_uf_corrections)
								num_uf_corrections = z;
							}
						if (cl[G] < uf)
							{
							z = underflow->CalcNumCorrectionsNeeded(cl[G]);
							if (z > num_uf_corrections)
								num_uf_corrections = z;
							}
						if (cl[T] < uf)
							{
							z = underflow->CalcNumCorrectionsNeeded(cl[T]);
							if (z > num_uf_corrections)
								num_uf_corrections = z;
							}
						cl += 4;
						}
					else	// if (has_underflowed) ...
						{
						if (child == nd->GetLeftChild())
							{
							*cl++ = fromA;
							*cl++ = fromC;
							*cl++ = fromG;
							*cl++ = fromT;
							}
						else
							{
							*cl++ *= fromA;
							*cl++ *= fromC;
							*cl++ *= fromG;
							*cl++ *= fromT;
							}
						}	// if (has_underflowed) ... else ...
					}	// loop over rate categories

					if (num_uf_corrections > 0)
						{
						// Need to go back and bounce up each element of the conditional likelihood
						// array that pertains to this pattern
						//
						cl = cl_start;
						double f = underflow->CalcCorrectionFactor(num_uf_corrections);
						for (z = 0; z < last; z++)
							*cl++ *= f;

						// Keep track of the number of times this pattern has been bounced so that 
						// the correction factors can be removed after the log of the site likelihood
						// has been computed
						//
						underflow->BounceSite(pattern, num_uf_corrections);
						}
				}	// loop over patterns
			}	// else child is internal node
		}	// loop over children
	}
#endif	// #if defined(POL_UNDERFLOW_CORRECTION)

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes `prMatrices[r]' for the TreeNodeAttribute object owned by node `nd' for all r = 0, 1, ..., ncat - 1.
|	Note that because data partitioning is not implemented, only the first ncat elements of prMatrices are used.
*/
void HKYAdHocEvaluator::RecalcPrMatrix(
  TreeNode *nd)		/**< is the node owning the TreeNodeAttribute for which the `prMatrices' array is to be recomputed. */
	{
	TreeNodeAttribute *attrib = (nd->attr)[0];
	attrib->pMatIsDirty = false;

	double ***p = attrib->prMatrices;
	double edgelength = nd->edgelen.f;
	for (unsigned r = 0; r < ncat; ++r)
		{
		RecalcSpecificPrMatrix(*p++, edgelength, relrate[r]);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes a specific transition probability matrix `p' using supplied relative rate `relrate'. Should only be
|	called from the RecalcPrMatrix(TreeNode *) member function.
*/
void HKYAdHocEvaluator::RecalcSpecificPrMatrix(
  double **p,			/**< is the 4x4 transition probability matrix to recalculate */
  double edgelen,		/**< is the edge length to use (not already multiplied by relative rate) */
  double relrate)		/**< is the relative rate to use for recalculating `p' */
	{
	double t = (relrate*edgelen);

	double Pi[4];
	Pi[A] = pi[A] + pi[G];
	Pi[C] = pi[C] + pi[T];
	Pi[G] = pi[A] + pi[G];
	Pi[T] = pi[C] + pi[T];

	double bigPiInv[4];
	bigPiInv[A] = 1.0 / Pi[A];
	bigPiInv[C] = 1.0 / Pi[C];
	bigPiInv[G] = 1.0 / Pi[G];
	bigPiInv[T] = 1.0 / Pi[T];

	double ta, tb, tc, td, y;
	double beta = 0.5 / ((pi[A] + pi[G])*(pi[C] + pi[T]) + kappa*((pi[A]*pi[G]) + (pi[C]*pi[T])));
	double x = exp(-beta*t);

	// changes to base A
	td = -beta*(1 + Pi[A]*(kappa - 1.0));
	y = exp(t*td);
	ta = pi[A]*(bigPiInv[A] - 1.0);
	tb = (Pi[A] - pi[A])*bigPiInv[A];
	tc = pi[A]*bigPiInv[A];
	p[A][A] = pi[A] + (x*ta) + (y*tb);
	p[C][A] = pi[A]*(1.0 - x);
	p[G][A] = pi[A] + (x*ta) - (y*tc);
	p[T][A] = p[C][A];

	// changes to base C
	td = -beta*(1 + Pi[C]*(kappa - 1.0));
	y = exp(t*td);
	ta = pi[C]*(bigPiInv[C] - 1.0);
	tb = (Pi[C] - pi[C])*bigPiInv[C];
	tc = pi[C]*bigPiInv[C];
	p[A][C] = pi[C]*(1.0 - x);
	p[C][C] = pi[C] + (x*ta) + (y*tb);
	p[G][C] = p[A][C];
	p[T][C] = pi[C] + (x*ta) - (y*tc);

	// changes to base G
	td = -beta*(1 + Pi[G]*(kappa - 1.0));
	y = exp(t*td);
	ta = pi[G]*(bigPiInv[G] - 1.0);
	tb = (Pi[G] - pi[G])*bigPiInv[G];
	tc = pi[G]*bigPiInv[G];
	p[A][G] = pi[G] + (x*ta) - (y*tc);
	p[C][G] = pi[G]*(1.0 - x);
	p[G][G] = pi[G] + (x*ta) + (y*tb);
	p[T][G] = p[C][G];

	// changes to base T
	td = -beta*(1 + Pi[T]*(kappa - 1.0));
	y = exp(t*td);
	ta = pi[T]*(bigPiInv[T] - 1.0);
	tb = (Pi[T] - pi[T])*bigPiInv[T];
	tc = pi[T]*bigPiInv[T];
	p[A][T] = pi[T]*(1.0 - x);
	p[C][T] = pi[T] + (x*ta) - (y*tc);
	p[G][T] = p[A][T];
	p[T][T] = pi[T] + (x*ta) + (y*tb);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates TreeNodeAttribute objects for each node in `tree' in preparation for computing likelihoods. Assumes `tree' 
|	is non-NULL. Note that because data partitioning is not implemented, only the 0th. element of the `attr' array is 
|	used. \todo Deletes and reallocates pre-existing attribute objects each time it is called, which should be 
|	avoided. 
*/
void HKYAdHocEvaluator::PrepareTree()
	{
	assert(tree != NULL);
	assert(npatterns > 0);

	// Create attributes structures for each node in the tree
	//
	for (TreeNode *nd = tree->GetFirstPreorder(); nd != NULL; nd = nd->GetNextPreorder())
		{
		//@POL shouldn't be deleting and recreating attributes each time likelihood is calculated
		//@POL note however that some functions (e.g. SetNumRateCategories) depend on this
		if (nd->attr.size() > 0)
			{
			assert(nd->attr.size() == 1);
			delete (nd->attr)[0];
			}

		(nd->attr).clear();
		TreeNodeAttribute *attrib = new TreeNodeAttribute(npatterns, ncat, NULL);
		(nd->attr).push_back(attrib);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Saves the patterns to the stream `o' along with information about the number of sites in the original data
|	matrix showing each pattern, and the status of the pattern (constant vs. variable, and for variable patterns, 
|	parsimony uninformative vs. informative).
*/
void HKYAdHocEvaluator::DebugShowPatterns(
  std::ostream &o,			/**< is the output stream on which to print the pattern summary */
  bool showOrigIndices)		/**< causes the `count' original character indices for each pattern to be displayed (can lead to very long lines) */
	{
	assert(pattern);

	unsigned i, j, z;

	// pos[i] is the current position along the states vector for taxon i
	//
	std::vector<unsigned> pos;
	pos.reserve(ntaxa);
	pos.assign(ntaxa, 0);

	string line = "   pattern \tfrequency \tpattern\t";
	o << line << std::endl;
	for (j = 0; j < npatterns; ++j)
		{
		unsigned pattern_freq = patfrq[(std::ptrdiff_t) j];
		std::string state_str;

		for (i = 0; i < ntaxa; ++i)
			{
			StateArray states = (StateArray) pattern[(std::ptrdiff_t) i];
			char s = states[pos[i]++];
			if (s == -4)
				state_str << '?';
			else if (s < 0)
				{
				unsigned ns = (unsigned)(-s);
				state_str = '{';
				for (z = 0; z < ns; ++z)
					{
					s = states[pos[i]++];
					if (z > 0)
						state_str << ',';
					switch (s)
						{
						case 0:		state_str << 'A'; break;
						case 1:		state_str << 'C'; break;
						case 2:		state_str << 'G'; break;
						default:	state_str << 'T';
						}
					}
				state_str << '}';
				}
			else
				{
				switch (s)
					{
					case 0:		state_str << 'A'; break;
					case 1:		state_str << 'C'; break;
					case 2:		state_str << 'G'; break;
					default:	state_str << 'T';
					}
				}
			}	// loop i over taxa

		line.clear();
		StrPrintF(line, "%10d \t%10d \t%s\t", j, pattern_freq, state_str.c_str());
		if (showOrigIndices && matrixShPtr)
			{
			const DataPatternInfo &info = matrixShPtr->GetPatternInfo(j);
			unsigned num_orig_indices = info.GetCount();
			StrPrintF(line, "%d", info.GetOrigIndex(0));
			for (unsigned z = 1; z < num_orig_indices; ++z)
				StrPrintF(line, ",%d", info.GetOrigIndex(z));
			}
		o << line << std::endl;
		}	// j loop over patterns
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Saves data in the form of a NEXUS TAXA block followed by a CHARACTERS block. Does not introduce the initial
|	#NEXUS keyword. Does not open the supplied stream `o', nor does it close `o'. Assumes `pattern' is not empty.
*/
void HKYAdHocEvaluator::DebugSaveData(std::ostream &o)
	{
	assert(pattern);

	unsigned i, j, k, z;
	unsigned longest_taxon_label = 0;

	o << "begin taxa;" << std::endl;
	o << "  dimensions ntax=" << ntaxa << ';' << std::endl;
	o << "  taxlabels" << std::endl;
	for (i = 0; i < ntaxa; ++i)
		{
		unsigned this_label_len = (unsigned)taxonLabels[i].length();
		if (this_label_len > longest_taxon_label)
			longest_taxon_label = this_label_len;
		o << "    '" << taxonLabels[i] << '\'' << std::endl;
		}
	o << "  ;" << std::endl;
	o << "end;" << std::endl;

	o << std::endl;

	o << "begin characters;" << std::endl;
	o << "  dimensions nchar=" << nchar << ';' << std::endl;
	o << "  format datatype=dna missing=? gap=-;" << std::endl;
	o << "  matrix" << std::endl;

	string format_str = MakeStrPrintF("%%%ds  ", longest_taxon_label + 4);

	for (i = 0; i < ntaxa; ++i)
		{
		const string quotedTaxonLabel = GetAsNexusToken(taxonLabels[i]);
		o << MakeStrPrintF(format_str.c_str(), quotedTaxonLabel.c_str());
		StateArray states = (StateArray) pattern[(std::ptrdiff_t) i];
		k = 0;	// k is current position in states vector
		for (j = 0; j < npatterns; ++j)
			{
			string state_str;
			char s = states[k++];
			if (s == -4)
				{
				state_str = '?';
				}
			else if (s < 0)
				{
				unsigned ns = (unsigned)(-s);
				state_str = '{';
				for (z = 0; z < ns; ++z)
					{
					s = states[k++];
					if (z > 0)
						state_str << ',';
					switch (s)
						{
						case 0:		state_str << 'A'; break;
						case 1:		state_str << 'C'; break;
						case 2:		state_str << 'G'; break;
						default:	state_str << 'T';
						}
					}
				state_str << '}';
				}
			else
				{
				switch (s)
					{
					case 0:		state_str = 'A'; break;
					case 1:		state_str = 'C'; break;
					case 2:		state_str = 'G'; break;
					default:	state_str = 'T';
					}
				}

			unsigned pattern_freq = patfrq[(std::ptrdiff_t) j];
			for (z = 0; z < pattern_freq; ++z)
				{
				o << state_str;
				}
			}	// j loop over states vector

		o << std::endl;
		}	// i loop over taxa

	o << "  ;" << std::endl;
	o << "end;" << std::endl;
	}

#if 0 //defined(POL_TEMP)

#	if defined(POL_UNDERFLOW_CORRECTION)
#		define OPEN_SITELIKE_LOGFILE	std::ofstream doof("sitelikes_uf.txt");
#	else
#		define OPEN_SITELIKE_LOGFILE	std::ofstream doof("sitelikes_nouf.txt");
#	endif

#	define OUTPUT_SITELIKE(site,siteLike,logSiteLike,frq)	\
		doof << (site) << '\t';								\
		if ((siteLike) == 0.0)								\
			doof << "***";									\
		else												\
			doof << MakeStrPrintF("%.12f", (logSiteLike));	\
		doof << '\t' << (frq) << std::endl;

#	define OUTPUT_COMMENT(s)								\
		doof << (s) << std::endl;
		
#	define CLOSE_SITELIKE_LOGFILE_AND_EXIT					\
		/* DebugShowPatterns(doof); */						\
		doof.close();										\
		exit(0);

#else

#	define OPEN_SITELIKE_LOGFILE
#	define OUTPUT_SITELIKE(site,siteLike,logSiteLike,frq)
#	define OUTPUT_COMMENT(s)
#	define CLOSE_SITELIKE_LOGFILE_AND_EXIT

#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the log-likelihood of `tree' under the HKY85 model unless data member `ignoreData' is true, in which case
|	this function simply returns 0.0. Returns -DBL_MAX if `paramOutOfBounds' is true. This variable will be set to true
|	when the slice sampler tries values of a parameter that are outside its domain. Calls CalcLogLikelihoodImpl to do 
|	the actual work.
*/
double HKYAdHocEvaluator::CalcLogLikelihood()
	{
	assert(tree != NULL);

	if (ignoreData)
		return 0.0;

	if (paramOutOfBounds)
		return -DBL_MAX;

	double lnL;
	try
		{
		// Try using no underflow protection first
		//
		lnL = CalcLogLikelihoodImpl();
		}
	catch(HKYAdHocEvaluator::XSiteLikeUnderflow)
		{
		// Underflow was a problem first time, so need to recalc log-likelihood.
		// Protection will be in place this time for sites that needed it.
		//
		lnL = CalcLogLikelihoodImpl();
		}

	return lnL;
	}

#if defined(POL_UNDERFLOW_CORRECTION)
#	if defined(POL_TEMP)
	std::string HKYAdHocEvaluator::DebugShowBouncesVector()
		{
		return MakeStrPrintF("%.1f", underflow->DebugGetAvgNumBounces());
		//return underflow->DebugCreateBounceVectorRepresentation();
		}
#	endif
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Does the actual work of computing the log-likelihood of `tree' under the HKY85 model. This function should not be
|	called directly; rather, CalcLogLikelihood() should be called and it will call this function. This function will
|	throw an XSiteLikeUnderflow exception if any site likelihood is zero. This is a signal to the calling function 
|	that the log-likelihood is not valid and must be recalculated, this time with underflow protection enabled for
|	those sites that underflowed the first time. This implementation does not work out from a focal edge, so it will 
|	be inefficient if the dirty parts of the tree are a long way from the root.
*/
double HKYAdHocEvaluator::CalcLogLikelihoodImpl()
	{
	bool underflow_happened = false;

#if defined(POL_UNDERFLOW_CORRECTION)
	// Bounces vector initially is UINT_MAX for each site. Sites that underflow get reset to 0
	// for next attempt. Underflow protection is only used for sites that have demonstrated a need
	// for it. (Note: "sites" here really means "site patterns".)
	//
	underflow->InitBouncesVector(npatterns);
	bool UF_protect = false;
	if (underflow->AnySiteHasUnderflowed())
		{
		UF_protect = true;
		}
#endif

	OPEN_SITELIKE_LOGFILE

	// Perform a postorder traversal, making sure all conditional likelihood arrays are up to date.
	//
	for (TreeNode *nd = tree->GetFirstPostorder(); nd != NULL; nd = nd->GetNextPostorder())
		{
		TreeNodeAttribute *attrib = GetNodeAttr(nd);

		if (!nd->IsRoot() && (entireTreeDirty || attrib->pMatIsDirty))
			{
			RecalcPrMatrix(nd);

			// Make sure conditional likelihood array is recomputed when 
			// the parent of this node is visited
			//
			assert(nd->GetParent() != NULL);
			TreeNodeAttribute *par_attr = GetNodeAttr(nd->GetParent());
			par_attr->clDirty = true;
			}

		if (!nd->IsShootTip() && (entireTreeDirty || attrib->clDirty))
			{
#if defined(POL_UNDERFLOW_CORRECTION)
			if (UF_protect)
				UFProtectedCondLike(nd);
			else
				RecalcCondLike(nd);
#else
			RecalcCondLike(nd);
#endif

			// This next part would not be necessary if we calculated likelihood out from
			// a focal node rather than always collecting down to root
			//
			if (nd->GetParent() != NULL)
				{
				TreeNodeAttribute *par_attr = GetNodeAttr(nd->GetParent());
				par_attr->clDirty = true;
				}
			}
		} // end of postorder traversal of tree

	TreeNode *rootNode = tree->GetRoot();
	unsigned rootNodeNum = rootNode->GetNodeNumber();
	assert(rootNodeNum == 0);

	TreeNodeAttribute *root_attr = GetNodeAttr(rootNode);
	double lnL = 0.0;
	double *cl = root_attr->condLike;

	// Visit all states observed in the root node (which is actually a leaf) and, for each, construct the site likelihood 
	// as the relative frequency of that state times the conditional likelihood of the entire tree above the root for that state.
	// If there is ambiguity, must sum over all possible states. If there is rate heterogeneity (e.g. ncat > 1), then 
	// assuming root node has state A, the site likelihood L_k is
	//
	//   L_k = pi_A Pr(data | A)
	//       = pi_A [Pr(data | A, r_1] Pr(r_1) + Pr(data | A, r_2] Pr(r_2) + Pr(data | A, r_3] Pr(r_3)]
	//       = pi_A [Pr(data | A, r_1] + Pr(data | A, r_2] + Pr(data | A, r_3]] / 3
	//
	// Here, 'data' refers to observed states at all other leaves besides the root leaf, and all three rate categories are 
	// equiprobable.
	//
	// Conditional likelihood arrays at interior nodes are laid out as follows:
	//
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	// | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T |
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	// |    rate 1     |    rate 2     |    rate 3     |    rate 4     |    rate 1     |    rate 2     |    rate 3     |    rate 4     |
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	// |                          pattern 1                            |                         pattern 2                             |
	// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
	//
	StateArray states = (StateArray) pattern[(std::ptrdiff_t) rootNodeNum];
	unsigned len = (unsigned)states.size();
	unsigned i = 0;	// index into states array
	unsigned k = 0;	// pattern index
	char *s = &states[0];

#if defined(POL_UNDERFLOW_CORRECTION)
	unsigned *x = NULL;
	if (UF_protect)
		x = underflow->GetBouncesArray();
#endif

	unsigned r;		// for looping through rate categories
	for (; i < len; ++s, ++i)
		{
		double siteLike = 0.0;
		char state = *s;
		if (state == -4)
			{
			double sum[] = {0.0, 0.0, 0.0, 0.0};
			for (r = 0; r < ncat; ++r)
				{
				sum[A] += cl[A];
				sum[C] += cl[C];
				sum[G] += cl[G];
				sum[T] += cl[T];
				cl += 4;
				}
			siteLike  = pi[A]*sum[A];
			siteLike += pi[C]*sum[C];
			siteLike += pi[G]*sum[G];
			siteLike += pi[T]*sum[T];
			siteLike /= (double)ncat;
			}
		else if (state < 0)
			{
			unsigned ns = (unsigned)(-state);
			double sum[] = {0.0, 0.0, 0.0, 0.0};
			assert(i + ns < len);
			char *s0 = s;
			unsigned i0 = i;
			for (r = 0; r < ncat; ++r)
				{
				s = s0;
				i = i0;
				for (unsigned z = 0; z < ns; ++z)
					{
					++s;
					++i;
					state = *s;
					assert(state == A || state == C || state == G || state == T);
					sum[(int) state] += cl[(int)state];
					}
				cl += 4;
				}
			// Note: some of the elements of sum array will still be 0.0, but those will not
			// contribute to siteLike so it doesn't hurt to sum over all bases
			//
			siteLike  = pi[A]*sum[A];
			siteLike += pi[C]*sum[C];
			siteLike += pi[G]*sum[G];
			siteLike += pi[T]*sum[T];
			siteLike /= (double)ncat;
			}
		else
			{
			assert(state == A || state == C || state == G || state == T);
			double sum = 0.0;
			for (r = 0; r < ncat; ++r)
				{
				sum += cl[state];	//@#POL should precalculate cl[state]
				cl += 4;
				}
			siteLike = pi[(int) state]*sum/(double)ncat;	//@#POL should multiply by precalculated inverse of (double)ncat
			}

		double logSiteLike = 0.0;
#if defined(POL_UNDERFLOW_CORRECTION)
		if (siteLike == 0.0)
			{
			underflow_happened = true;
			underflow->ProtectSite(k);
			}
		else
			{
			if (UF_protect)
				{
				if (*x++ < UINT_MAX)
					logSiteLike = underflow->GetCorrectedLnLike(siteLike, k);
				else
					logSiteLike = std::log(siteLike);
				}
			else
				logSiteLike = std::log(siteLike);
			}
#else
		logSiteLike = (siteLike > 0.0 ? std::log(siteLike) : 0.0);
#endif
		
		OUTPUT_SITELIKE(k,siteLike,logSiteLike,patfrq[k])

		double freq = (double)patfrq[(std::ptrdiff_t)k]; //@#POL should walk down patfrq rather than seek index each time
		logSiteLike *= freq;
		lnL += logSiteLike;

		k++;
		}	// loop across patterns

	if (underflow_happened)
		{
		entireTreeDirty = true;
		throw HKYAdHocEvaluator::XSiteLikeUnderflow();
		}
	else
		{
		entireTreeDirty = false;
		}

	CLOSE_SITELIKE_LOGFILE_AND_EXIT

	return lnL;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `ncat' and causes the base class member function CalculateRates() to be called to refresh the
|	`relrate' data member.
*/
void HKYAdHocEvaluator::SetNumRateCategories(
  unsigned n)	/**< is the number of rate categories */
	{
	ncat = n;

	// Tell the base class PhoRateManager about the change in rate categories
	//
	SetNumVariableRateCategories(ncat);

	// Refresh the relrate array
	//
	relrate.clear();
	relrate = GetRepRates();
	if (tree != NULL)
		PrepareTree();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If the `attr' vector for `nd' is empty, creates a new one and returns a pointer to it; otherwise, returns pointer to
|	the existing TreeNodeAttributes object at nd->attr[0]. Note that HKYAdHocEvaluator only uses the first element in
|	the `attr' vector because it assumes data are not partitioned.
*/
TreeNodeAttribute *HKYAdHocEvaluator::GetNodeAttr(TreeNode *nd)
	{
	if (nd->attr.empty())
		nd->attr.push_back(new TreeNodeAttribute(npatterns, 1, NULL));
	return nd->attr[0];
	}

