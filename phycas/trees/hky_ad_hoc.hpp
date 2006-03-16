#ifndef HKY_ADHOC_H
#define HKY_ADHOC_H

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include "phycas/likelihood/rate_manager.hpp"
#include "ncl/output/temp_basic_output_operators.hpp"

class CompressedMatrix;
class GibbsParameter;
class PhoTaxaManager;
class PhoCharactersManager;
class TreeNodeAttribute;
class Tree;
class UnderflowManager;
typedef boost::shared_ptr<const CompressedMatrix>	CompressedMatrixShPtr;
typedef boost::shared_ptr<GibbsParameter>			GibbsParameterShPtr;
typedef std::vector<GibbsParameterShPtr>			VecGibbsParameterShPtr;
typedef boost::shared_array<unsigned>				UintShArr;
typedef std::vector<char>							StateArray; 
typedef boost::shared_array<StateArray>				DataMatrix;
typedef boost::shared_ptr<UnderflowManager>			UnderflowShPtr;

class TreeNode;

/*----------------------------------------------------------------------------------------------------------------------
|	This class is a temporary replacement for Mark's likelihood computing machinery, which is currently out for repairs 
|	and upgrades. The HKYAdHocEvaluator, as the name suggests, can only compute likelihoods for the HKY85 model. Call 
|	the SetTree member function to provide a tree, the CopyData member function to copy data from a PhoCharactersManager
|	object, PrepareTree to construct all the node attribute structures containing conditional likelihood arrays, and 
|	CalcLogLikelihood to compute the log-likelihood of the tree. 
|	
|	Current limitations include:
|	- no underflow correction
|	- some output uses std::cerr and std::cout
*/
class HKYAdHocEvaluator : public PhoRateManager
	{
	public:

		class XSiteLikeUnderflow {};	/**< Thrown if underflow occurs when calculating site likelihoods in CalcLogLikelihoodImpl */

		/** Enumeration storing the Index of each model parameter */
		enum ParamIndex
			{
			kHKYFreqAIndex = 0,		/**< Index of frequency of A in parameter vector */
			kHKYFreqCIndex,			/**< Index of frequency of C in parameter vector */
			kHKYFreqGIndex,			/**< Index of frequency of G in parameter vector */
			kHKYFreqTIndex,			/**< Index of frequency of T in parameter vector */
			kHKYKappaIndex,			/**< Index of the transition/transversion rate ratio kappa in parameter vector */
			kHKYRateVarIndex,		/**< Index of the discrete gamma rate heterogeneity variance parameter in parameter vector */
			kHKYPInvarIndex,		/**< Index of the pinvar parameter in parameter vector */
			kHKYMeanRateIndex,		/**< Index of the mean rate parameter in parameter vector */
			kHKYLastParam			/**< Numerical value equals number of defined parameters */
			};

							HKYAdHocEvaluator(
								GibbsParameterShPtr freqAParam, 
								GibbsParameterShPtr freqCParam, 
								GibbsParameterShPtr freqGParam, 
								GibbsParameterShPtr freqTParam, 
								GibbsParameterShPtr kappaParam, 
								GibbsParameterShPtr rateVarParam, 
								GibbsParameterShPtr pInvarParam, 
								GibbsParameterShPtr meanRatesParam, 
								unsigned ncateg);
							~HKYAdHocEvaluator();

		void				CopyData(PhoTaxaManager & tm, PhoCharactersManager & cm);
		void				Clear();
		void				PrepareTree();
		void				InvalidateTree();
		TreeNodeAttribute * GetNodeAttr(TreeNode * nd);

		double				GetGammaShape();
		double				GetRateVariance();

		double				GetKappa();
		double				CalcTRatio();

		void				ParameterChanged();

		void				DebugSaveData(std::ostream &outf);
		void				DebugShowPatterns(std::ostream &outf, bool showOrigIndices = false);
		void				DebugShowCondLikeArrays(TreeNode *nd);

#if defined(POL_UNDERFLOW_CORRECTION)
		void				UFProtectedCondLike(TreeNode *nd);
#	if defined(POL_TEMP)
		std::string			DebugShowBouncesVector();
#	endif
#endif
		void				RecalcCondLike(TreeNode *nd);
		void				RecalcPrMatrix(TreeNode *nd);

		double				CalcLogLikelihood();
		std::vector<double> GetNucleotideFrequencies();

		void				SetTree(Tree * t);
		unsigned			GetNTaxa();
		unsigned			GetNPatterns();
		unsigned			GetNChar();

		void				SetIgnoreData(bool ignore);
		void				SetNumRateCategories(unsigned n);

		unsigned			GetNumUpdatableParameters();
		void				ShowUpdatableParameters(NxsOutputStream &out) const;
		void				ResetParameterPtr(ParamIndex pCode, GibbsParameterShPtr p);

		const std::vector<std::string> & GetTaxonLabels() const
			{
			return taxonLabels;
			}
	private:

		double				CalcLogLikelihoodImpl();
		void				RecalcSpecificPrMatrix(double **p, double edgelen, double relrate);

	private:

		std::vector<std::string>	taxonLabels;		/**< Vector of taxon label strings */
		std::vector<std::string>	patternStr;			/**< Vector of patterns represented as strings */

		CompressedMatrixShPtr		matrixShPtr;		/**< pointer to compressed matrix, set in CopyData and used in DebugShowPatterns */

		// Currently supporting two ways of storing data for tip nodes
		//
		DataMatrix				pattern;				/**< Pattern matrix: row i is vector of ord. coded states for taxon i */

		bool					entireTreeDirty;		/**< True if InvalidateTree called. Overrides pMatDirty and clDirty flags in TreeAttributes structures. Reset to false in CalcLogLikelihood(). */
		bool					paramOutOfBounds;		/**< True if slice sampler tries a value that is out of bounds for any parameter; causes likelihood to equal 0.0 */

		UintShArr				patfrq;					/**< Pattern frequencies */

		bool					ignoreData;				/**< If true, CalcLogLikelihood() always returns 0.0. */

		Tree			*		tree;					/**< Tree on which likelihoods are needed */
		unsigned				ntaxa;					/**< No. taxa */
		unsigned				nchar;					/**< No. original sites */
		unsigned				npatterns;				/**< No. sites after crunching */

		double					pi[4];					/**< Empirical base frequencies */
		double					kappa;					/**< Value of the transition/transversion rate ratio parameter */
		double					rateVar;				/**< Value of the discrete gamma rate heterogeneity variance parameter */
		unsigned				ncat;					/**< Discrete gamma rate heterogeneity number of categories */

#if defined(POL_UNDERFLOW_CORRECTION)
		UnderflowShPtr			underflow;				/**< Object that ensures site likelihoods do not underflow */
		//UnderflowManager		*underflow;
#endif

		std::vector<double>		relrate;				/**< The representative relative rates for each rate category */

		VecGibbsParameterShPtr  modelParams;			/**< Array of pointers to the model's parameters */

		std::ofstream			outf;
	};
	
template<class OUT_STREAM>
class GenericPrinterClass<kVerboseModelParameterOutStyle, HKYAdHocEvaluator, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const HKYAdHocEvaluator & hkyAHE)
			{
			outStream << "\n\nParameter summary:\n";
			hkyAHE.ShowUpdatableParameters(outStream);
			}
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of `ignoreData' data member. If `ignoreData' is set to true, the CalcLogLikelihood() function always 
|	returns 0.0; otherwise, it actually calculates the log-likelihood.
*/
inline void HKYAdHocEvaluator::SetIgnoreData(bool ignore)
	{
	ignoreData = ignore;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reassigns model parameter at index `pCode' to `p', where `pCode' is one of the constants in the ParamIndex enum.
*/
inline void HKYAdHocEvaluator::ResetParameterPtr(
  ParamIndex pCode,			/**< is the postition of the parameter in the model's param array */
  GibbsParameterShPtr p) 	/**< is a shared pointer to the new parameter */
	{
	assert(pCode < HKYAdHocEvaluator::kHKYLastParam);
	modelParams[pCode] = p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor calls Clear() to do all the work.
*/
inline HKYAdHocEvaluator::~HKYAdHocEvaluator() 
	{
	Clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows access to the value of `rateVar'.
*/
inline double HKYAdHocEvaluator::GetRateVariance()
	{
	return rateVar;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the inverse of `rateVar', which equals the shape parameter of the gamma distribution of relative rates 
|	across sites.
*/
inline double HKYAdHocEvaluator::GetGammaShape()
	{
	assert(rateVar > 0.0);
	return (1.0/rateVar);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows access to the value of `kappa'.
*/
inline double HKYAdHocEvaluator::GetKappa()
	{
	return kappa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calculates the transition/transversion ratio given the transition/transversion rate ratio `kappa' and the relative
|	base frequencies, which are stored in the `pi' array. Here are the details of the calculation:
|>
|	Parameters: b = transversion rate, k = kappa
|	Pr(any transition | dt)   = Pr(AG) + Pr(CT) + Pr(GA) + Pr(TC) 
|	                          = (piA piG k b dt) + (piC piT k b dt) + (piG piA k b dt) + (piT piC k b dt)
|	                          = 2 k b dt (piA piG + piC piT)
|	
|	Pr(any transversion | dt) = Pr(AC) + Pr(AT) + Pr(CA) + Pr(CG) + Pr(GC) + Pr(GT) + Pr(TA) + Pr(TG)
|	                          = (piA piC b dt) + (piA piT b dt) + (piC piA b dt) + (piC piG b dt)
|	                            + (piG piC b dt) + (piG piT b dt) + (piT piA b dt) + (piT piG b dt)
|	                          = 2 b dt (piA + piG) (piC + piT)
|	
|	          2 k b dt (piA piG + piC piT)     k (piA piG + piC piT)
|	TRatio = ------------------------------ = -----------------------
|	         2 b dt (piA + piG) (piC + piT)   (piA + piG) (piC + piT)
|<
*/
inline double HKYAdHocEvaluator::CalcTRatio()
	{
	double numerator   = kappa * ((pi[0] * pi[2]) + (pi[1] * pi[3]));
	double denominator = (pi[0] + pi[2]) * (pi[1] + pi[3]);
	double tratio = numerator / denominator;
	return tratio;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a new vector from the relative base frequencies stored in the `pi' array and returns this vector.
*/
inline std::vector<double> HKYAdHocEvaluator::GetNucleotideFrequencies()
	{
	return std::vector<double> (pi, (pi + 4));
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the data member `tree' to point to the supplied Tree pointer `t'.
*/
inline void HKYAdHocEvaluator::SetTree(
  Tree *t)	/**< is a pointer to the new Tree */
	{
	tree = t;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows access to the value of `ntaxa'.
*/
inline unsigned HKYAdHocEvaluator::GetNTaxa()
	{
	return ntaxa;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows access to the value of `npatterns'.
*/
inline unsigned HKYAdHocEvaluator::GetNPatterns()
	{
	return npatterns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allows access to the value of `nchar*'.
*/
inline unsigned HKYAdHocEvaluator::GetNChar()
	{
	return nchar;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Invalidate transition probability matrices and conditional likelihood arrays for all nodes in `tree'.
*/
inline void HKYAdHocEvaluator::InvalidateTree()
	{
	entireTreeDirty = true;
	}

#endif
