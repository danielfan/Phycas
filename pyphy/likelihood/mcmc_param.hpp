#if ! defined(MCMC_PARAM_HPP)
#define MCMC_PARAM_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "pyphy/likelihood/mcmc_updater.hpp"		// for MCMCUpdater base class
//#include "phycas/modules/mcmc/slice_sampler.hpp"	// for AdHocDensity, SliceSampler, SliceSamplerShPtr, LotShPtr, ProbDistShPtr

//struct CIPRES_Matrix;

//namespace CipresNative
//{
//class DiscreteMatrix;
//}
	
//class ProbabilityDistribution;

namespace phycas
{

class HKY;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the transition/transversion rate ratio parameter of models such as the HKY model. This parameter is
|	commonly symbolized using the Greek letter kappa, hence the name.
*/
class KappaParam : public MCMCUpdater
	{
	public:
						KappaParam();
						virtual ~KappaParam() 
							{
							//std::cerr << "KappaParam dying..." << std::endl;
							}

		virtual void	setModel(ModelShPtr p);
		virtual void	update();				// override virtual from MCMCUpdater base class
		virtual double	operator()(double k);	// override pure virtual from AdHocDensity (base class of MCMCUpdater)

	private:

		HKY *			hky;
	};

class GTR;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates one of the relative rate parameters in the GTR model. 
*/
class GTRRateParam : public MCMCUpdater
	{
	public:
						GTRRateParam(unsigned w);
						virtual ~GTRRateParam() 
							{
							//std::cerr << "GTRRateParam dying..." << std::endl;
							}

		virtual void	setModel(ModelShPtr p);
		virtual void	update();				// override virtual from MCMCUpdater base class
		virtual double	operator()(double k);	// override pure virtual from AdHocDensity (base class of MCMCUpdater)

	private:

		GTR *			gtr;		/**< */
		unsigned		which;		/**< */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the discrete gamma shape parameter.
*/
class DiscreteGammaShapeParam : public MCMCUpdater
	{
	public:
						DiscreteGammaShapeParam(bool invert);
						virtual ~DiscreteGammaShapeParam() 
							{
							//std::cerr << "DiscreteGammaShapeParam dying..." << std::endl;
							}

		virtual void	update();				// override virtual from MCMCUpdater base class
		virtual double	operator()(double k);	// override pure virtual from AdHocDensity (base class of MCMCUpdater)

	private:

		bool			invert_shape;	/**< If true, inverse of shape parameter (rather than the shape parameter itself) is managed by this updater */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the proportion of invariable sites parameter.
*/
class PinvarParam : public MCMCUpdater
	{
	public:
						PinvarParam();
						virtual ~PinvarParam()
							{
							//std::cerr << "PinvarParam dying..." << std::endl;
							}

		virtual void	update();				// override virtual from MCMCUpdater base class
		virtual double	operator()(double k);	// override pure virtual from AdHocDensity (base class of MCMCUpdater)
	};

#if POLPY_NEWWAY
/*----------------------------------------------------------------------------------------------------------------------
|	
*/
class FlexRateParam : public MCMCUpdater
{
	public:
									FlexRateParam(unsigned w, std::vector<double> & rr);
									virtual ~FlexRateParam()
										{
										//std::cerr << "FlexRateParam dying..." << std::endl;
										}
	
		virtual void				update();				// override virtual from MCMCUpdater base class
		virtual std::string			getPriorDescr() const;	// override virtual from MCMCUpdater base class
		virtual double				recalcPrior();			// override virtual from MCMCUpdater base class
		virtual void				setChainManager(ChainManagerWkPtr p);	// override virtual from MCMCUpdater base class
		virtual double				operator()(double r);

	private:

		void						refreshLeftRightValues();
	private:

		std::vector<double>	&		rel_rates;
		double						left_value;
		double						right_value;
		unsigned					which;
};

/*----------------------------------------------------------------------------------------------------------------------
|	
*/
class FlexProbParam : public MCMCUpdater
{
	public:
								FlexProbParam(unsigned w);
								virtual ~FlexProbParam()
									{
									//std::cerr << "FlexProbParam dying..." << std::endl;
									}
	
		virtual void			update();				// override virtual from MCMCUpdater base class
		virtual double			operator()(double f);

	private:

		unsigned				which;
};
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates a base frequency parameter. More precisely, this class encapsulates a parameter that governs a base
|	frequency. Because all four base frequencies must add to 1.0, the value of this parameter and the other three 
|	parameter values must be normalized so that the model has a correct set of base frequencies.
*/
class BaseFreqParam : public MCMCUpdater
{
	public:
							BaseFreqParam(unsigned w);
							virtual ~BaseFreqParam()
								{
								//std::cerr << "BaseFreqParam dying..." << std::endl;
								}
	
		virtual void		update();				// override virtual from MCMCUpdater base class
		virtual double		operator()(double k);

	private:

		unsigned			which;
};

/*----------------------------------------------------------------------------------------------------------------------
|	Represents the mean of the edge length prior. This parameter is used in models employing a hyperprior governing the
|	mean of the edge length priors. Whenever the HyperPriorParam value changes, the means of the priors affecting all
|	EdgeLenParam objects in the MCMCChainManager are reset to this new value.
*/
class HyperPriorParam : public MCMCUpdater
	{
	public:
							HyperPriorParam();
							virtual ~HyperPriorParam()
								{
								//std::cerr << "HyperPriorParam dying..." << std::endl;
								}

		virtual void		update();				// override virtual from MCMCUpdater base class
		virtual double		operator()(double k);	
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Represents all edge lengths in the associated Tree; i.e., and edge length "master" parameter. An edge length master 
|	parameter does not update any particular edge length (that is left up to one of the defined moves, such as a 
|	`LargetSimonMove); instead, its	job is to compute the joint prior over all edge lengths. The update() member 
|	function and operator() of an EdgeLenMasterParam are not overridden because they do nothing.
*/
class EdgeLenMasterParam : public MCMCUpdater
	{
	public:
							EdgeLenMasterParam();
							virtual ~EdgeLenMasterParam()
								{
								//std::cerr << "EdgeLenMasterParam dying..." << std::endl;
								}

		virtual double		recalcPrior();
		virtual void		setPriorMeanAndVariance(double m, double v);

	protected:

		double				lnPriorOneEdge(TreeNode & nd) const;
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Represents the length of the edge associated with the TreeNode `nd'.
*/
class EdgeLenParam : public EdgeLenMasterParam
	{
	public:

							EdgeLenParam(TreeNode * node);
							virtual ~EdgeLenParam()
								{
								//std::cerr << "EdgeLenParam dying..." << std::endl;
								}

		void				update();  // override virtual from MCMCUpdater base class
		virtual double		operator()(double k);

	private:

		TreeNode *			nd;		/**< The node whose edge is being manipulated by this parameter */
	};

} // namespace phycas

#include "pyphy/likelihood/mcmc_param.inl"

#endif
