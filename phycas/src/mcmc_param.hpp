/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
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

#if ! defined(MCMC_PARAM_HPP)
#define MCMC_PARAM_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "phycas/src/mcmc_updater.hpp"				// for MCMCUpdater base class

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

class Codon;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the nonsynonymous/synonymous rate ratio parameter of a codon model. This parameter is commonly 
|	symbolized using the Greek letter omega, hence the name.
*/
class OmegaParam : public MCMCUpdater
	{
	public:
						OmegaParam();
						virtual ~OmegaParam() 
							{
							//std::cerr << "OmegaParam dying..." << std::endl;
							}

		virtual void	setModel(ModelShPtr p);
		virtual void	update();				// override virtual from MCMCUpdater base class
		virtual double	operator()(double k);	// override pure virtual from AdHocDensity (base class of MCMCUpdater)

	private:

		Codon *			codon;
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

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates an among-sites relative rate parameter (component of the FLEXCAT model). More precisely, this class 
|	encapsulates a parameter that governs an unnormalized FLEXCAT relative rate. Before being used in the likelihood
|	function, all FlexRateParam values are normalized to have mean 1.0.
*/
class FlexRateParam : public MCMCUpdater
{
	public:
									FlexRateParam(unsigned & s, double & ub, std::vector<double> & rr);
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

		unsigned &					nspacers;
		double &					upper_bound;
		std::vector<double>	&		rel_rates;
		double						left_value;
		double						right_value;
		unsigned					which;
};

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates an among-sites relative rate category probability parameter (component of the FLEXCAT model). More 
|	precisely, this class encapsulates a parameter that governs an unnormalized FLEXCAT rate category probability. 
|	Before being used in the likelihood function, all FlexProbParam values are normalized to have sum 1.0.
*/
class FlexProbParam : public MCMCUpdater
{
	public:
								FlexProbParam(std::vector<double> & rp);
								virtual ~FlexProbParam()
									{
									//std::cerr << "FlexProbParam dying..." << std::endl;
									}
	
		virtual void			update();				// override virtual from MCMCUpdater base class
		virtual void			setLot(LotShPtr r);
		virtual double			operator()(double f);

	private:

		unsigned				which;
		std::vector<double>	&	rate_probs;
};

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates a state frequency parameter. More precisely, this class encapsulates a parameter that governs a state
|	frequency. Because all state frequencies must add to 1.0, the value of this parameter and the other parameter values
|	must be normalized so that the model has a correct set of state frequencies.
*/
class StateFreqParam : public MCMCUpdater
{
	public:
							StateFreqParam(unsigned w);
							virtual ~StateFreqParam()
								{
								//std::cerr << "StateFreqParam dying..." << std::endl;
								}
	
		virtual void		update();				// override virtual from MCMCUpdater base class
		virtual double		operator()(double k);

	private:

		unsigned			which;
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
        enum                EdgeLenType
                                {
                                internal,   /**> only computes prior for internal edges */
                                external,   /**> only computes prior for external edges */
                                both        /**> computes prior for all edges */
                                };

                            EdgeLenMasterParam(EdgeLenMasterParam::EdgeLenType t = EdgeLenMasterParam::both);
							virtual ~EdgeLenMasterParam()
								{
								//std::cerr << "EdgeLenMasterParam dying..." << std::endl;
								}
		virtual double		recalcPrior();
		virtual void		setPriorMeanAndVariance(double m, double v);

	protected:

		double				lnPriorOneEdge(TreeNode & nd) const;

    private:

        EdgeLenType         edgeLenType;    /**> holds the edge length type, which determines for which edge lengths the prior is computed when recalcPrior is called */
	};

typedef boost::shared_ptr<EdgeLenMasterParam> EdgeLenMasterParamShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Represents the mean of the edge length prior. This parameter is used in models employing a hyperprior governing the
|	mean of the edge length priors. Whenever the HyperPriorParam value changes, the means of the priors affecting all
|	EdgeLenParam objects in the MCMCChainManager are reset to this new value.
*/
class HyperPriorParam : public MCMCUpdater
	{
	public:
							HyperPriorParam();
							HyperPriorParam(EdgeLenMasterParamShPtr p);
							virtual ~HyperPriorParam()
								{
								//std::cerr << "HyperPriorParam dying..." << std::endl;
								}

		virtual void		update();				// override virtual from MCMCUpdater base class
		virtual double		operator()(double k);	

    private:

        EdgeLenMasterParamShPtr    edgelen_master_param;   /**> is the edge length master parameter whose prior this parameter controls */
	};

} // namespace phycas

#include "phycas/src/mcmc_param.inl"

#endif
