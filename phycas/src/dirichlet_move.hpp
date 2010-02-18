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

#if ! defined(STATE_FREQ_MOVE_HPP)
#define STATE_FREQ_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
#include "phycas/src/mcmc_updater.hpp"		// for base class MCMCUpdater

namespace phycas
{

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

typedef boost::shared_ptr<DirichletDistribution>    DirichletShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	A DirichletMove proposes new parameter values that are slightly different than the current parameter values by 
|   sampling from a Dirichlet distribution with parameters equal to the current frequencies multiplied by a large 
|   value (the tuning parameter 'psi').
*/
class DirichletMove : public MCMCUpdater
	{
	public:
									DirichletMove();
									virtual ~DirichletMove() 
										{
										//std::cerr << "DirichletMove dying..." << std::endl;
										}

		// Accessors
		double						getPsi() const;
		unsigned					getDimension() const;


		// Modifiers
		void						setPsi(double x);
		void						setMaxPsi(double x);
		void						setMinPsi(double x);
		void						setDimension(unsigned d);

		// Utilities
		void						reset();
		virtual void				sendCurrValuesToModel(const double_vect_t & v) {}
		virtual void				getCurrValuesFromModel(double_vect_t & v) const {}
		virtual double_vect_t		listCurrValuesFromModel();
        virtual void                getParams() {}
        virtual void                setParams(const std::vector<double> & v) {}

		// These are virtual functions in the MCMCUpdater base class
        virtual void                setPosteriorTuningParam(double x);
        virtual void                setPriorTuningParam(double x);
		virtual void				setBoldness(double x);
		virtual bool				update();
		virtual double				getLnHastingsRatio() const;
		virtual double				getLnJacobian() const;
		virtual void				proposeNewState();
		virtual void				revert();
		virtual void				accept();
		
		virtual void				educateWorkingPrior();
		virtual void				finalizeWorkingPrior();

	private:

		DirichletMove &				operator=(const DirichletMove &);	// never use - don't define

	protected:

		unsigned					dim;			/**< The number of parameters governed by this move */
		double						boldness;		/**< Ranges from 0 to 100 and determines the boldness of the move */
		double						psi;			/**< Larger values result in changes of smaller magnitude */
		double						min_psi;		/**< Smallest allowed value of psi (used when exploring the prior) */
		double						max_psi;		/**< Largest allowed value of psi (used when exploring the posterior) */
		std::vector<double> 		new_params;	    /**< Proposed new parameter values */
		std::vector<double> 		orig_params;	/**< Saved parameter values (in case revert is necessary) */
		std::vector<double> 		c_forward;	    /**< Dirichlet parameter vector used to propose new frequencies */
		std::vector<double> 		c_reverse;	    /**< Dirichlet parameter vector used to propose original frequencies (only used to compute Hastings ratio) */
        DirichletShPtr              dir_forward;    /**< Points to an ad hoc Dirichlet distribution object used to assess the forward move density and to propose a new frequency vector */
        DirichletShPtr              dir_reverse;    /**< Points to an ad hoc Dirichlet distribution object used to assess the reverse move density */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	A StateFreqMove proposes new state frequencies that are slightly different than the current frequencies by sampling 
|   from a Dirichlet distribution with parameters equal to the current frequencies multiplied by a large value (the 
|   tuning parameter 'psi').
*/
class StateFreqMove : public DirichletMove
	{
	public:
									StateFreqMove();
                                    virtual ~StateFreqMove() {}

		virtual void				educateWorkingPrior();

		virtual void				sendCurrValuesToModel(const double_vect_t & v);
		virtual void				getCurrValuesFromModel(double_vect_t & v) const;
		virtual double_vect_t		listCurrValuesFromModel();
        virtual void                getParams();
        virtual void                setParams(const std::vector<double> & v);

	private:

		StateFreqMove &				operator=(const StateFreqMove &);	// never use - don't define
    };

/*----------------------------------------------------------------------------------------------------------------------
|	A RelRatesMove proposes new GTR relative rates that are slightly different than the current rates by sampling 
|   from a Dirichlet distribution with parameters equal to the current rates multiplied by a large value (the 
|   tuning parameter 'psi').
*/
class RelRatesMove : public DirichletMove
	{
	public:
									RelRatesMove();
                                    virtual ~RelRatesMove() {}

		virtual void				educateWorkingPrior();

		virtual void				sendCurrValuesToModel(const double_vect_t & v);
		virtual void				getCurrValuesFromModel(double_vect_t & v) const;
		virtual double_vect_t		listCurrValuesFromModel();
        virtual void                getParams();
        virtual void                setParams(const std::vector<double> & v);

	private:

		RelRatesMove &				operator=(const RelRatesMove &);	// never use - don't define
    };
	
/*----------------------------------------------------------------------------------------------------------------------
|	A SubsetRelRatesMove proposes new subset relative rates that are slightly different than the current rates by  
|   sampling from a Dirichlet distribution with parameters equal to the current rates multiplied by a large value (the 
|   tuning parameter 'psi').
*/
class SubsetRelRatesMove : public DirichletMove
	{
	public:
									SubsetRelRatesMove();
                                    virtual ~SubsetRelRatesMove() {}

		virtual void				educateWorkingPrior();
		double 						recalcWorkingPrior() const;

		virtual void				sendCurrValuesToModel(const double_vect_t & v);
		virtual void				getCurrValuesFromModel(double_vect_t & v) const;
		virtual double_vect_t		listCurrValuesFromModel();
        virtual void                getParams();
        virtual void                setParams(const std::vector<double> & v);
		void						setPartitionModel(PartitionModelShPtr m);

	private:

		PartitionModelShPtr			partition_model;
		SubsetRelRatesMove &		operator=(const SubsetRelRatesMove &);	// never use - don't define
    };

} // namespace phycas

#include "phycas/src/dirichlet_move.inl"

#endif
