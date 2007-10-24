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

#if ! defined(NCAT_MOVE_HPP)
#define NCAT_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
//#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "phycas/src/mcmc_updater.hpp"		// for base class MCMCUpdater

namespace phycas
{

class ProbabilityDistribution;
typedef boost::shared_ptr<ProbabilityDistribution>	ProbDistShPtr;

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the Larget-Simon local move.
*/
class NCatMove : public MCMCUpdater
	{
	public:
										NCatMove();
										virtual ~NCatMove() 
											{
											//std::cerr << "NCatMove dying..." << std::endl;
											}

		// Accessors
		//
		bool							addCatMoveProposed() const;

		// Modifiers
		//
		double							getPhi() const;
		void							setPhi(double new_phi);

		unsigned						getNCatMax() const;
		void							setNCatMax(unsigned new_ncat_max);

		double							getLambda() const;
		void							setLambda(double new_lambda);

		double							getL() const;
		void							setL(double new_L);

		unsigned						getS() const;
		void							setS(unsigned new_s);

		ProbDistShPtr					getCatProbPrior() const;
		void							setCatProbPrior(ProbDistShPtr new_cat_prob_prior);

		// Utilities
		//

		// These are virtual functions in the MCMCUpdater base class
		//
		virtual bool					update();
		virtual double					getLnHastingsRatio() const;
		virtual double					getLnJacobian() const;
		virtual void					proposeNewState();
		virtual void					revert();
		virtual void					accept();

	private:

		void							reset();
		void							proposeAddCatMove();
		void							proposeDelCatMove();
		NCatMove &						operator=(const NCatMove &);	// never use - don't define

	private:

		unsigned						ncat_max;				/**< Maximum number of rate categories encountered. If this value is exceeded, the tree must be re-equipped (transition matrices and conditional likelihood arrays need to be reallocated) to accommodate the increased number of rate categories */
		double							lambda;					/**< Parameter of the poisson distribution determining the ncat prior */
		ProbDistShPtr					cat_prob_prior;			/**< Probability distribution being used for the prior of category probabilities */
		unsigned						s;						/**< Number of imaginary "spacer" rates between each actual rate */
		double							u1;						/**< Uniform(0,1) random deviate chosen during addcat move to select position of new rate (new rate is at u1*L) */
		std::vector<double>::iterator	tmp_rate_iter;			/**< Iterator used to insert the new rate in an addcat move (and to delete the new rate again if the move is reverted) */		
		std::vector<double>::iterator	tmp_prob_iter;			/**< Iterator used to insert the new probability in an addcat move (and to delete the new probability again if the move is reverted) */
		double							tmp_prob;				/**< Holds new probability chosen during an addcat move (in case revert is necessary) */
		double							tmp_rate;				/**< Holds new rate chosen during an addcat move (in case revert is necessary) */
		unsigned						ncat_before;			/**< Holds number of rate categories before the proposed move */
		unsigned						ncat_after;				/**< Holds number of rate categories after the proposed move */
		double							rate_left;				/**< Holds value of rate to the left of the inserted rate after an addcat move */
		double							rate_right;				/**< Holds value of rate to the right of the inserted rate after an addcat move */
		double							L;						/**< Unnormalized relative rates range from 0 to L */
		double							phi;					/**< Probability that an addcat move is proposed; 1 - phi is the probability that a delcat move will be proposed instead */
		bool							addcat_move_proposed;	/**< True if last move proposed was an addcat move; false otherwise */

		double							ln_jacobian;			/**< The natural log of the Jacobian for the move last proposed */
		double							ln_hastings;			/**< The natural log of the Hastings ratio for the move last proposed */

		unsigned total_updates;
		std::vector<unsigned> 			ncat_distr;

		CDF								cdf;					/**< CDF object needed for its LnGamma function */
	};

} // namespace phycas

#include "phycas/src/ncat_move.inl"

#endif
