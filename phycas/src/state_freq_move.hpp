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
|	A StateFreqMove proposes new state frequencies that are slightly different than the current state frequencies by 
|   sampling from a Dirichlet distribution with parameters equal to the current frequencies multiplied by a large 
|   value (the tuning parameter 'psi').
*/
class StateFreqMove : public MCMCUpdater
	{
	public:
									StateFreqMove();
									virtual ~StateFreqMove() 
										{
										//std::cerr << "StateFreqMove dying..." << std::endl;
										}


		// Accessors
		double						getPsi() const;

		// Modifiers
		void						setPsi(double x);

		// Utilities
		void						reset();

		// These are virtual functions in the MCMCUpdater base class
		virtual bool				update();
		virtual double				getLnHastingsRatio() const;
		virtual double				getLnJacobian() const;
		virtual void				proposeNewState();
		virtual void				revert();
		virtual void				accept();

	private:

		StateFreqMove &				operator=(const StateFreqMove &);	// never use - don't define

	private:

		double						psi;			/**< Larger values result in changes of smaller magnitude */
		std::vector<double> 		new_freqs;	    /**< Proposed new frequencies */
		std::vector<double> 		orig_freqs;	    /**< Saved frequencies (in case revert is necessary) */
		std::vector<double> 		c_forward;	    /**< Dirichlet parameter vector used to propose new frequencies */
		std::vector<double> 		c_reverse;	    /**< Dirichlet parameter vector used to propose original frequencies (only used to compute Hastings ratio) */
        DirichletShPtr              dir_forward;    /**< Points to an ad hoc Dirichlet distribution object used to assess the forward move density and to propose a new frequency vector */
        DirichletShPtr              dir_reverse;    /**< Points to an ad hoc Dirichlet distribution object used to assess the reverse move density */
	};

} // namespace phycas

#include "phycas/src/edge_move.inl"

#endif
