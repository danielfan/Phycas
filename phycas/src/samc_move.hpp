/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis						  |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford	  |
|																			  |
|  This program is free software; you can redistribute it and/or modify		  |
|  it under the terms of the GNU General Public License as published by		  |
|  the Free Software Foundation; either version 2 of the License, or		  |
|  (at your option) any later version.										  |
|																			  |
|  This program is distributed in the hope that it will be useful,			  |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of			  |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the			  |
|  GNU General Public License for more details.								  |
|																			  |
|  You should have received a copy of the GNU General Public License along	  |
|  with this program; if not, write to the Free Software Foundation, Inc.,	  |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.				  |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#if defined(SAMC_ONE)

#if ! defined(SAMC_MOVE_HPP)
#define SAMC_MOVE_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
//#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "phycas/src/mcmc_updater.hpp"		// for base class MCMCUpdater

//struct CIPRES_Matrix;

//class ProbabilityDistribution;

namespace phycas
{

class ExponentialDistribution;
typedef boost::shared_ptr<ExponentialDistribution>	ExponentialDistributionShPtr;

class MCMCChainManager;
typedef boost::weak_ptr<MCMCChainManager>			ChainManagerWkPtr;

//class HKY;

typedef std::map<unsigned, std::vector<double> > PolytomyDistrMap;
typedef std::vector<double> VecPolytomyDistr;
typedef std::vector<double> DistanceMatrixRow;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the SAMC dimension-change move
*/
class SamcMove : public MCMCUpdater
	{
	public:
									SamcMove(ProbDistShPtr);
									virtual ~SamcMove() 
										{
										}
		// Modifiers
		//
		void						setEdgeLenDistMean(double mean);

		// Utilities
		//
		//void						finalize();

		bool						extrapolate(unsigned leaf_num, double theta_diff, double ln_proposal_ratio);
		bool						project(unsigned leaf_num, double theta_diff, double ln_proposal_ratio);

		// These are virtual functions in the MCMCUpdater base class
		//
		virtual bool				update();
		virtual void				revert();
		virtual void				accept();

        unsigned                    getNTax();
        void                        setNTax(unsigned ntax);

        void                        setDistanceMatrixRow(const DistanceMatrixRow & row);
        double                      calcLeafEdgeLen(unsigned leaf_k, TreeNode * nd);

        void                        setTemperature(double temp);
        double                      getTemperature() const;

	private:

		void						reset();
		SamcMove &					operator=(const SamcMove &);	// never use - don't define

		void						calcPk(unsigned leaf_k);
		double						getPkl(unsigned leaf_k, TreeNode * nd_l);
		TreeNode *					chooseRandomAttachmentNode(unsigned leaf_k);

    private:

        bool                            infinite_temperature;               /**< if true, uses a discrete uniform distribution for selecting attachment node during extrapolate moves; if false, uses less efficient but more reasonable method based on least squares length of edge of attached leaf */
        const double                    temperature_threshhold;             /**< temperatures above this value trigger setting infinite_temperature to true */
        double                          temperature;                        /**< assuming infinite_temperatore is false, determines the temperature used for computing insertion probabilities for extrapolate moves */
        double                          prev_ln_like;
        std::vector<DistanceMatrixRow>  distance_matrix;
		unsigned						num_taxa;							/**< The number of taxa */
		TreeNode *						leaf_sib;
		TreeNode *						leaf;
		TreeNode *						parent;
		double							leaf_sib_orig_edgelen;				/**<  */
		double							orig_edgelen;						/**< Length of deleted edge saved (in case revert of delete-edge move is necessary) */
		double							par_orig_edgelen;					/**< Length of deleted parent edge saved (in case revert of delete-edge move is necessary) */
		bool							last_move_projection;				/**< true if project was the last move*/
		std::vector<double>				pvect;								/**< vector of pkl values in preorder sequence */
		TreeNode *						new_leaf_sib_parent;				/**< node serving as likelihood root after projection */
		ProbDistShPtr					term_edge_dist;						/**< Probability distribution for proposing new edge lengths during the extrapolate move*/
	};

} // namespace phycas

#include "phycas/src/samc_move.inl"

#endif

#endif //if defined(SAMC_ONE)
