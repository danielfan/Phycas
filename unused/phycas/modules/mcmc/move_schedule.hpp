#if !defined (MOVE_SCHEDULE_HPP) 
#define MOVE_SCHEDULE_HPP
#include <boost/shared_ptr.hpp>
#include "ncl/output/temp_basic_output_operators.hpp"
class MCMCMove;
class Lot;
class LargetSimonMove;
class BushMaster;
class EdgeMove;
class TreeScalerMove;
typedef boost::shared_ptr<Lot> LotShPtr;
typedef boost::shared_ptr<LargetSimonMove> LargetSimonMoveShPtr;
typedef boost::shared_ptr<BushMaster> BushMasterShPtr;
typedef boost::shared_ptr<EdgeMove> EdgeMoveShPtr;
typedef boost::shared_ptr<TreeScalerMove> TreeScalerMoveShPtr;
class Tree;
class MCMCSettings;
class MoveSchedule
	{
	public:
				/** Specifies the index of each possible Metropolis-Hastings move */
		enum TreeMoveIndex
			{
			TreeMove_Local  = 0,	/**< is the index of the Larget-Simon LOCAL move */
			TreeMove_Bush   = 1,	/**< is the index of the Bush move, which attempts to either add or remove an edge */
			TreeMove_Scaler = 2		/**< is the index of the Scaler move, which scales all edge lengths up or down proportionally */
			};
		enum MoveCategories
			{
			kMetropHastings =   0x01,
			kGibbsMove		=   0x02
			};

		MoveSchedule(const MCMCSettings & settings)
			:vProbTopoMoves()
			{
			Reset(settings);
			}
		bool			DoGibbsThisRound(unsigned ts) const;
		int				GetMoveCategory(unsigned ts) const;
		MCMCMove	  * GetNextTreeMove(Lot & r, const Tree & tree) const;
		void			Reset(const MCMCSettings & settings);
		void			ResetMoves(Tree * tree, LotShPtr r);

		std::vector<double>		vProbTopoMoves;
			// should be protected
		LargetSimonMoveShPtr	lsm;					/**< handles proposals that change branch lengths and topology */
		BushMasterShPtr			bushMaster;			/**< proposes dimension-changing moves (i.e. removal or insertion of edges) to trees */
		EdgeMoveShPtr			edgeMove;				/**< used on the star tree when LargetSimonMove is not possible; just changes one edge at a time */
		TreeScalerMoveShPtr		treeScalerMove;		/**< attempt to scale all edges in tree up or down proportionally */
		bool					doGibbs;
		unsigned				gibbsEvery;
	};

inline bool MoveSchedule::DoGibbsThisRound(unsigned ts) const
	{
	return (doGibbs && (ts % gibbsEvery == 0));
	}

inline int MoveSchedule::GetMoveCategory(unsigned ts) const
	{
	if (!vProbTopoMoves.empty())
		return (DoGibbsThisRound(ts) ? (kMetropHastings | kGibbsMove) : kMetropHastings);
	return  (DoGibbsThisRound(ts) ? kGibbsMove : 0);
	}
	
template<int FORMAT_HINT, class OUT_STREAM>
class GenericPrinterClass<FORMAT_HINT, MoveSchedule, OUT_STREAM>
	{
	public:
		GenericPrinterClass(OUT_STREAM & outStream, const MoveSchedule & ms)
		{
		if (!ms.vProbTopoMoves.empty())
			{
			assert(ms.vProbTopoMoves.size() >  2); 
			std::string tmp;
			outStream << "\n\nTree move probabilities:\n";
			StrPrintF(tmp, "  %.2f  Larget-Simon LOCAL move\n", ms.vProbTopoMoves[MoveSchedule::TreeMove_Local]);
			StrPrintF(tmp, "  %.2f  Bush move\n", ms.vProbTopoMoves[MoveSchedule::TreeMove_Bush]);
			StrPrintF(tmp, "  %.2f  Scaler move\n", ms.vProbTopoMoves[MoveSchedule::TreeMove_Scaler]);
			outStream << tmp;
			}
		}
	};

#endif

