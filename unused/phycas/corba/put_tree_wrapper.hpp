#if ! defined (PUT_TREE_WRAPPER_HPP)
#define PUT_TREE_WRAPPER_HPP
#include <boost/function.hpp>

class Tree;
namespace cipresCORBA
{

class CorbaWrapper;
/**
 * Class that insulates code (such as MCMC) which can be used a source of 
 *	trees for a wide variety of analyses from knowing about all possible
 *	CIPRes services that take a const reference to a tree as an argument and
 *	return void (or a value that can be ignored).
 * By using PutTreeWrapper the MCMC code does not have to include headers for 
 *	the service classes (such as TreeDrawer, TreeDatabase, etc.).
 *
 */
class PutTreeWrapper
	{
	public:
		typedef boost::function1<void, const Tree &> TreeAcceptorFunc;
		
		PutTreeWrapper(const VecString & taxLabels);
		
		void AddTreeFunc(const TreeAcceptorFunc & treeFunc)
			{
			vecTreeAcceptorFunc.push_back(treeFunc);
			}
		void PutTree(const Tree &t);
			
	protected:
		typedef std::vector<TreeAcceptorFunc> VecTreeAcceptorFunc;
		
		VecTreeAcceptorFunc vecTreeAcceptorFunc;
	};
} // namespace cipresCORBA
#endif