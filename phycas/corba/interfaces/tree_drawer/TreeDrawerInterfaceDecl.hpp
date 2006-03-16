#if ! defined(TREE_DRAWER_INTERFACE_DECL_HPP)
#define TREE_DRAWER_INTERFACE_DECL_HPP

namespace cipresCORBA
{

class TreeDrawerInterface
	{
	public:
		virtual char * execute(const char * command) = 0;
		
		virtual void drawTree (const char * newickTree) = 0;
	};

} //namespace cipresCORBA

#endif

    