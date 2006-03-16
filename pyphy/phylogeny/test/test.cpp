#include <string>
#include <iostream>
#include "pyphy/phylogeny/basic_tree.hpp"

int main()
	{
	Tree t;
	t.BuildFromString("((a:0.1,b:0.2)x:0.6,c:0.3,(d:0.4, e:0.5)y:0.7)z");
	std::cout << "  " << t.DebugWalkTree() << std::endl;
	std::cout << "  " << t.MakeNewick() << std::endl;

	std::cerr << "Press return to quit" << std::endl;
	int ch = std::cin.get();

	return 0;
	}
