#include "test.hpp"

// ************************* begin modifiable *****************************************
#include <cstdlib>
#include <iostream>
// *************************** end modifiable *****************************************

int main()
	{
	// ************************* begin modifiable *****************************************
	std::cout << "Python 2.4.1" << std::endl;
	std::exit(0);

	// *************************** end modifiable *****************************************

	char answer = AskUser("\nPress y to quit...");
	if (answer == '\0')
		std::cerr << "Oops, something bad happened in the AskUser function." << std::endl;
	else if (answer == 'y')
		std::cerr << "\nYou pressed y so I'm quitting" << std::endl;
	else
		std::cerr << "\nYou pressed " << answer << " but I'm quitting anyway!" << std::endl;
 
	return 0;
	}
