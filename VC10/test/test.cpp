#include "test.hpp"
#include <boost/format.hpp>

// ************************* begin modifiable *****************************************
#include <iostream>
#include <vector>

// *************************** end modifiable *****************************************

int main()
	{
	// ************************* begin modifiable *****************************************

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
