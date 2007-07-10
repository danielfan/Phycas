#include "test.hpp"

// ************************* begin modifiable *****************************************
#include <iostream>
#include <process.h>

// *************************** end modifiable *****************************************

int main()
	{
	// ************************* begin modifiable *****************************************

	int retcode = system("python -V");
	if (retcode)
		std::cout << "nopython";
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
