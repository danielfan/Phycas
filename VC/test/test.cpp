#include "test.hpp"
#include <boost/format.hpp>

// ************************* begin modifiable *****************************************
#include <iostream>
#include <vector>

// *************************** end modifiable *****************************************

int main()
	{
	// ************************* begin modifiable *****************************************
    std::cerr << "sizeof(int)    = " << sizeof(int) << std::endl;
    std::cerr << "sizeof(int8_t) = " << sizeof(signed char) << std::endl;
    std::cerr << "sizeof(float)  = " << sizeof(float) << std::endl;
    std::cerr << "sizeof(double) = " << sizeof(double) << std::endl;

    unsigned n = 100;
    unsigned s = 1000;
    unsigned p = 500;
    unsigned k = 61;

    double nn = (double)n;
    double ss = (double)s;
    double pp = (double)p;
    double kk = (double)k;

    double std_bytes    = 8.0*(nn - 2.0)*kk*pp + 8.0*(2.0*nn - 3.0)*kk*kk;
    double std_kbytes = std_bytes/1024.0;
    double std_mbytes = std_kbytes/1024.0;

    double unimap_bytes = 5.0*(2.0*nn - 3.0)*ss;
    double unimap_kbytes = unimap_bytes/1024.0;
    double unimap_mbytes = unimap_kbytes/1024.0;

    double m = std_bytes/unimap_bytes;

    std::cerr << "number of taxa                   = " << n << std::endl;
    std::cerr << "number of patterns (traditional) = " << p << std::endl;
    std::cerr << "number of states (unimap)        = " << s << std::endl;
    std::cerr << "number of states                 = " << k << std::endl;
    std::cerr << std::endl;
    std::cerr << "traditional (mbytes) = " << std_bytes << std::endl;
    std::cerr << "traditional (kbytes) = " << std_kbytes << std::endl;
    std::cerr << "traditional (mbytes) = " << std_mbytes << std::endl;
    std::cerr << std::endl;
    std::cerr << "unimap (bytes)       = " << unimap_bytes << std::endl;
    std::cerr << "unimap (kbytes)      = " << unimap_kbytes << std::endl;
    std::cerr << "unimap (mbytes)      = " << unimap_mbytes << std::endl;
    std::cerr << std::endl;
    std::cerr << "m           = " << m << std::endl;

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
