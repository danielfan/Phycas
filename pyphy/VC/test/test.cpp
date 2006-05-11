#include "test.hpp"

// ************************* begin modifiable *****************************************
#include <vector>

// *************************** end modifiable *****************************************

int main()
	{
	// ************************* begin modifiable *****************************************
	double lower = 0.0;
	double upper = 1.0;
	std::vector<double> v;
	v.push_back(0.3);
	v.push_back(0.6);
	v.push_back(0.8);

	// (0.1)  0.3         0.6        0.8          lower_bound = begin()  0.0,     *it
	//        0.3  (0.4)  0.6        0.8          lower_bound = 0.6      *(--it), *it
	//        0.3         0.6 (0.7)  0.8          lower_bound = 0.8      *(--it), *it
	//        0.3         0.6        0.8  (0.9)   lower_bound = end()    *(--it), L

	double val = 0.0;

	std::vector<double>::iterator p;

	for (p = v.begin(); p != v.end(); ++p)
		std::cerr << *p << "  " << std::endl;

	val = 0.1;
	p = std::lower_bound(v.begin(), v.end(), val);
	std::cerr << "val = " << val << std::endl;
	if (p == v.begin())
		std::cerr << "p = begin()" << std::endl;
	else if (p == v.end())
		std::cerr << "p = end()" << std::endl;
	else
		std::cerr << "p = " << *p << std::endl;

	val = 0.4;
	p = std::lower_bound(v.begin(), v.end(), val);
	std::cerr << "val = " << val << std::endl;
	if (p == v.begin())
		std::cerr << "p = begin()" << std::endl;
	else if (p == v.end())
		std::cerr << "p = end()" << std::endl;
	else
		std::cerr << "p = " << *p << std::endl;

	val = 0.7;
	p = std::lower_bound(v.begin(), v.end(), val);
	std::cerr << "val = " << val << std::endl;
	if (p == v.begin())
		std::cerr << "p = begin()" << std::endl;
	else if (p == v.end())
		std::cerr << "p = end()" << std::endl;
	else
		std::cerr << "p = " << *p << std::endl;

	val = 0.9;
	p = std::lower_bound(v.begin(), v.end(), val);
	std::cerr << "val = " << val << std::endl;
	if (p == v.begin())
		std::cerr << "p = begin()" << std::endl;
	else if (p == v.end())
		std::cerr << "p = end()" << std::endl;
	else
		std::cerr << "p = " << *p << std::endl;

#if 0
	typedef std::pair<std::vector<double>::iterator, std::vector<double>::iterator> MyPair;
	MyPair p = std::equal_range(v.begin(), v.end(), val);
	std::cerr << "val = " << val << std::endl;
	if (p.first == v.end() && p.second == v.end())
		{
		std::cerr << "first  = end()" << std::endl;
		std::cerr << "second = end()" << std::endl;
		//std::cerr << "prev = " << v[v.size() - 1] << std::endl;
		//std::cerr << "next = " << upper << std::endl;
		}
	else if (p.first != v.end() && p.second == v.end())
		{
		std::cerr << "first  = " << *(p.first) << std::endl;
		std::cerr << "second = end()" << std::endl;
		//std::cerr << "prev = " << *(p.first) << std::endl;
		//std::cerr << "next = " << upper << std::endl;
		}
	else if (p.first != v.end() && p.second != v.end())
		{
		if (p.first == p.second)
			{
			if (p.first == v.begin())
				{
				std::cerr << "first  = begin()" << std::endl;
				std::cerr << "second = begin()" << std::endl;
				//std::cerr << "prev = " << lower << std::endl;
				//std::cerr << "next = " << *(p.second) << std::endl;
				}
			else
				{
				std::cerr << "first  = " << *(p.first) << std::endl;
				std::cerr << "second = " << *(p.second) << std::endl;
				//std::cerr << "prev = " << *(--p.first) << std::endl;
				//std::cerr << "next = " << *(p.second) << std::endl;
				}
			}
		else
			{
			std::cerr << "first  = " << *(p.first) << std::endl;
			std::cerr << "second = " << *(p.second) << std::endl;
			//std::cerr << "prev = " << *(p.first) << std::endl;
			//std::cerr << "next = " << *(p.second) << std::endl;
			}
		}
	else
		std::cerr << "oops!" << std::endl;
#endif		

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
