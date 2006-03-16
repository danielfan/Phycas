#include "test.hpp"

// ************************* begin modifiable *****************************************
//#include <vector>
//#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
//#include <cmath>

void foo(int k)
	{
	std::cerr << "foo(int) called with value " << k << std::endl;
	}

double bar(double d)
	{
	//std::cerr << "bar(double) called with value " << d << std::endl;
	return sqrt(d);
	}

// *************************** end modifiable *****************************************

int main()
	{
	// ************************* begin modifiable *****************************************
	std::vector<double> v;
	v.push_back(1.0);
	v.push_back(2.0);
	v.push_back(3.0);
	v.push_back(4.0);

	std::vector<double> v2(4);

	//std::transform(v.begin(), v.end(), v2.begin(), boost::lambda::bind(psqrt, boost::lambda::_1));
	//std::for_each(v2.begin(), v2.end(), std::cout << boost::lambda::_1 << constant('\n'));

	using namespace boost::lambda;
	//double i = 9.0; 

	//std::for_each(v.begin(), v.end(), std::cerr << bind(&bar, _1) << '\n');
	//std::transform(v.begin(), v.end(), v2.begin(), bind(&bar, _1));
	std::transform(v.begin(), v.end(), v2.begin(), bind(static_cast<double(*)(double)>(&sqrt), _1));


	std::for_each(v2.begin(), v2.end(), std::cout << boost::lambda::_1 << constant('\n'));

	//double sqrt_i = bind(&bar, _1)(i);
	//std::cerr << "Square root of the value " << i << " is " << sqrt_i << std::endl;

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
