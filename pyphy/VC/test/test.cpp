#include "test.hpp"

// ************************* begin modifiable *****************************************
#include "boost/shared_ptr.hpp"

class IntegerHolder
	{
	public:

		IntegerHolder(int v) : value(v) {}
		~IntegerHolder()
			{
			std::cerr << "IntegerHolder(" << value << ") is being destroyed..." << std::endl;
			}

		int value;
	};

typedef boost::shared_ptr<IntegerHolder> IntegerHolderShPtr;
typedef boost::shared_ptr<const IntegerHolder> ConstIntegerHolderShPtr;

void func(IntegerHolderShPtr p)
	{
	ConstIntegerHolderShPtr q = p;
	std::cerr << "Inside func: value of q is " << q->value << std::endl;
	std::cerr << "Now trying to change p..." << std::endl;
	p->value = 2;
	std::cerr << "Inside func: value of p is " << p->value << std::endl;
	std::cerr << "Inside func: value of q is " << q->value << std::endl;
	std::cerr << "Leaving func..." << std::endl;
	}

void run()
	{
	IntegerHolderShPtr p(new IntegerHolder(1));
	std::cerr << "Inside run, before calling func: value of p is " << p->value << std::endl;
	func(p);
	std::cerr << "Inside run, after calling func: value of p is " << p->value << std::endl;
	std::cerr << "Leaving run..." << std::endl;
	}

// *************************** end modifiable *****************************************

int main()
	{
	// ************************* begin modifiable *****************************************
	run();

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
