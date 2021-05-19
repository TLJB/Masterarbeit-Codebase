#include<iostream>
#include"toplevel.h"
#include<boost/exception/diagnostic_information.hpp>

int main()
{
	using namespace fem;
	// This error catch doesn't make all that much sense
	try 
	{
		TopLevel<3,3> test;
		test.run();
		std::cout << "test run" << std::endl;
	}
	catch (std::exception &exc)
	{
		std::cerr << "execution failed" << std::endl;
    std::cerr << boost::diagnostic_information(exc) << std::endl;
		return 1;
	}
		return 0;
}
