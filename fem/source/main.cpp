#include<iostream>
#include"toplevel.h"

int main()
{
	using namespace fem;
	// This error catch doesn't make all that much sense
	try 
	{
		TopLevel<2,2> test;
		test.run();
		std::cout << "test run" << std::endl;
	}
	catch (...)
	{
		std:: cout << "execution failed" << std::endl;
		return 1;
	}
		return 0;
}
