/**
 * @file main.cpp
 * @author Till Budde (tilljanis.budde@tu-dortmund.de)
 * @brief Main File of the FEM Code
 * @version 0.1
 * @date 2021-06-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include<iostream>
#include"toplevel.h"
#include<boost/exception/diagnostic_information.hpp>

/**
 * @brief Execute the FEM Code
 * 
 * @return int 0 if succesfull, else 1
 */
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
	catch (std::exception &exc)
	{
		std::cerr << "execution failed" << std::endl;
    std::cerr << boost::diagnostic_information(exc) << std::endl;
		return 1;
	}
		return 0;
}
