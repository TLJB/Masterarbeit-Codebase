/**
 * @file main.cpp
 * @author Till Budde (tilljanis.budde@tu-dortmund.de)
 * @brief Main File of the FEM Code
 * @version 0.1
 * @date 16.09.2021
 * 
 * @copyright Copyright (c) 2021
 * 
 * This is the main file of the FEM Code for the master thesis 
 * "Implementation und Analyse generalisierter Koh√§sivzonenelemente in die Finite-Elemente- Bibliothek deal.II".
 * It's only use is to create the fem::TopLevel Object and execute the program.
 * 
 */
#include<iostream>
#include"toplevel.h"
#include "timer.h"
#include<boost/exception/diagnostic_information.hpp>

/**
 * @brief Create and run the fem::TopLevel Object
 * 
 * @return int 0 if successful, else 1
 */
int main()
{
	using namespace fem;
	
	try 
	{
		// the template parameter decide wether the calculation is 2D or 3D
		// and wether a numerical cell_tangent is calculated (only interface)
		Timer t;
		TopLevel<2,2, true> fem;
		fem.run();
		std::cout << "Programm completed in " << t.elapsed() << "seconds" << std::endl;
	}
	catch (std::exception &exc)
	{
		std::cerr << "execution failed" << std::endl;
    std::cerr << boost::diagnostic_information(exc) << std::endl;
		return 1;
	}
		return 0;
}
