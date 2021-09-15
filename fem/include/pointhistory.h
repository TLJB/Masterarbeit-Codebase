/**
 * @file pointhistory.h
 * @author Till Budde
 * @brief Implements structure to save state dependent variables at quadrature points
 * @version 0.1
 * @date 2021-06-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef POINHISTORY_H
#define POINHISTORY_H

#include <deal.II/base/symmetric_tensor.h>

namespace fem {

  using namespace dealii;

	template<int dim>
	struct PointHistoryBulk {
		SymmetricTensor<2,dim> strain_pl;
		double alpha; 
		SymmetricTensor<2,dim> old_stress;
	};

	template<int dim> 
	struct PointHistoryInter {
		double kappa;
	};
	/**
	 * @brief Save state dependent variables of quadrature point
	 * 
	 * @tparam dim  number of dimensions
	 */
	template<int dim>
	struct PointHistory
	{
		PointHistoryBulk<dim>  bulk;
		PointHistoryInter<dim> inter;
	};

}

#endif