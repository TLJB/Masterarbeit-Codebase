#ifndef POINHISTORY_H
#define POINHISTORY_H

#include <deal.II/base/symmetric_tensor.h>

namespace fem {

    using namespace dealii;
	// PointHistory contains internal variables of the previous timestep
	template<int dim>
	struct PointHistory
	{
		SymmetricTensor<2,dim> strain_pl;
		double alpha; 
		SymmetricTensor<2,dim> old_stress;
	};
}

#endif