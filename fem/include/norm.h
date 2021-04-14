#include <deal.II/base/symmetric_tensor.h>

namespace fem {
    using namespace dealii;
	// deal does have a native function for the norm of a FullMatrix but 
	// not of a SymmetricTensor of rank 2
	template <int dim>
	double norm (SymmetricTensor<2,dim> tensor) 
	{
		double tmp = 0;
		for (unsigned int i=0; i<dim; ++i)
		{
			for (unsigned int j=0; j<dim; ++j)
			{
				tmp += sqrt(tensor[i][j]*tensor[i][j]);
			}
		}
		return tmp;
	}
}