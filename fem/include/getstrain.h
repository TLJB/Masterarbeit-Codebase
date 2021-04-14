#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/fe_values.h>

namespace fem{

    using namespace dealii;

	// This version of get_strain calculates the strain of dof of a node 
	// at quadrature point - to be used during the assembly of the 
	// system matrix 
	template <int dim>
	SymmetricTensor<2,dim>
	get_strain (const FEValues<dim> &fe_values,
				const unsigned int   shape_func,
				const unsigned int   q_point)
	{
		SymmetricTensor<2,dim> strain;
		for (unsigned int i=0; i<dim; ++i)
		strain[i][i] = fe_values.shape_grad_component (shape_func,q_point,i)[i];
		for (unsigned int i=0; i<dim; ++i)
		for (unsigned int j=i+1; j<dim; ++j)
		strain[i][j]
			= (fe_values.shape_grad_component (shape_func,q_point,i)[j] +
			   fe_values.shape_grad_component (shape_func,q_point,j)[i]) / 2;
		return strain;
	}

	// This version of get_strain uses the solution of the timestep to 
	// to calculate the strain - it can therefore only be used to update 
	// the quadrature point history	
	template <int dim>
	SymmetricTensor<2,dim>
	get_strain (const std::vector<Tensor<1,dim> > &grad)
	{
		Assert (grad.size() == dim, ExcInternalError());
		SymmetricTensor<2,dim> strain;
		for (unsigned int i=0; i<dim; ++i)
		{
			strain[i][i] = grad[i][i];
		}
		for (unsigned int i=0; i<dim; ++i)
		{
			for (unsigned int j=i+1; j<dim; ++j)
			{
				strain[i][j] = (grad[i][j] + grad[j][i]) / 2;
			}
		}
		return strain;
	}	

}