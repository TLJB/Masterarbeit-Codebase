#ifndef SMALL_STRAIN_H
#define SMALL_STRAIN_H


#include<deal.II/lac/full_matrix.h>
#include<deal.II/numerics/vector_tools.h>
#include"kronecker.h"
#include"bodyforce.h"
#include"pointhistory.h"
#include"getstrain.h"
#include"norm.h"

namespace SSLinHard{

    using namespace dealii;
    using namespace fem;

	// Material contains all material parameters and solves the local 
	// System of equations at each quadrature point	
	template<int dim,int spacedim>
	class Material
	{
		public:
		Material();
		~Material();
		// calc_cell_matrix returns the local stiffnes matrix at each 
		// quadrature point in a format that deal can use. I.e. a
		// dofs_per_cell*dof_per_cell matrix
		FullMatrix<double> calc_cell_matrix(FESystem<dim,spacedim> &fe,						
							  typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
							  QGauss<dim>  quadrature_formula);
		// Similar to calc_cell_matrix, calculates the right hand side 
		// per cell in a dofs_per_cell*1 format							  
		Vector<double> calc_cell_rhs(FESystem<dim,spacedim> &fe,
									typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
									QGauss<dim>  quadrature_formula);
									
		SymmetricTensor<2,dim> calc_stress(SymmetricTensor<2,dim>);
		
		private:
		// predictor_corrector checks wether the step is elastic or 
		// plastic calculates the algorithmic tangenr matrix and 
		// internal state variable and saves these	
		void predictor_corrector(FEValues<dim,spacedim> &fe_values,
							  typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
							  QGauss<dim>  quadrature_formula,
							  unsigned int q_point);	
		// E_tensor is the classic linear elasticity stress-strain 
		// tensor whereas stress-strain-tensor is either the E_tensor or
		// the algorithmic tangent - calculation in predictor corrector							  				  
		SymmetricTensor<4,dim> E_tensor;
		SymmetricTensor<4,dim> stress_strain_tensor;
		double lambda;
		double mu;
		double yield_stress = 0.3;
		double hardening = 1000;
		double youngs_modulus = 210;
		double poisson_ratio = 0.33;
	};


// Material Functions --------------------------------------------------

	// the constructor calculates the E_tensor
	template <int dim, int spacedim>
	Material<dim,spacedim>::Material()
	{
		lambda = (youngs_modulus*poisson_ratio)/((1+poisson_ratio)*(1-2*poisson_ratio));
		mu	 = youngs_modulus/(2*(1+poisson_ratio));
		for (unsigned int i=0; i<dim;++i)
		{
			for (unsigned int j=0; j<dim;++j)
			{
				for (unsigned int k=0; k<dim;++k)
				{
					for (unsigned int l=0; l<dim;++l)
					{
					 E_tensor[i][j][k][l] = lambda*kronecker(i,j)*kronecker(k,l)
											+ mu*(kronecker(i,k)*kronecker(j,l) 
											+ kronecker(i,l)*kronecker(j,k));
					}
				}
			}
		}
	}

	template<int dim, int spacedim>
	Material<dim,spacedim>::~Material()
	{}

	// calculates the cell matrix,
	// as input the fe object, the current cell, and the 
	// quadrature_formula are needed	
	template<int dim, int spacedim>
	FullMatrix<double> Material<dim,spacedim>::calc_cell_matrix(FESystem<dim,spacedim> &fe,
									typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
									QGauss<dim>  quadrature_formula)
	{
		// the constructor of fe_values needs fe & quadrature_formula,
		// the other arguments determine which values are calculated by 
		// fe_values(reinit)		
		FEValues<dim,spacedim> fe_values(fe, quadrature_formula,
								update_values | update_gradients | 
								update_quadrature_points | update_JxW_values);
		const unsigned int dofs_per_cell = fe.dofs_per_cell;
		const unsigned int n_q_points = quadrature_formula.size();
		FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);		

		// fe_values.reinit calculates the displacement, gradient, 
		// Jacobian etc.											   
		fe_values.reinit(cell);		
		
		// loop over all quadrature_points
		for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
		{	
			// predictor corrector checks wether the step is elastic or
			// plastic
			predictor_corrector(fe_values,cell,quadrature_formula,q_point);	
			
			// loop over dofs - note that i=0 is the first dof of the 
			// first node and i=1 is the second dof of the first node
			// deal needs this format to rearrange the cell_matrix into
			// the global stiffness matrix
			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{	
					// assemble the cell matrix
					const SymmetricTensor<2,dim>
					strain_i = get_strain (fe_values, i, q_point),
					strain_j = get_strain (fe_values, j, q_point);
					cell_matrix(i,j) = strain_i*stress_strain_tensor*strain_j*fe_values.JxW(q_point);
					
					// There are a number of Assert type functions in 
					// deal to catch errors 
					AssertIsFinite(cell_matrix(i,j));
				}
			}
		}
		
		return cell_matrix;												   
	}
	
	// calc_cell_rhs is very similiar to calc_cell_matrix - both could 
	// be merged into one function
	template<int dim, int spacedim>
	Vector<double> Material<dim,spacedim>::calc_cell_rhs(FESystem<dim,spacedim> &fe,
									typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
									QGauss<dim>  quadrature_formula)
	{
		// See calc_cell_matrix
		FEValues<dim,spacedim> fe_values(fe, quadrature_formula,
								update_values | update_gradients | 
								update_quadrature_points | update_JxW_values);
		fe_values.reinit(cell);						
		const unsigned int n_q_points = quadrature_formula.size();
		const unsigned int dofs_per_cell = fe.dofs_per_cell;
		Vector<double> cell_rhs(dofs_per_cell);
		
		// this section gets the bodyforce values 
		BodyForce<dim> body_force;
		std::vector<Vector<double> > body_force_values (n_q_points,
														Vector<double>(dim));
		body_force.vector_value_list (fe_values.get_quadrature_points(),
										body_force_values);

		// assemble rhs by looping over dofs and q_points
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
			const unsigned int
			component_i = fe.system_to_component_index(i).first;
			for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
			{
				cell_rhs(i) += body_force_values[q_point](component_i) *
								fe_values.shape_value (i,q_point)
								*	
								fe_values.JxW (q_point);
				AssertIsFinite(cell_rhs(i));				
			}   
		}
		
		return cell_rhs;
	}
	
	template<int dim, int spacedim>
	void Material<dim,spacedim>::predictor_corrector(FEValues<dim,spacedim> &fe_values,
							  typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
							  QGauss<dim>  quadrature_formula,
							  unsigned int q_point)
	{
		SymmetricTensor<2,dim> strain;
		SymmetricTensor<2,dim> strain_pl;
		SymmetricTensor<2,dim> trial_stress;
		SymmetricTensor<2,dim> trial_stress_dev;
		double alpha;
		double phi;
		double delta_lambda;
		unsigned int dofs_per_cell;
		PointHistory<dim> *local_quadrature_points_history
			= reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
		dofs_per_cell = quadrature_formula.size();	
		
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
			strain +=  get_strain (fe_values, i, q_point);
		}	
		strain_pl = local_quadrature_points_history[q_point].strain_pl;
		alpha = local_quadrature_points_history[q_point].alpha;
		trial_stress = E_tensor*(strain - strain_pl);
		for (unsigned int i=0; i<dim; ++i)
		{
			for (unsigned int j=0; j<dim; ++j)
			{
				trial_stress_dev[i][j] = trial_stress[i][j] 
										- 1/3*trace(trial_stress)
										*kronecker(i,j);
			}
		}
		phi = norm(trial_stress_dev) - yield_stress - hardening*alpha;
		if (phi <= 0)
			stress_strain_tensor = E_tensor;
		else
		{
			SymmetricTensor<4,dim> E_algo_tensor;
			SymmetricTensor<2,dim> direction;
			SymmetricTensor<2,dim> direction_inverse;
			double nenner = hardening;
			direction = trial_stress_dev/norm(trial_stress_dev);
			direction_inverse = invert(direction);
			for (unsigned int i=0;i<dim;++i)
			{	  
				for (unsigned int j=0; j<dim; ++j)
				{
					for (unsigned int k=0; k<dim; ++k)
					{
						for (unsigned int l=0; l<dim; ++l)
						{					  					  
							nenner += 1/3 * direction_inverse[i][j]
									*E_tensor[i][j][k][l]*direction[k][l]; 
							
						}
					}
				}
			}
			AssertIsFinite(nenner);			
			for (unsigned int i=0; i<dim; ++i)
			{	  
				for (unsigned int j=0; j<dim; ++j)
				{
					for (unsigned int k=0; k<dim; ++k)
					{
						for (unsigned int l=0; l<dim; ++l)
						{	
							for (unsigned int m=0; k<dim; ++k)
							{
								for (unsigned int n=0; n<dim; ++n )				  					  
								{	
									for (unsigned int o=0; k<dim; ++k)
									{
										for (unsigned int p=0; n<dim; ++n )				  					  
										{	
											E_algo_tensor[m][n][o][p] = E_tensor[m][n][o][p] + (E_tensor[i][j][k][l]*direction[k][l]*E_tensor[i][j][m][n])*direction[o][p]
																		/nenner;						
										}
									}
								}
							}
						}
					}
				}
			}
			stress_strain_tensor = E_algo_tensor;	
			delta_lambda = phi/nenner;	
			local_quadrature_points_history[q_point].strain_pl += delta_lambda*direction;
			local_quadrature_points_history[q_point].alpha += delta_lambda;					
		}	
	}
	
	// calculates the stress, input is the elastic strain
	template<int dim, int spacedim>
	SymmetricTensor<2,dim> Material<dim,spacedim>::calc_stress(SymmetricTensor<2,dim> strain_el)
	{
		return E_tensor*strain_el;
	}

// Material Functions End ----------------------------------------------

}

#endif