/**
 * @file small-strain.h
 * @author Till Budde (tilljanis.budde@tu-dortmund.de)
 * @brief Implements small strain elasticity + linear hardening
 * @version 0.1
 * @date 2021-06-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef SMALL_STRAIN_H
#define SMALL_STRAIN_H


#include<deal.II/lac/full_matrix.h>
#include<deal.II/numerics/vector_tools.h>
#include<deal.II/base/tensor.h>
#include"bodyforce.h"
#include"pointhistory.h"

/**
 * @brief Namespace for small-strain elasticity + linear-hardening
 * 
 */
namespace SSLinEla{

    using namespace dealii;
    using namespace fem;

  /**
   * @brief Small Strain Linear Elasticity + Linear Hardening Material Model
   * 
   * @tparam dim number of dimensions 
   * @tparam spacedim number of spatial dimensions
   */
	template<int dim,int spacedim>
	class Material
	{
		public:
		/**
		 * @brief Construct a new Material object
		 * 
		 */
		Material();
		/**
		 * @brief Destroy the Material object
		 * 
		 */
		~Material();
		/**
		 * @brief Calculate the elemntal stiffness matrix
		 * 
		 * @param fe FESystem object (masterElement)
		 * @param cell pointer to cell
		 * @param quadrature_formula quadrature formula object
		 * @param Ue displacement of element nodes
		 * @return FullMatrix<double> elemental stiffness matrix
		 */
		FullMatrix<double> calc_cell_matrix(FESystem<dim,spacedim> &fe,						
							  typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
							  QGauss<dim>  quadrature_formula,
                Vector<double> Ue);
		/**
		 * @brief Calculate the elemntal right hand side vector
		 * 
		 * @param fe FESystem object (masterElement)
		 * @param cell pointer to cell
		 * @param quadrature_formula quadrature formula object
		 * @param Ue displacement of element nodes
		 * @return Vector<double> elemental right-hand-side vector
		 */
		Vector<double> calc_cell_rhs(FESystem<dim,spacedim> &fe,
									typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
									QGauss<dim>  quadrature_formula,
                  Vector<double> Ue);
									
		/**
		 * @brief 
		 * 
		 * @param fe 
		 * @param cell 
		 * @param quadrature_formula 
		 * @param Ue 
		 * @return SymmetricTensor<2,dim> 
		 */
		SymmetricTensor<2,dim> calc_stress(FESystem<dim,spacedim> &fe,
									typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
									QGauss<dim>  quadrature_formula,
                  Vector<double> Ue,
									unsigned int q_point);

		private:
		// E_tensor is the classic linear elasticity stress-strain 
		// tensor whereas stress-strain-tensor is either the E_tensor or
		// the algorithmic tangent - calculation in predictor corrector							  				  
		SymmetricTensor<4,dim> E_tensor;
		double lambda;
		double mu;
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
						E_tensor[i][j][k][l] = lambda*(i==j)*(k==l)
																	+ mu*( (i==k)*(j==l) + (i==l)*(j==k));
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
									QGauss<dim>  quadrature_formula,
                  Vector<double> Ue)
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
		unsigned int nodes_per_cell = cell->n_vertices();

		// fe_values.reinit calculates the displacement, gradient, 
		// Jacobian etc.											   
		fe_values.reinit(cell);		
		
		// loop over all quadrature_points
		for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
		{	
			
			// loop over dofs - note that i=0 is the first dof of the 
			// first node and i=1 is the second dof of the first node
			// deal needs this format to rearrange the cell_matrix into
			// the global stiffness matrix
			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
		    const unsigned int component_i =
 			     fe.system_to_component_index(i).first;
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{	
					const unsigned int component_j =
						fe.system_to_component_index(j).first;
					// assemble the cell matrix
					cell_matrix(i,j) += (fe_values.shape_grad(i,q_point)
														*E_tensor
														*fe_values.shape_grad(j,q_point)
														*fe_values.JxW(q_point))[component_i][component_j];
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
									QGauss<dim>  quadrature_formula,
                  Vector<double> Ue)
	{
		// See calc_cell_matrix
		FEValues<dim,spacedim> fe_values(fe, quadrature_formula,
								update_values | update_gradients | 
								update_quadrature_points| update_inverse_jacobians | update_JxW_values);
		fe_values.reinit(cell);						
		const unsigned int n_q_points = quadrature_formula.size();
		const unsigned int dofs_per_cell = fe.dofs_per_cell;
		Vector<double> cell_rhs(dofs_per_cell), finte(dofs_per_cell), fvole(dofs_per_cell);
		
		// this section gets the bodyforce values 
		BodyForce<dim> body_force;
		std::vector<Vector<double> > body_force_values (n_q_points,
														Vector<double>(dim));
		body_force.vector_value_list (fe_values.get_quadrature_points(),
										body_force_values);

		// // assemble rhs by looping over dofs and q_points
		for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
		{

			SymmetricTensor<2,dim> strain; strain=0;
			for (unsigned int i=0; i<dofs_per_cell; ++i){
				const unsigned int
				component_i = fe.system_to_component_index(i).first;
				for (unsigned int j=0; j<dim; ++j) {

					strain[component_i][j] += Ue[i]*(fe_values.shape_grad(i,q_point))[j];

				}

			}
			strain = 0.5*(strain + transpose(strain));

			SymmetricTensor<2,dim> stress = E_tensor*strain;

			finte=0;
			fvole=0;
			for (unsigned int i=0; i<dofs_per_cell; ++i)	{
				const unsigned int
				component_i = fe.system_to_component_index(i).first;
					fvole(i) += - body_force_values[q_point](component_i) 
									*	fe_values.shape_value (i,q_point)
									*	fe_values.JxW (q_point);
					finte(i) += (fe_values.shape_grad(i,q_point)
											* stress
											* fe_values.JxW(q_point))[component_i];
				}   
				cell_rhs.add(+1,finte, -1, fvole);
				for (auto rhs_val : cell_rhs)
					AssertIsFinite(rhs_val);				
		}
		return cell_rhs;
	}

	template<int dim, int spacedim>
	SymmetricTensor<2,dim> Material<dim,spacedim>::calc_stress(FESystem<dim,spacedim> &fe,
								typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
								QGauss<dim>  quadrature_formula,
								Vector<double> Ue,
								unsigned int q_point)
	{

		// See calc_cell_matrix
		FEValues<dim,spacedim> fe_values(fe, quadrature_formula,
								update_values | update_gradients | 
								update_quadrature_points| update_inverse_jacobians | update_JxW_values);
		fe_values.reinit(cell);						
		const unsigned int n_q_points = quadrature_formula.size();
		const unsigned int dofs_per_cell = fe.dofs_per_cell;
		Vector<double> cell_rhs(dofs_per_cell), finte(dofs_per_cell), fvole(dofs_per_cell);
		
		// this section gets the bodyforce values 
		BodyForce<dim> body_force;
		std::vector<Vector<double> > body_force_values (n_q_points,
														Vector<double>(dim));
		body_force.vector_value_list (fe_values.get_quadrature_points(),
										body_force_values);

		// // assemble rhs by looping over dofs and q_points
		SymmetricTensor<2,dim> strain; strain=0;
		for (unsigned int i=0; i<dofs_per_cell; ++i){
			const unsigned int
			component_i = fe.system_to_component_index(i).first;
			for (unsigned int j=0; j<dim; ++j) {

				strain[component_i][j] += Ue[i]*(fe_values.shape_grad(i,q_point))[j];

			}

		}
		strain = 0.5*(strain + transpose(strain));

		SymmetricTensor<2,dim> stress = E_tensor*strain;
		return stress;
	}
	
// Material Functions End ----------------------------------------------

}

#endif