/**
 * @file LinElaInter.h
 * @author Till Janis Budde
 * @brief  Linear Elastic Interface Element for testing purposes
 * @version 0.1
 * @date 2021-05-19
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include<deal.II/lac/full_matrix.h>
#include<deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include"bodyforce.h"
#include"pointhistory.h"
#include<boost/exception/diagnostic_information.hpp>
#include<deal.II/base/point.h>
#include<vector>

/**
 * @brief Namespace of linear elastic Interface Element
 * 
 */
namespace LinElaInter {

  using namespace dealii;
  using namespace fem;

	/**
	 * @brief Spring-like interface element
	 * 
	 * @tparam dim number of dimensions
	 * @tparam spacedim  number of spatial dimensions
	 */
	template<int dim,int spacedim>
	class Material {

    public:
		/**
		 * @brief Construct a new Material object
		 * 
		 */
		Material(){};
		/**
		 * @brief Destroy the Material object
		 * 
		 */
		~Material(){};

    /**
     * @brief Calculate Cell Matrix for Interface Element
     * 
     * @param fe FESystem / Finite Element
     * @param cell active_cell_iterator
     * @param quadrature_formula quadrature formula object
     * @param Ue Elemental Displacements
     * @return FullMatrix<double> 
     */
		FullMatrix<double> calc_cell_matrix(FESystem<dim,spacedim> &fe,						
							  typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
							  QGauss<dim>  quadrature_formula,
                Vector<double> Ue
                );

    
    /**
     * @brief Calculate Right Hand Side per cell
     * 
     * @param fe 
     * @param cell 
     * @param quadrature_formula 
     * @return Vector<double> 
     */
		Vector<double> calc_cell_rhs(FESystem<dim,spacedim> &fe,
									typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
									QGauss<dim>  quadrature_formula,
                  Vector<double> Ue
                  );

		/**
		 * @brief Calculate the stress Tensor
		 * 
		 * @param SymmetricTensor<2,dim> Strain Tensor
		 * @return SymmetricTensor<2,dim> Stress Tensor 
		 */
		SymmetricTensor<2,dim> calc_stress(SymmetricTensor<2,dim> strain);

    private:
		double stiffness = 210;
  };

	template<int dim, int spacedim>
	FullMatrix<double> Material<dim,spacedim>::calc_cell_matrix(FESystem<dim,spacedim> &fe,
									typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
									QGauss<dim>  quadrature_formula,
                  Vector<double> Ue
                  )
	{
		// the constructor of fe_values needs fe & quadrature_formula,
		// the other arguments determine which values are calculated by 
		// fe_values(reinit)		
		FEValues<dim,spacedim> fe_values(fe, quadrature_formula,
								update_values | update_gradients | 
								update_quadrature_points);
		const unsigned int dofs_per_cell = fe.dofs_per_cell;
		const unsigned int n_q_points = quadrature_formula.size()/2;
		const unsigned int nodes_per_cell = cell->n_vertices();
		FullMatrix<double> cell_matrix(dofs_per_cell);		

		// fe_values.reinit calculates the displacement, gradient, 
		// Jacobian etc.											   
		try {
			fe_values.reinit(cell);		
		}
		catch (std::exception &exc) {
    	std::cerr << boost::diagnostic_information(exc) << std::endl;
		}
		
    auto identity =  Physics::Elasticity::StandardTensors<2>::I;


		Vector<double> qp(n_q_points), N(nodes_per_cell);
		qp[0] = -1/sqrt(3);
		qp[1] = 1/sqrt(3);

		auto C = stiffness*identity;
		for (unsigned int q_point =0; q_point<n_q_points;++q_point){

			N[0]  = 0.5*(1-qp[q_point]);
			N[1]  = 0.5*(1+qp[q_point]);
			N[2]  = 0.5*(1-qp[q_point]);
			N[3]  = 0.5*(1+qp[q_point]);

			for (unsigned int i=0; i<dofs_per_cell; ++i) {
				const unsigned int 
					component_i = fe.system_to_component_index(i).first;
				auto N_i = N[i/dim];
				for (unsigned int j=0; j<dofs_per_cell; ++j) {	
					const unsigned int 
						component_j = fe.system_to_component_index(j).first;
					auto N_j = N[j/dim];
					cell_matrix(i,j) += N_i * C[component_i][component_j] * N_j
															* ( (i<5) ? 1:-1)
															* ( (j<5) ? 1:-1);
					AssertIsFinite(cell_matrix(i,j));
				}
			}
		}
		cell_matrix.print(std::cout,10,3);	
		return cell_matrix;												   
	};

	template<int dim, int spacedim>
	Vector<double> Material<dim,spacedim>::calc_cell_rhs(FESystem<dim,spacedim> &fe,
									typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
									QGauss<dim>  quadrature_formula,
                  Vector<double> Ue
                  )
	{
		// See calc_cell_matrix
		FEValues<dim,spacedim> fe_values(fe, quadrature_formula,
								update_values | update_gradients | 
								update_quadrature_points );
		try {
			fe_values.reinit(cell);						
		}
		catch (std::exception &exc) {
    	std::cerr << boost::diagnostic_information(exc) << std::endl;
		}
		const unsigned int n_q_points = quadrature_formula.size()/2;
		const unsigned int dofs_per_cell = fe.dofs_per_cell;
		Vector<double> cell_rhs(dofs_per_cell);
		double A, dA;
		// Tensor<1,spacedim> jump;
		// this section gets the bodyforce values 
		BodyForce<dim> body_force;
		std::vector<Vector<double> > body_force_values (n_q_points*2,
														Vector<double>(dim));
		body_force.vector_value_list (fe_values.get_quadrature_points(),
										body_force_values);


		Vector<double> qp(n_q_points),w(n_q_points), N(dofs_per_cell/dim);
		qp[0] = -1/sqrt(3);
		qp[1] = 1/sqrt(3);
		w[0] = 0.5;
		w[1] = 0.5;

		for (auto q_point = 0; q_point < n_q_points; ++q_point){

			N[0]  = 0.5*(1-qp[q_point]);
			N[1]  = 0.5*(1+qp[q_point]);
			N[2]  = 0.5*(1-qp[q_point]);
			N[3]  = 0.5*(1+qp[q_point]);

			Tensor<1,dim,double> jump;
			jump[0] = (Ue[4] - Ue[0])*N[0] + (Ue[6] - Ue[2])*N[1];
			jump[1] = (Ue[5] - Ue[1])*N[0] + (Ue[7] - Ue[3])*N[1];
			Tensor<1,dim,double> traction;
			for (unsigned int i=0; i<dim;++i) {
				traction[i] = stiffness*jump[i];
			}
			for (auto i=0; i< dofs_per_cell; ++i){
				const unsigned int 
					component_i = fe.system_to_component_index(i).first;
				auto fvol_i = N[i/dim] * body_force_values[q_point](component_i);
				auto fint_i = N[i/dim] * traction[component_i]* ( (i<5) ? 1: -1);
				// cell_rhs(i) += fint_i - fvol_i;
				cell_rhs(i) += -fint_i + fvol_i;
				AssertIsFinite(cell_rhs(i));				
			}
		}

		// assemble rhs by looping over dofs and q_points
		// for (unsigned int i=0; i<dofs_per_cell/2; ++i)
		// {
		// 	const unsigned int
		// 	component_i = fe.system_to_component_index(i).first;
    //   // jump.clear();
		// 	for (unsigned int q_point=0; q_point<n_q_points/2; ++q_point) {
    //     auto N_i = fe_values.shape_value(i,q_point);
    //     auto N_ip = fe_values.shape_value(i+ dofs_per_cell/2,q_point+n_q_points/2);
    //     auto G_i = fe_values.shape_grad(i,q_point);
    //     auto jump = N_ip * Ue[q_point+n_q_points/2] - N_i * Ue[q_point];
    //     double fint_i = stiffness * jump;
    //     double fvol_i = N_i * body_force_values[q_point](component_i);
		// 		// cell_rhs(i) += body_force_values[q_point](component_i) *
		// 		// 				fe_values.shape_value (i,q_point)
		// 		// 				*	
		// 		// 				fe_values.JxW (q_point);
    //     cell_rhs(i) += (fint_i - fvol_i) ;
    //     // cell_rhs(i) += (fint_i - fvol_i) * fe_values.JxW(q_point);
		// 		// AssertIsFinite(cell_rhs(i));				
		// 		if (!isfinite(cell_rhs(i))){
		// 			auto a = 1.;
		// 		}
		// 	}   
		// }
		return cell_rhs;
	}

}