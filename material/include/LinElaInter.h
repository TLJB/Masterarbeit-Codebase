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
#include <deal.II/grid/tria.h>
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
							  QGauss<dim-1>  quadrature_formula,
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
									QGauss<dim-1>  quadrature_formula,
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
									QGauss<dim-1>  quadrature_formula,
                  Vector<double> Ue
                  )
	{
		// the constructor of fe_values needs fe & quadrature_formula,
		// the other arguments determine which values are calculated by 
		// fe_values(reinit)		
		FEFaceValues<dim,spacedim> fe_values(fe, quadrature_formula,
								update_values | update_JxW_values | 
								update_quadrature_points);
		const unsigned int n_q_points = quadrature_formula.size();
		const unsigned int n_faces = cell->n_faces();
		const unsigned int dofs_per_cell = fe.dofs_per_cell/n_faces*2;
		FullMatrix<double> cell_matrix(dofs_per_cell);		

		typename Triangulation<dim,spacedim>::active_cell_iterator t=cell;

    auto identity =  Physics::Elasticity::StandardTensors<dim>::I;
		auto C = stiffness*identity;

		for (unsigned int q_point = 0; q_point != n_q_points; ++q_point) {

			for (unsigned int face_i : {n_faces-2,n_faces-1}) {
				unsigned int dofs_per_face = fe.n_dofs_per_face(face_i,0);	
				for (unsigned int face_dof_i=0; face_dof_i != dofs_per_face; ++face_dof_i) {
					fe_values.reinit(t,face_i);	
					auto component_i = fe.face_system_to_component_index(face_dof_i,face_i).first;
					auto node_i = fe.face_system_to_component_index(face_dof_i,face_i).second;
					auto index_i = dofs_per_face*(face_i - (n_faces-2)) + node_i*fe.n_components() + component_i;
					auto cell_index = fe.face_to_cell_index(
												face_dof_i,face_i,cell->face_orientation(face_i),false,false);
					double N_i = fe_values.shape_value(
											fe.face_to_cell_index(
												face_dof_i,face_i,cell->face_orientation(face_i),false,false),q_point);

					for (unsigned int face_j : {n_faces-2,n_faces-1}) {
						fe_values.reinit(t,face_j);	
						for (unsigned int face_dof_j=0; face_dof_j != dofs_per_face; ++face_dof_j) {
							auto component_j = fe.face_system_to_component_index(face_dof_j,face_j).first;
							auto node_j = fe.face_system_to_component_index(face_dof_j,face_j).second;
							auto index_j = dofs_per_face*(face_j - (n_faces-2)) + node_j*fe.n_components() + component_j;
							double N_j = fe_values.shape_value(
													fe.face_to_cell_index(
														face_dof_j,face_j,cell->face_orientation(face_j),false,false),q_point);
							
							cell_matrix(index_i,index_j) +=
								N_i * C[component_i][component_j] * N_j
									* fe_values.JxW(q_point)
									* ( (face_i == n_faces -2) ? 1 : -1)
									* ( (face_j == n_faces -2) ? 1 : -1);
							AssertIsFinite(cell_matrix(index_i,index_j));
						}
					}
				}
			}
		}

		return cell_matrix;												   
	};

	template<int dim, int spacedim>
	Vector<double> Material<dim,spacedim>::calc_cell_rhs(FESystem<dim,spacedim> &fe,
									typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
									QGauss<dim-1>  quadrature_formula,
                  Vector<double> Ue
                  )
	{
		// See calc_cell_matrix
		FEFaceValues<dim,spacedim> fe_values(fe, quadrature_formula,
								update_values | update_JxW_values | 
								update_quadrature_points);

		const unsigned int n_q_points = quadrature_formula.size();
		const unsigned int n_faces = cell->n_faces();
		const unsigned int dofs_per_cell = fe.dofs_per_cell/n_faces*2;
		Vector<double> cell_rhs(dofs_per_cell);


		typename Triangulation<dim,spacedim>::active_cell_iterator t=cell;

		for (auto q_point = 0; q_point < n_q_points; ++q_point){

			Tensor<1,dim,double> jump;
			for (unsigned int face : {n_faces-2,n_faces-1}) {
				fe_values.reinit(t,face);	
				unsigned int dofs_per_face = fe.n_dofs_per_face(face,0);	
				for (unsigned int face_dof=0; face_dof != dofs_per_face; ++face_dof) {
					auto component_i = fe.face_system_to_component_index(face_dof,face).first;
					auto node_i = fe.face_system_to_component_index(face_dof,face).second;
					auto vector_index = dofs_per_face*(face - (n_faces-2)) + node_i*fe.n_components() + component_i;
					jump[component_i] += Ue[vector_index]
															* fe_values.shape_value(
																fe.face_to_cell_index(
																	face_dof,face,cell->face_orientation(face),false,false),q_point)
															* ( (face == n_faces -2) ? -1 : 1);
				}
			}

    	auto identity =  Physics::Elasticity::StandardTensors<dim>::I;
			auto C=stiffness*identity;
			Tensor<1,dim,double> traction=C*jump;

			for (unsigned int face : {n_faces-2,n_faces-1}) {
				typename Triangulation<dim,spacedim>::active_cell_iterator t=cell;
				fe_values.reinit(t,face);	
				unsigned int dofs_per_face = fe.n_dofs_per_face(face,0);	
				for (unsigned int face_dof=0; face_dof != dofs_per_face; ++face_dof) {
					auto component_i = fe.face_system_to_component_index(face_dof,face).first;
					auto node_i = fe.face_system_to_component_index(face_dof,face).second;
					auto vector_index = dofs_per_face*(face - (n_faces-2)) + node_i*fe.n_components() + component_i;

					auto fint_i = fe_values.shape_value(
													fe.face_to_cell_index(
														face_dof,face,cell->face_orientation(face),false,false),q_point)
												* traction[component_i]
												* fe_values.JxW(q_point)
												* ( ( face == n_faces -2) ? -1 : +1);

					cell_rhs(vector_index) += fint_i;
					AssertIsFinite(cell_rhs(vector_index));
				}
			}
		}

		return cell_rhs;
	}

}