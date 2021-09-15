/**
 * @file toplevel.h
 * @author Till Budde
 * @brief FEM Code for the masterthesis "Implementation und Analyse generalisierter Kohäsivzonenelemente in die Finite-Elemente- Bibliothek deal.II"
 * @version 0.1
 * @date 2021-06-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef TOPLEVEL_H
#define TOPLEVEL_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/base/symmetric_tensor.h>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <exception>
#include<vector>
#include "femtime.h"
#include "small-strain.h"
#include "LinElaInter.h"
#include "pointhistory.h"
#include "boundaryvalues.h"
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/numerics/data_out_dof_data.h>

#include "CustomExceptions.h"
#include <boost/throw_exception.hpp>

/**
 * @brief Namespace for the FEM Simulation
 * 
 */
namespace fem {

    using namespace dealii;

	/**
	 * @brief TopLevel implements the FEM Simultion, including meshing, assembly and solving
	 * 
	 * @tparam dim Number of Dimensions 
	 * @tparam spacedim of spatial dimensions
	 */
	template<int dim, int spacedim>
	class TopLevel
	{
		public:
		TopLevel();
		~TopLevel();
		void run ();
		
		private:
		/**
		 * @brief create / read the grid
		 * 
		 */
		void make_grid();
		
		// setup system distributes the dofs over the grid and 
		// intialises the system matrix/right-hand-side(rhs)/solution
		// with the appropriate sparsity patern	
		/**
		 * @brief distributes dofs, initialise global matrix, vectors
		 * 
		 */
		void setup_system();
		
		/**
		 * @brief call material routines and sort into global system of equations
		 * 
		 */
		void assemble_system();
		
		/**
		 * @brief Solve the global system of equations
		 * 
		 */
		void solve();
				
		void output_results();
		
		/**
		 * @brief Set the up quadrature point history object / initialise sdv's
		 * 
		 */
		void setup_quadrature_point_history();	
		
		/**
		 * @brief update state dependent variables
		 * 
		 */
		void update_quadrature_point_history();
		
		/**
		 * @brief Do the initial timestep
		 * 
		 */
		void do_initial_timestep();
		/**
		 * @brief do timestep
		 * 
		 */
		void do_timestep();
		
    // check_for_distorted_cells = true necessary for zero volume elements
    const bool check_for_distorted_cells = false;
    typename Triangulation<dim,spacedim>::MeshSmoothing smooth_mesh;
		/**
		 * @brief Contains information about the grid
		 * 
		 */
		Triangulation<dim,spacedim>   triangulation;
		
		/**
		 * @brief Distributes Dofs over grid
		 * 
		 */
		DoFHandler<dim,spacedim>	  dof_handler;
		
		// FESystem contains the basic information of the finite element 
		// system such as the type of interpolation used
		/**
		 * @brief contains information about the fe-system such as type of interpolation etc.
		 * 
		 */
		FESystem<dim,spacedim>		fe_bulk;
		FESystem<dim,spacedim>		fe_inter;
		
		/**
		 * @brief global stiffness matrix
		 * 
		 */
		SparseMatrix<double> system_matrix;
		/**
		 * @brief sparsity pattern of the stiffness matrix
		 * 
		 */
		SparsityPattern		 sparsity_pattern;
		/**
		 * @brief global displacement field
		 * 
		 */
		Vector<double>	   solution;
		Vector<double>	   solution_update;
		Vector<double>	   residual;
		/**
		 * @brief right-hand-side vector of SoE (fint,fvole,fdyn etc.)
		 * 
		 */
		Vector<double>	   system_rhs;
		
		AffineConstraints<double> constraints;
		
		/**
		 * @brief Time object implements total-time, dt number of steps etc.
		 * 
		 */
		Time time;
		/**
		 * @brief Bulk material model (Small-Strain Linear elasticity)
		 * 
		 */
		SSLinEla::Material<dim,spacedim> bulk;
		/**
		 * @brief Interface material model (linear spring)
		 * 
		 */
		LinElaInter::Material<dim,spacedim> inter;
		
		/**
		 * @brief state dependent variables at quadrature points
		 * 
		 */
		std::vector<PointHistory<dim>> quadrature_point_history;
		
		/**
		 * @brief solution / displacement field of timestep
		 * 
		 */
		Vector<double> incremental_displacement;
		
		/**
		 * @brief type of gauss quadrature used
		 * 
		 */
		QGauss<dim>  quadrature_formula_bulk;

		QGauss<dim-1>  quadrature_formula_inter;
		
		/**
		 * @brief max displacement (reached through linear loadcurve)
		 * 
		 */
		double total_displacement = 0.05;

	};


// FE Functions --------------------------------------------------------											


	template<int dim, int spacedim>
	TopLevel<dim,spacedim>::TopLevel()
		:
  triangulation(smooth_mesh, check_for_distorted_cells),
	// dof_handler handles the global numbering of degrees of freedom
	dof_handler (triangulation),
	// FE_Q<dim>(1) is a lagrangian polynomial of degree 1
	fe_bulk (FE_Q<dim,spacedim>(1), dim),
	fe_inter (FE_FaceQ<dim,spacedim>(1),dim),
	// quadrature_formula(2) is a Gauss Formula with 2 q_points in each
	// direction
	quadrature_formula_bulk(2),
	quadrature_formula_inter(2)
	{
		// Initialise the other objects here
		Time time;
		SSLinEla::Material<dim,spacedim> bulk;	
		LinElaInter::Material<dim,spacedim> inter;	
	}

	template <int dim,int spacedim>
	TopLevel<dim,spacedim>::~TopLevel()
	{
	dof_handler.clear ();
	system_matrix.clear();
	}
		
	template<int dim,int spacedim>
	void TopLevel<dim,spacedim>::make_grid()
	{
    if (dim == 2) {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream input_file("../mesh/mesh_manual-2d.msh");
      if (!input_file.is_open()){
        cexc::file_read_error exc;
        BOOST_THROW_EXCEPTION(exc);
      }
      // exception catch necessary to allow for zero volume elements 
      // also (check_for_distorted_elements = true)
      try {
        grid_in.read_msh(input_file);
      }
      catch (std::exception &exc) {
        // ignore
        std::cerr << boost::diagnostic_information(exc) << std::endl;
      }

    } else if (dim == 3){
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream input_file("../mesh/mesh_manual-3d.msh");
      if (!input_file.is_open()) {
        cexc::file_read_error exc;
        BOOST_THROW_EXCEPTION(exc);
      }
      // exception catch necessary to allow for zero volume elements 
      // also (check_for_distorted_elements = true)
      try {
        grid_in.read_msh(input_file);
      }
      catch (std::exception &exc) {
        // ignore
        // std::cerr << exc.what();
      }
    }
    // triangulation.refine_global(5);
    // try {
    // std::ofstream out("../output/grid-1.vtk");
    // if (!out.is_open()) {
    //   cexc::file_read_error exc;
    //   BOOST_THROW_EXCEPTION(exc);
    // }
    // GridOut       grid_out;
    // grid_out.write_vtk(triangulation, out);
    // std::cout << "Grid written to grid-1.svg" << std::endl;
    // }
    // catch (std::exception &exc) {
    //   std::cerr << boost::diagnostic_information(exc) << std::endl;
    // }

		// initialise the internal variables	   	
		setup_quadrature_point_history();	

		// This section names the different sides of the mesh, so that we
		// can assign different boundary condiditons
		for (typename Triangulation<dim,spacedim>::active_cell_iterator cell=triangulation.begin_active(); cell !=triangulation.end(); ++cell)
		{
			for (unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; ++f)
			{
				if (cell->face(f)->at_boundary())
				{
					const Point<dim> face_center = cell->face(f)->center();
					
					if (dim==2) {
						if (face_center[1] == -1)
							cell->face(f)->set_boundary_id (1);
						else if (face_center[1] == 3)
							cell->face(f)->set_boundary_id (2);
						else if (face_center[0] == -1)
							cell->face(f)->set_boundary_id (3);
						else
							cell->face(f)->set_boundary_id (0);			  
					} else if (dim==3) {
							if (face_center[1] == -1)
								cell->face(f)->set_boundary_id (1);
							else if (face_center[1] == 3)
								cell->face(f)->set_boundary_id (2);
							else if (face_center[0] == -1)
								cell->face(f)->set_boundary_id (3);
							else if (face_center[2] == -1)
									cell->face(f)->set_boundary_id(4);
							else
								cell->face(f)->set_boundary_id (0);			  
					}
				}
			}
		}
	}
	
	
	template <int dim,int spacedim>
	void TopLevel<dim,spacedim>::setup_system ()
	{
		dof_handler.distribute_dofs (fe_bulk);
		// use knowledge of the distribution of the dofs to create a 
		// sparsitiy pattern fo the system of equations
		DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
		DoFTools::make_sparsity_pattern (dof_handler, dsp);
		sparsity_pattern.copy_from(dsp);
		system_matrix.reinit(sparsity_pattern);
		solution.reinit(dof_handler.n_dofs());
		system_rhs.reinit(dof_handler.n_dofs());
		residual.reinit(dof_handler.n_dofs());
		incremental_displacement.reinit (dof_handler.n_dofs());
	}
	
	template<int dim, int spacedim>
	void TopLevel<dim,spacedim>::assemble_system()
	{

		system_matrix.reinit(sparsity_pattern);
		system_rhs.reinit(dof_handler.n_dofs());
		residual.reinit(dof_handler.n_dofs());

		const unsigned int  dofs_per_cell = fe_bulk.dofs_per_cell;
		FullMatrix<double>  cell_matrix(dofs_per_cell,dofs_per_cell);
		Vector<double>		cell_rhs(dofs_per_cell);
		std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
		
		// initialise the counter for the loop over cells
		typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof_handler.begin_active(),
													   endc = dof_handler.end();
		// loop over cells
		for (;cell<endc; ++cell)
		{	
			cell->get_dof_indices (local_dof_indices);
      Vector<double> Ue(dofs_per_cell); 
			for (unsigned int i=0; i<dofs_per_cell; ++i){
        Ue[i] = solution(local_dof_indices[i]);
      }
      if (cell->material_id() == 1) {
        cell_matrix = bulk.calc_cell_matrix(fe_bulk,cell,quadrature_formula_bulk,Ue); 
        cell_rhs = bulk.calc_cell_rhs(fe_bulk, cell, quadrature_formula_bulk,Ue);
      }
      else if (cell->material_id() == 2) {
        cell_rhs = inter.calc_cell_rhs(fe_inter, cell, quadrature_formula_inter,Ue);
        cell_matrix = inter.calc_cell_matrix(fe_inter,cell,quadrature_formula_inter,Ue); 
      }
      else {
        cexc::not_mat_error exc;
        BOOST_THROW_EXCEPTION(exc);
      }

			// this section arranges the cell matrix/rhs into the global
			// system of equations
			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{
					system_matrix.add (local_dof_indices[i],
										local_dof_indices[j],
										cell_matrix(i,j));
				}			
				system_rhs(local_dof_indices[i]) += cell_rhs(i);
			}
		} 
		// Here the vertical displacement of face(1) (the lower end)
		// is clamped
		std::vector<bool> xmask,ymask,zmask;
		if (dim==2) {
			ymask = std::vector<bool>{false,true}, xmask=std::vector<bool>{true,false};
		} else if (dim==3) {
			xmask=std::vector<bool>{true,false,false}, ymask=std::vector<bool>{false,true,false}, zmask=std::vector<bool>{false,false,true};
		}
			std::map<types::global_dof_index,double> boundary_values;

			constraints.clear();
			VectorTools::interpolate_boundary_values (
					dof_handler,
					1,
					Functions::ZeroFunction<dim>(dim),
					constraints,
					ComponentMask(ymask)
					);
			// Here we assigne the BC's calculated in BoundaryValues to
			// face(2) (the upper side)
			VectorTools::interpolate_boundary_values (
					dof_handler,
					2,
					BoundaryValues<dim>(time.get_timestep(), time.get_no_timesteps(), total_displacement),
					constraints,
					ComponentMask(ymask)
					); 

			VectorTools::interpolate_boundary_values (
					dof_handler,
					3,
					Functions::ZeroFunction<dim>(dim),
					constraints,
					ComponentMask(xmask)
					);

			if (dim == 3) {
				VectorTools::interpolate_boundary_values (
						dof_handler,
						4,
						Functions::ZeroFunction<dim>(dim),
						constraints,
						ComponentMask(zmask)
						);
			}

			constraints.close();

		// No boundary values are set for face(0) therefore a zero
		// force / stress BC is applied
											
		// Save the incremental displacement for the update of 
		// quadrature point history									
		incremental_displacement = solution;
	}
	
	
	template<int dim, int spacedim>
	void TopLevel<dim,spacedim>::solve()
	{
		SparseDirectUMFPACK solver;
		solution_update.reinit(dof_handler.n_dofs());
		solver.initialize(system_matrix);
		// solver.factorize(system_matrix);
		solver.vmult(solution_update,residual);
		solution -= solution_update;
		constraints.distribute(solution);
	}
	
	template <int dim, int spacedim>
	void TopLevel<dim,spacedim>::output_results()
	{
		// This section creates an output
		DataOut<dim> data_out;
		// data_out.set_flags(write_higher_order_elements=false);
		data_out.attach_dof_handler (dof_handler);
		std::vector<DataComponentInterpretation::DataComponentInterpretation> component_type;
		if (dim==2) {
			component_type = {DataComponentInterpretation::component_is_part_of_vector, DataComponentInterpretation::component_is_part_of_vector};
		} else {
			component_type = {DataComponentInterpretation::component_is_part_of_vector, DataComponentInterpretation::component_is_part_of_vector, DataComponentInterpretation::component_is_part_of_vector};
		}
		data_out.add_data_vector(dof_handler, solution, "DSPL", component_type);
		data_out.add_data_vector(dof_handler, system_rhs, "RHS", component_type);
		data_out.build_patches();
		
		// This section determines the format of the output
		std::string filename = "../output/solution_";
		std::string timestep = std::to_string(time.get_timestep());
		int prec_zeros = 4 - timestep.length();
		for (auto i=0 ; i<prec_zeros; ++i){
			filename += "0";
		}
		filename += timestep + ".vtk";
		std::ofstream output_vtk (filename);
		if (!output_vtk.is_open()) {
        cexc::file_error exc;
        BOOST_THROW_EXCEPTION(exc);
		}
			
		// std::ofstream output_grid ("../output/grid.svg");
    // DataOut::VtkFlags
		data_out.write_vtk (output_vtk);
		// grid_out.write_svg(triangulation,output_grid);

	}

	template<int dim, int spacedim>
	void TopLevel<dim, spacedim>::setup_quadrature_point_history()
	{
		// This section makes sure that quadrature_point_history is 
		// empty and of the correct size
		unsigned int size = 0;
		for (auto	cell : dof_handler.active_cell_iterators())
			{
				triangulation.clear_user_data();
				if (cell->material_id()==1) {
					size+= quadrature_formula_bulk.size();
				} else if (cell->material_id() == 2) {
					size+= quadrature_formula_inter.size();
				} else {
					cexc::not_imp_error exc;
					BOOST_THROW_EXCEPTION(exc);
				}
			}
		std::vector<PointHistory<dim> > tmp;
		tmp.swap (quadrature_point_history);

		quadrature_point_history.resize(size);
		
		// This section the quadratrue_point_history to the quadrature 
		// points of each cell	
		unsigned int history_index = 0;
		for (auto cell : dof_handler.active_cell_iterators())
		{			
			cell->set_user_pointer (&quadrature_point_history[history_index]);
			if (cell->material_id()==1) {
				history_index += quadrature_formula_bulk.size();
			} else if (cell->material_id() == 2) {
				history_index += quadrature_formula_inter.size();
			} else {
				cexc::not_imp_error exc;
				BOOST_THROW_EXCEPTION(exc);
			}
		}
		Assert (history_index == quadrature_point_history.size(),ExcInternalError());			
	}

	// This section demonstrates how to update the
	// quadrature_point_history after the solution of the global system
	// has been obtained  it has no real use in the current form of this
	// code
	template<int dim, int spacedim>
	void TopLevel<dim,spacedim>::update_quadrature_point_history()
	{
		
		FEValues<dim,spacedim> fe_values (fe_bulk, quadrature_formula_bulk,
								update_values | update_gradients);
		std::vector<std::vector<Tensor<1,dim>>>
		displacement_increment_grads (quadrature_formula_bulk.size(), 
									std::vector<Tensor<1,dim>>(dim));
		// loop over cells							
		for (auto cell : dof_handler.active_cell_iterators())
		{
			// get the local_quadrature_points_history of each cell
			PointHistory<dim> *local_quadrature_points_history
				= reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());	
			Assert (local_quadrature_points_history >= &quadrature_point_history.front(),
					ExcInternalError());
			Assert (local_quadrature_points_history < &quadrature_point_history.back(),
					ExcInternalError());
			
			// this section calculates the updated internal variables		
			// fe_values.reinit(cell);
			// fe_values.get_function_gradients (incremental_displacement,
			// 								  displacement_increment_grads);
			if (cell->material_id()==1) {
				for (unsigned int q=0; q<quadrature_formula_bulk.size(); ++q) {
					unsigned int dofs_per_cell = fe_bulk.dofs_per_cell;
					std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
					cell->get_dof_indices (local_dof_indices);
					Vector<double> Ue(dofs_per_cell); 
					for (unsigned int i=0; i<dofs_per_cell; ++i){
						Ue[i] = solution(local_dof_indices[i]);
					}
					SymmetricTensor<2,dim> stress = bulk.calc_stress(fe_bulk,cell,quadrature_formula_bulk, Ue, q);
					local_quadrature_points_history[q].old_stress = stress;
				}
			} else if (cell->material_id() == 2) {
				for (unsigned int q=0; q<quadrature_formula_inter.size(); ++q) {
				// local_quadrature_points_history[q].old_stress = new_stress;
				}
			} else {
				cexc::not_imp_error exc;
				BOOST_THROW_EXCEPTION(exc);
			}
		}
	}

	template<int dim,int spacedim>
	void TopLevel<dim,spacedim>::do_initial_timestep()
	{
		std::cout << "doing initial timestep" << std::endl;
		make_grid();
		std::cout << "grid made" << std::endl;
		setup_system();
		// assemble_system();
		// solve();
		std::cout << "initial step completed" << std::endl;
		output_results(); 
	}
	
	template<int dim,int spacedim>
	void TopLevel<dim,spacedim>::do_timestep()
	{
		time.increment();
		std::cout << "Timestep No. " << time.get_timestep() << " time " << time.get_current() << std::endl;
		double rsn = 1.;
		unsigned int iter = 0;
		assemble_system();
		constraints.condense(system_matrix,system_rhs);
		system_matrix.vmult(residual,solution);
		residual.add(-1,system_rhs);
		while (rsn > 1e-12) {
			solve();
			iter++;
			assemble_system();
			system_matrix.vmult(residual,solution);
			residual.add(-1,system_rhs);
			rsn = residual.l2_norm();
			std::cout << "  iter = " << iter << ",  residual = " << rsn << std::endl;
			constraints.condense(system_matrix,system_rhs);
			if (iter > 15) {
        cexc::convergence_error exc;
        BOOST_THROW_EXCEPTION(exc);
			}
		}
		std::cout << "Step completed" << std::endl;
		update_quadrature_point_history();
		output_results(); 
	}
	   
	template<int dim, int spacedim>
	void TopLevel<dim,spacedim>::run()
	{
		 do_initial_timestep();
		 while (time.get_current() < time.get_end())
		 {
			do_timestep();
		 }			
		//  output_results(); 
		 std::cout << "output results" << std::endl;
	 }
	 
// FE Functions End ----------------------------------------------------

}

#endif