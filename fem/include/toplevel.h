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
#include <deal.II/lac/precondition.h>
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
#include <deal.II/base/symmetric_tensor.h>
#include <fstream>
#include <iostream>
#include <assert.h>
#include "femtime.h"
#include "small-strain.h"
#include "pointhistory.h"
#include "boundaryvalues.h"
#include <deal.II/fe/fe_interface_values.h>

namespace fem {

    using namespace dealii;

	// The TopLevel controlls the fem simulation	
	template<int dim, int spacedim>
	class TopLevel
	{
		public:
		TopLevel();
		~TopLevel();
		void run ();
		
		private:
		// Make grid creates the grid		
		void make_grid();
		
		// setup system distributes the dofs over the grid and 
		// intialises the system matrix/right-hand-side(rhs)/solution
		// with the appropriate sparsity patern	
		void setup_system();
		
		// assemble system calls material.calc_cell_matrix and 
		// material.calc_cell_rhs and arranges them into the system 
		// matrix/rhs it also calls the boundary function		
		void assemble_system();
		
		//	solves the system build by assemble system	
		void solve();
				
		void output_results();
		
		// setup_quadrature_point _history sets up the internal 
		// variables at each quadrature point - in this programms all
		// are initialized as zero		
		void setup_quadrature_point_history();	
		
		// updates quadrature points after the timestep has been solved 		
		void update_quadrature_point_history();
		
		void do_initial_timestep();
		void do_timestep();
		
		// Triangulation contains information of the grid
		// such as the number of active cells and their iterators
		Triangulation<dim,spacedim>   triangulation;
		
		// The DoFHandler distributes the dofs over the grid
		DoFHandler<dim,spacedim>	  dof_handler;
		
		// FESystem contains the basic information of the finite element 
		// system such as the type of interpolation used
		FESystem<dim,spacedim>		fe;
		
		SparseMatrix<double> system_matrix;
		SparsityPattern		 sparsity_pattern;
		Vector<double>	   solution;
		Vector<double>	   system_rhs;
		
		Time time;
		Material<dim,spacedim> material;
		
		// A Vector of type PointHistory to save the internal variables 
		// at each quadrature point
		std::vector<PointHistory<dim>> quadrature_point_history;
		
		// incremental displacement saves the solution of the 
		// displacement field during each step for further use in
		// update_quadrature_point_history
		Vector<double> incremental_displacement;
		
		// QGauss is a class for gauss type integration
		QGauss<dim>  quadrature_formula;
		
		// Could be outsourced to a user input function
		double total_displacement = 0.05;

	};


// FE Functions --------------------------------------------------------											


	template<int dim, int spacedim>
	TopLevel<dim,spacedim>::TopLevel()
		:
	// dof_handler handles the global numbering of degrees of freedom
	dof_handler (triangulation),
	// FE_Q<dim>(1) is a lagrangian polynomial of degree 1
	fe (FE_Q<dim,spacedim>(2), dim),
	// quadrature_formula(2) is a Gauss Formula with 2 q_points in each
	// direction
	quadrature_formula(3)
	{
		// Initialise the other objects here
		Time time;
		Material<dim,spacedim> material;	
	}

	template <int dim,int spacedim>
	TopLevel<dim,spacedim>::~TopLevel ()
	{
	dof_handler.clear ();
	system_matrix.clear();
	}
		
	template<int dim,int spacedim>
	void TopLevel<dim,spacedim>::make_grid()
	{
    assert(dim==2);
		// create a cube or quadrat depending on dim with edge 
		// coordinates 1, -1 	
		// GridGenerator::hyper_cube (triangulation, -1, 1);
		// refine the mesh
		// triangulation.refine_global (1);

    // create two-element mesh

    if (dim == 2) {
      const std::vector<Point<dim>> vertices = {
        {-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {1.0, 3.0}, {-1.0, 3.0}, {-1.0, 1.0}
      };
      const std::vector<std::array<int, GeometryInfo<dim>::vertices_per_cell>>
        cell_vertices = {
          {0,1,5,2},
          {5,2,4,3}
        };
      const unsigned int n_cells = cell_vertices.size();
      std::vector<CellData<dim>> cells(n_cells , CellData<dim>());
      for (unsigned int itercell = 0; itercell < n_cells; ++itercell) {
        for (unsigned int itervertice = 0; itervertice < cell_vertices[itercell].size(); ++itervertice) {
          cells[itercell].vertices[itervertice] = cell_vertices[itercell][itervertice];
        }
        cells[itercell].material_id = 0;
      }  
      triangulation.create_triangulation(vertices,cells,SubCellData());
    } else if (dim == 3){
    }

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
					
					if (face_center[dim-1] == -1)
						cell->face(f)->set_boundary_id (1);
					else if (face_center[dim-1] == 3)
						cell->face(f)->set_boundary_id (2);
					else
						cell->face(f)->set_boundary_id (0);			  
				}
			}
		}
	}
	
	
	template <int dim,int spacedim>
	void TopLevel<dim,spacedim>::setup_system ()
	{
		dof_handler.distribute_dofs (fe);
		// use knowledge of the distribution of the dofs to create a 
		// sparsitiy pattern fo the system of equations
		DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
		DoFTools::make_sparsity_pattern (dof_handler, dsp);
		sparsity_pattern.copy_from(dsp);
		system_matrix.reinit(sparsity_pattern);
		solution.reinit(dof_handler.n_dofs());
		system_rhs.reinit(dof_handler.n_dofs());
		incremental_displacement.reinit (dof_handler.n_dofs());
	}
	
	template<int dim, int spacedim>
	void TopLevel<dim,spacedim>::assemble_system()
	{
		const unsigned int  dofs_per_cell = fe.dofs_per_cell;
		FullMatrix<double>  cell_matrix(dofs_per_cell,dofs_per_cell);
		Vector<double>		cell_rhs(dofs_per_cell);
		std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
		
		// initialise the counter for the loop over cells
		typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof_handler.begin_active(),
													   endc = dof_handler.end();
		// loop over cells
		for (;cell<endc; ++cell)
		{	
			cell_matrix = material.calc_cell_matrix(fe,cell,quadrature_formula); 
			cell_rhs = material.calc_cell_rhs(fe, cell, quadrature_formula);
			
			// this section arranges the cell matrix/rhs into the global
			// system of equations
			cell->get_dof_indices (local_dof_indices);
			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{
					system_matrix.add (local_dof_indices[i],
										local_dof_indices[j],
										cell_matrix(i,j));
					system_rhs(local_dof_indices[i]) += cell_rhs(i);
				}			
			}
		} 
		// Here the vertical displacement of face(1) (the lower end)
		// is clamped
		std::map<types::global_dof_index,double> boundary_values;
		VectorTools::interpolate_boundary_values (dof_handler,
        1,
        Functions::ZeroFunction<dim>(dim),
        boundary_values);
												  
		// Here we assigne the BC's calculated in BoundaryValues to
		// face(2) (the upper side)
		VectorTools::interpolate_boundary_values (dof_handler,
        2,
        BoundaryValues<dim>(time.get_timestep(), time.get_no_timesteps(), total_displacement),
        boundary_values); 
		// No boundary values are set for face(0) therefore a zero
		// force / stress BC is applied
		
		// apply the boundary values here
		MatrixTools::apply_boundary_values (boundary_values,
        system_matrix,
        solution,
        system_rhs);
											
		// Save the incremental displacement for the update of 
		// quadrature point history									
		incremental_displacement = solution;
	}
	
	
	template<int dim, int spacedim>
	void TopLevel<dim,spacedim>::solve()
	{
		// Max number of iterations is 1000, the tolerance is 1e-12
		SolverControl		   solver_control (1000, 1e-12);
		// Use a conjugated gradient solver
		SolverCG<>			  cg (solver_control);
		// No preconditioner is used here
		cg.solve (system_matrix, solution, system_rhs,
					PreconditionIdentity());	
	}
	
	template <int dim, int spacedim>
	void TopLevel<dim,spacedim>::output_results()
	{
		// This section creates an output
		DataOut<dim> data_out;
		GridOut grid_out;
		data_out.attach_dof_handler (dof_handler);
		std::vector<std::string> solution_names;
		solution_names.emplace_back("x_displacement");
		solution_names.emplace_back("y_displacement");
		data_out.add_data_vector (solution, solution_names);
		data_out.build_patches ();
		
		// This section determines the format of the output
		std::ofstream output_vtk ("output/solution.vtk");
		std::ofstream output_grid ("output/grid.svg");
		data_out.write_vtk (output_vtk);
		grid_out.write_svg(triangulation,output_grid);

	}

	template<int dim, int spacedim>
	void TopLevel<dim, spacedim>::setup_quadrature_point_history()
	{
		// This section makes sure that quadrature_point_history is 
		// empty and of the correct size
		unsigned int ncell = 0;
		for (typename Triangulation<dim>::active_cell_iterator
			cell = triangulation.begin_active();
			cell != triangulation.end(); ++cell)
			{
				triangulation.clear_user_data();
				ncell++;
			}
		std::vector<PointHistory<dim> > tmp;
		tmp.swap (quadrature_point_history);

		quadrature_point_history.resize(ncell*quadrature_formula.size());
		
		// This section the quadratrue_point_history to the quadrature 
		// points of each cell	
		unsigned int history_index = 0;
		for (typename Triangulation<dim,spacedim>::active_cell_iterator
			cell = triangulation.begin_active();
			cell != triangulation.end(); ++cell)
		{			
			cell->set_user_pointer (&quadrature_point_history[history_index]);
			history_index += quadrature_formula.size();
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
		
		FEValues<dim,spacedim> fe_values (fe, quadrature_formula,
								update_values | update_gradients);
		std::vector<std::vector<Tensor<1,dim>>>
		displacement_increment_grads (quadrature_formula.size(), 
									std::vector<Tensor<1,dim>>(dim));
		// loop over cells							
		for (typename DoFHandler<dim,spacedim>::active_cell_iterator
			cell = dof_handler.begin_active(); cell != dof_handler.end();
			++cell)
		{
			// get the local_quadrature_points_history of each cell
			PointHistory<dim> *local_quadrature_points_history
				= reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());	
			Assert (local_quadrature_points_history >= &quadrature_point_history.front(),
					ExcInternalError());
			Assert (local_quadrature_points_history < &quadrature_point_history.back(),
					ExcInternalError());
			
			// this section calculates the updated internal variables		
			fe_values.reinit(cell);
			fe_values.get_function_gradients (incremental_displacement,
											  displacement_increment_grads);
			
			for (unsigned int q=0; q<quadrature_formula.size(); ++q)
			{
				const SymmetricTensor<2,dim> 
				new_stress = material.calc_stress( (get_strain(displacement_increment_grads[q]) - local_quadrature_points_history[q].strain_pl) );
				local_quadrature_points_history[q].old_stress = new_stress;
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
		assemble_system();
		solve();
		std::cout << "initial step completed" << std::endl;
	}
	
	template<int dim,int spacedim>
	void TopLevel<dim,spacedim>::do_timestep()
	{
		time.increment();
		std::cout << "Timestep No. " << time.get_timestep() << " time " << time.get_current() << std::endl;
		assemble_system();
		solve();
		std::cout << "Step completed" << std::endl;
		update_quadrature_point_history();
	}
	   
	template<int dim, int spacedim>
	void TopLevel<dim,spacedim>::run()
	{
		 do_initial_timestep();
		 while (time.get_current() < time.get_end())
		 {
			do_timestep();
		 }			
		 output_results(); 
		 std::cout << "output results" << std::endl;
	 }
	 
// FE Functions End ----------------------------------------------------

}

#endif