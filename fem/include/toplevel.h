/**
 * @file toplevel.h
 * @author Till Budde
 * @brief FEM Code for the masterthesis "Implementation und Analyse
 * generalisierter Kohäsivzonenelemente in die Finite-Elemente- Bibliothek
 * deal.II"
 * @version 0.1
 * @date 2021-06-28
 *
 * @copyright Copyright (c) 2021
 *
 */
#ifndef TOPLEVEL_H
#define TOPLEVEL_H

#include "DamInter.h"
#include "GenElaInter.h"
#include "LinElaInter.h"
#include "NeoHooke.h"
#include "boundaryvalues.h"
#include "femtime.h"
#include "pointhistory.h"
#include "small-strain.h"
#include <assert.h>
#include <cstdlib>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_dof_data.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string.h>
#include <tuple>
#include <vector>

#include "CustomExceptions.h"
#include <boost/throw_exception.hpp>

#include "timer.h"

/**
 * @brief Main namespace for the FEM calculation
 *
 * Within this namespace all classes, functions etc. are included that
 * are not specific to the different materials
 *
 */
namespace fem {

using namespace dealii;

/**
 * @brief TopLevel implements the FEM Simultion, including meshing, assembly and
 * solving
 *
 * @tparam dim Number of Dimensions
 * @tparam spacedim of spatial dimensions
 */
template <int dim, int spacedim, bool num = false> class TopLevel {
public:
  TopLevel();
  ~TopLevel();
  void run();

private:
  /**
   * @brief create / read the grid
   *
   */
  void make_grid();

  /**
   * @brief distributes dofs, initialise global matrix, vectors
   *
   * This function read the gmsh mesh and attaches it to the triangulation.
   * Afterwards the faces of the mesh are assigned an ID, that is used for
   * the application of Boundary Conditions
   * This Function also calls the setup_quadrature_point_history() function.
   *
   */
  void setup_system();

  /**
   * @brief call material routines and sort into global system of equations
   *
   * This function creates the linear system of equations that needs to be
   * solved (in this newton step).
   * To do this it first loops over all cells and calls the appropriate
   * material function (inter or bulk).
   * Aftterward the cell_matrix / cell_rhs is distributed to the system matrix
   * / rhs.
   *
   * After the system of equations is build, the boundary conditions are applied
   */
  void assemble_system();

  /**
   * @brief Create a constraints object
   *
   * This function clears the \ref constraints object and then
   * creates the new constraints for this timestep.
   *
   */
  void create_constraints();

  /**
   * @brief Solve the linear system of equations
   *
   * This function solves the linear system of equations of the Newton-Raphson
   * step.
   * It uses a direct solver to do so (SparseDirectUMFPACK).
   * afterwards the boundary conditions are distributed to the solution
   * (the displacement field)
   *
   */
  void solve(unsigned int iter);

  /**
   * @brief Writes the output to a vtk file
   *
   * This function writes the displacement field and RHS-Vector to a
   * vtk file
   *
   */
  void output_results();

  /**
   * @brief Set the up quadrature point history object / initialise sdv's
   *
   * Allocate Memory for the storage of interanal state dependent variables
   * in the quadrature_point_history vector according to the type of each
   * element.
   * All sdv's are initialized as 0.
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
   * Calls the following functions in order
   * - make_grid()	-> creates the grid
   * - setup_system() -> allocates memory, distributes dofs
   * - output_results() ->  writes output vtk files (initial state)
   *
   */
  void do_initial_timestep();

  /**
   * @brief Solves the nonlinear system of equations
   *
   * Increments the time and then uses a newton-raphson scheme to solve
   * the nonlinear system of equations.
   *
   * Calls the functions
   * - assemble_system() -> Assembles system matrix / RHS and creates BC'x
   * - solve() -> Solves linear system and applies BC's to displacement field
   *
   */
  void do_timestep();

  /**
   * @brief Check whether the material is of type bulk or interface
   *
   * @return int (1) if the material is of type bulk, (2) if of type interface
   * (-1) if neither
   */
  int check_material_id(int);

  /**
   * @brief check_for_distorted_cells = true necessary for zero volume elements
   * or was in version 9.2 of deal
   *
   */
  const bool check_for_distorted_cells = true;
  /**
   * @brief type of mesh smoothing
   *
   * None is required so because of zero-volume interface elements
   */
  typename Triangulation<dim, spacedim>::MeshSmoothing smooth_mesh =
      Triangulation<dim, spacedim>::MeshSmoothing::none;
  /**
   * @brief Contains information about the grid
   * (https://www.dealii.org/current/doxygen/deal.II/classTriangulation.html)
   *
   * This class stores information about the mesh, but does not know anything
   * about the degrees of freedom.
   * Also supplies an iterator over the cells of the mesh, much like
   * the \ref dof_handler.
   *
   */
  Triangulation<dim, spacedim> triangulation;

  /**
   * @brief Distributes Dofs over grid
   * (https://www.dealii.org/current/doxygen/deal.II/classDoFHandler.html)
   *
   * Using the \ref triangulation and the finite-element \ref fe_bulk this class
   * first distributes the degrees of freedom.
   * Afterwards the dof_handler supplies iterators over the degrees of freedom
   * and cells.
   *
   */
  DoFHandler<dim, spacedim> dof_handler;

  /**
   * @brief contains information about the fe-system of the bulk material
   * 				such as type of interpolation etc.
   *
   */
  FESystem<dim, spacedim> fe_bulk;
  /**
   * @brief contains information about the fe-system of the interface
   * 		such as type of interpolation etc.
   *
   */
  FESystem<dim, spacedim> fe_inter;

  /**
   * @brief global stiffness matrix
   *
   */
  SparseMatrix<double> system_matrix;
  /**
   * @brief sparsity pattern of the stiffness matrix
   *
   */
  SparsityPattern sparsity_pattern;
  /**
   * @brief global displacement field
   *
   */
  Vector<double> solution;
  /**
   * @brief Newton Update of the global displacement field
   *
   */
  Vector<double> solution_update;
  /**
   * @brief Residual of the Newton - Rapshon scheme
   *
   */
  Vector<double> residual;
  /**
   * @brief right-hand-side vector of SoE (fint,fvole,fdyn etc.)
   *
   */
  Vector<double> system_rhs;

  /**
   * @brief Implementation of BC's
   * (https://www.dealii.org/current/doxygen/deal.II/classAffineConstraints.html)
   *
   * Constrains the degrees of freedom of the \ref system_matrix / \ref
   * system_rhs by setting the matrix and vector entries to the proper values.
   * After the linear system is solved distributes the proper (Dirichlet)
   * Boundary Values to the \ref solution vector.
   *
   */
  AffineConstraints<double> constraints;

  /**
   * @brief Time object implements total-time, dt number of steps etc.
   *
   */
  Time time;
  /**
   * @brief Bulk material model (NeoHooke)
   *
   */
  NeoHooke::Material<dim, spacedim> bulk;
  /**
   * @brief Interface material model (Change Namespace accordingly)
   *
   */
  LinElaInter::Material<dim, spacedim> inter;

  /**
   * @brief state dependent variables at quadrature points
   *
   */
  std::vector<PointHistory<dim>> quadrature_point_history;

  /**
   * @brief type of gauss quadrature used for the bulk material
   *
   */
  QGauss<dim> quadrature_formula_bulk;

  /**
   * @brief type of gauss quadrature used for the interface material
   *
   */
  QGauss<dim - 1> quadrature_formula_inter;

  /**
   * @brief max displacement (reached through linear loadcurve)
   *
   */
  double total_displacement = 0.2;
  /**
   * @brief Damping parameter of Newton Raphson scheme
   *
   */
  double damping_parameter = 1;

  /**
   * @brief Vector of material-ids associated with bulk material
   *
   */
  std::vector<int> bulk_ids;
  /**
   * @brief Vector of material-ids associated with interface material
   *
   */
  std::vector<int> interface_ids;

  /**
   * @brief Filename of the .Msh-file
   *
   */
  std::string mesh_filename;

  /**
   * @brief Pertubation used for numerical tangent (if num)
   *
   */
  double pertubation_tangent = 1e-6;
  /**
   * @brief tolerance of residual of newton raphson scheme
   *
   */
  double tol_newton = 1e-8;
  /**
   * @brief Maximum number of Newton steps
   *
   */
  unsigned int max_iter_newton = 50;
};

// FE Functions --------------------------------------------------------

template <int dim, int spacedim, bool num>
TopLevel<dim, spacedim, num>::TopLevel()
    : // check_for_distorted_cells=false necessary, so that the initial mesh
      // with its zero volume interface elements can be read / stored.
      triangulation(smooth_mesh, check_for_distorted_cells),
      // link the dof_handler to the triangulation.
      // Allows for the distribution of dof's over the mesh stored in
      // triangulation
      dof_handler(triangulation),
      // FE_Q<dim>(1) is a lagrangian polynomial of degree 1
      // uses the basic (1D) lagrange polynomials to assemble a dim-dimensional
      // Finite Element
      fe_bulk(FE_Q<dim, spacedim>(1), dim),
      // See above but for Finite Elements that operate on faces.
      // Note: FE_FaceQ does not number dofs / behave in an expected way. Refer
      // to official documentation.
      fe_inter(FE_FaceQ<dim, spacedim>(1), dim),
      // quadrature_formula(2) is a Gauss Formula with 2 q_points in each
      // direction
      quadrature_formula_bulk(2), quadrature_formula_inter(2) {
  // Initialise the other objects here
  if (dim == 2) {
    mesh_filename = "../mesh/2d.msh";
  } else if (dim == 3) {
    mesh_filename = "../mesh/3d.msh";
  }
}

template <int dim, int spacedim, bool num>
TopLevel<dim, spacedim, num>::~TopLevel() {
  dof_handler.clear();
  system_matrix.clear();
}

template <int dim, int spacedim, bool num>
int TopLevel<dim, spacedim, num>::check_material_id(int id) {

  // if both vectors are empty probably a Version 1 .Msh file was used
  // test mesh used ids 2 and 1 so fill the vector with those
  if (bulk_ids.empty() && interface_ids.empty()) {
    interface_ids = {2};
    bulk_ids = {1};
  }
  if (std::find(bulk_ids.begin(), bulk_ids.end(), id) != bulk_ids.end()) {
    return 1;
  } else if (std::find(interface_ids.begin(), interface_ids.end(), id) !=
             interface_ids.end()) {
    return 2;
  } else {
    return -1;
  }
}

template <int dim, int spacedim, bool num>
void TopLevel<dim, spacedim, num>::make_grid() {

  std::ifstream input_file(mesh_filename);
  if (!input_file.is_open()) {
    cexc::file_read_error exc;
    BOOST_THROW_EXCEPTION(exc);
  }
  // Object that reads the input mesh file
  GridIn<dim> grid_in;
  // link the input reader to the triangulation
  grid_in.attach_triangulation(triangulation);
  // exception catch necessary to allow for zero volume elements
  try {
    grid_in.read_msh(input_file);
  } catch (std::exception &exc) {
    // ignore
    // std::cerr << boost::diagnostic_information(exc) << std::endl;
  }

  // Read which material ids correspond to bulk and interface
  input_file.close();
  input_file.open(mesh_filename);
  if (!input_file.is_open()) {
    cexc::file_read_error exc;
    BOOST_THROW_EXCEPTION(exc);
  }
  bool stop = false;
  std::string line;
  while (input_file.good() && !stop) {
    std::getline(input_file, line);
    if (line.find("$PhysicalNames") != std::string::npos) {
      while (input_file.good() && !stop) {
        std::getline(input_file, line);
        if (line.find("$EndPhysicalNames") != std::string::npos) {
          stop = true;
          break;
        }
        if (line.length() == 1) {
          // n_physical_names = int(line);
        } else {
          int ignore, tag;
          std::string name;
          std::istringstream iss(line);
          iss >> ignore >> tag >> name;
          if (name.find("interface") != std::string::npos) {
            interface_ids.push_back(tag);
          } else if (name.find("bulk") != std::string::npos) {
            bulk_ids.push_back(tag);
          }
        }
      }
    }
  }
  input_file.close();

  // initialise the internal variables
  setup_quadrature_point_history();

  // This section names the different sides of the mesh, so that we
  // can assign different boundary conditions
  for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
           triangulation.begin_active();
       cell != triangulation.end(); ++cell) {
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
      if (cell->face(f)->at_boundary()) {
        const Point<dim> face_center = cell->face(f)->center();

        if (dim == 2) {
          if (std::abs(face_center[1] - (-1)) < 1e-8)
            cell->face(f)->set_boundary_id(1);
          else if (std::abs(face_center[1] - 1) < 1e-8)
            cell->face(f)->set_boundary_id(2);
          else if (std::abs(face_center[0] - (-1)) < 1e-8)
            cell->face(f)->set_boundary_id(3);
          else if (std::abs(face_center[0] - 1) < 1e-8)
            cell->face(f)->set_boundary_id(5);
          else
            cell->face(f)->set_boundary_id(0);
        } else if (dim == 3) {
          if (std::abs(face_center[1] - (-1)) < 1e-4)
            cell->face(f)->set_boundary_id(1);
          else if (std::abs(face_center[1] - 1) < 1e-4)
            cell->face(f)->set_boundary_id(2);
          else if (std::abs(face_center[0] - (-1)) < 1e-4)
            cell->face(f)->set_boundary_id(3);
          else if (std::abs(face_center[2] - (-1)) < 1e-4)
            cell->face(f)->set_boundary_id(4);
          else if (std::abs(face_center[0] - 1) < 1e-4)
            cell->face(f)->set_boundary_id(5);
          else
            cell->face(f)->set_boundary_id(0);
        }
      }
    }
  }

  std::ofstream out("../output/grid.msh");
  GridOut grid_out;
  grid_out.write_msh(triangulation, out);
}

template <int dim, int spacedim, bool num>
void TopLevel<dim, spacedim, num>::setup_system() {
  // distribute degrees of freedom over the mesh
  // dof_handler assumes that all elements are of type fe_bulk, with higher
  // polynomial degrees this leads to redundant dofs in the interface elements
  // Telling the dof_handler that there are different kinds of elements is
  // possible but useless, since this does not work for FEFaceSystem and even if
  // it did FEFaceSystem does not behave as it should anyways.
  dof_handler.distribute_dofs(fe_bulk);
  // use knowledge of the distribution of the dofs to create a
  // sparsitiy pattern fo the system of equations
  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
  // Allocate memory for the matrix and vectors of the global system of
  // equations
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
  residual.reinit(dof_handler.n_dofs());
}

template <int dim, int spacedim, bool num>
void TopLevel<dim, spacedim, num>::create_constraints() {

  // create the component_masks, which describe in which spatial direction a
  // BC is applied
  std::vector<bool> xmask, ymask, zmask, xymask;
  if (dim == 2) {
    ymask = std::vector<bool>{false, true},
    xmask = std::vector<bool>{true, false};
    xymask = std::vector<bool>{true, true};
  } else if (dim == 3) {
    xmask = std::vector<bool>{true, false, false},
    ymask = std::vector<bool>{false, true, false},
    zmask = std::vector<bool>{false, false, true};
  }

  // clear constraints of last timestep before new ones are applied
  constraints.clear();
  // Apply BC to face 1 (clamp the lower end in vertical direction)
  VectorTools::interpolate_boundary_values(dof_handler, 1,
                                           Functions::ZeroFunction<dim>(dim),
                                           constraints, ComponentMask(ymask));
  VectorTools::interpolate_boundary_values(dof_handler, 3,
                                           Functions::ZeroFunction<dim>(dim),
                                           constraints, ComponentMask(xmask));
  // Apply BC's to face 2 (displacement on the upper end)
  VectorTools::interpolate_boundary_values(
      dof_handler, 2,
      BoundaryValues<dim>(time.get_timestep(), time.get_no_timesteps(),
                          total_displacement),
      constraints, ComponentMask(ymask));

  // Activate for Biaxial stress
  // VectorTools::interpolate_boundary_values(dof_handler, 5,
  //     BoundaryValues<dim>(time.get_timestep(), time.get_no_timesteps(),
  //                         total_displacement),
  //                                          constraints,
  //                                          ComponentMask(xmask));

  // Apply BC's to face 3 (clamp the left hand side in horizontal direction)
  // VectorTools::interpolate_boundary_values(dof_handler, 3,
  //                                          Functions::ZeroFunction<dim>(dim),
  //                                          constraints,
  //                                          ComponentMask(ymask));

  if (dim == 3) {
    VectorTools::interpolate_boundary_values(dof_handler, 4,
                                             Functions::ZeroFunction<dim>(dim),
                                             constraints, ComponentMask(zmask));
  }

  // close the constraint object after all BC's have been generated
  constraints.close();
}

template <int dim, int spacedim, bool num>
void TopLevel<dim, spacedim, num>::assemble_system() {

  // clear the system of equations (set matrix and vectors to zero)
  system_matrix.reinit(sparsity_pattern);
  system_rhs.reinit(dof_handler.n_dofs());
  residual.reinit(dof_handler.n_dofs());

  // get the number of dofs_per_cell, this is equal in the case of
  // Lagrange polynomials of degree one but not for higher degrees
  const unsigned int dofs_per_cell = fe_bulk.dofs_per_cell;
  // Allocate memory for the cell contributions
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);
  std::tuple<FullMatrix<double>, Vector<double>> cell_contrib;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // loop over cells
  for (typename DoFHandler<dim, spacedim>::active_cell_iterator cell :
       dof_handler.active_cell_iterators()) {
    // get the local dof indices from the connectivity list
    cell->get_dof_indices(local_dof_indices);
    // get the nodal displacements per element
    Vector<double> Ue(dofs_per_cell);
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      Ue[i] = solution(local_dof_indices[i]);
    }
    // call the material routine for the bulk or interface element
    if (check_material_id(cell->material_id()) == 1) {
      cell_contrib =
          bulk.calc_cell_contrib(fe_bulk, cell, quadrature_formula_bulk, Ue);
    } else if (check_material_id(cell->material_id()) == 2) {
      auto quad_backup_1 = quadrature_point_history;
      cell_contrib =
          inter.calc_cell_contrib(fe_inter, cell, quadrature_formula_inter, Ue,
                                  quadrature_point_history);
      auto quad_backup_2 = quadrature_point_history;
      if (num) {
        auto Ubackup = Ue;
        Vector<double> cell_rhs_pert(dofs_per_cell);
        FullMatrix<double> cell_matrix_ana(dofs_per_cell, dofs_per_cell);
        std::tie(cell_matrix_ana, cell_rhs) = cell_contrib;
        for (unsigned int _i = 0; _i != dofs_per_cell; ++_i) {
          quadrature_point_history = quad_backup_1;
          Ue[_i] += std::max(pertubation_tangent, pertubation_tangent * Ue[_i]);
          cell_contrib =
              inter.calc_cell_contrib(fe_inter, cell, quadrature_formula_inter,
                                      Ue, quadrature_point_history);
          std::tie(std::ignore, cell_rhs_pert) = cell_contrib;
          for (unsigned int _j = 0; _j != dofs_per_cell; ++_j) {
            cell_matrix(_j, _i) =
                (cell_rhs_pert[_j] - cell_rhs[_j]) /
                std::max(pertubation_tangent, pertubation_tangent * Ue[_i]);
          }
          Ue = Ubackup;
          quadrature_point_history = quad_backup_2;
        }
        cell_contrib = std::make_tuple(cell_matrix, cell_rhs);
      }
    } else {
      cexc::not_mat_error exc;
      BOOST_THROW_EXCEPTION(exc);
    }
    // unpack the tuple returned by the material routine and write to
    // cell_matrix or cell_rhs
    std::tie(cell_matrix, cell_rhs) = cell_contrib;

    // this section arranges the cell matrix/rhs into the global
    // system of equations
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      for (unsigned int j = 0; j < dofs_per_cell; ++j) {
        system_matrix.add(local_dof_indices[i], local_dof_indices[j],
                          cell_matrix(i, j));
      }
      system_rhs(local_dof_indices[i]) -= cell_rhs(i);
    }
  }
}

template <int dim, int spacedim, bool num>
void TopLevel<dim, spacedim, num>::solve(unsigned int iter) {
  // damping_parameter = 1;
  // Create the Linear Solver
  // Use a direct sparse matrix solver
  SparseDirectUMFPACK solver;
  // reinitialise the update of the displacement vector to zero
  solution_update.reinit(dof_handler.n_dofs());
  // associate the solver and the system matrix
  solver.initialize(system_matrix);
  // Solve the linear system of equations
  solver.vmult(solution_update, residual);
  if (iter == 0)
    constraints.distribute(solution_update);
  else
    constraints.set_zero(solution_update);
  // perform newton update
  solution.add(damping_parameter, solution_update);
}

template <int dim, int spacedim, bool num>
void TopLevel<dim, spacedim, num>::output_results() {
  // This section creates an output
  DataOut<dim> data_out;
  // data_out.set_flags(write_higher_order_elements=false);
  // link the output to the dof_handler
  data_out.attach_dof_handler(dof_handler);
  // component_type defines the entries of the displacement field as vectors of
  // dim 2 or dim 3 Vectors
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
      component_type;
  if (dim == 2) {
    component_type = {DataComponentInterpretation::component_is_part_of_vector,
                      DataComponentInterpretation::component_is_part_of_vector};
  } else {
    component_type = {DataComponentInterpretation::component_is_part_of_vector,
                      DataComponentInterpretation::component_is_part_of_vector,
                      DataComponentInterpretation::component_is_part_of_vector};
  }
  // Create the output for the displacement field and RHS vector
  // Output is either 2D or 3D vectors
  data_out.add_data_vector(dof_handler, solution, "DSPL", component_type);
  data_out.add_data_vector(dof_handler, system_rhs, "RHS", component_type);
  Vector<double> mat_id(triangulation.n_active_cells());
  // create vector of material ids per cell
  int iter = 0;
  for (auto cell : dof_handler.active_cell_iterators()) {
    mat_id[iter] = cell->material_id();
    ++iter;
  }
  data_out.add_data_vector(mat_id, "Material_ID", DataOut<dim>::type_cell_data,
                           {DataComponentInterpretation::component_is_scalar});
  // build the output
  data_out.build_patches();

  // Generate filename based on timestep
  std::string filename = "../output/solution_";
  std::string timestep = std::to_string(time.get_timestep());
  int prec_zeros = 4 - timestep.length();
  for (auto i = 0; i < prec_zeros; ++i) {
    filename += "0";
  }
  filename += timestep + ".vtk";
  // write to vtk file
  std::ofstream output_vtk(filename);
  if (!output_vtk.is_open()) {
    cexc::file_error exc;
    BOOST_THROW_EXCEPTION(exc);
  }
  data_out.write_vtk(output_vtk);
}

template <int dim, int spacedim, bool num>
void TopLevel<dim, spacedim, num>::setup_quadrature_point_history() {
  // This section makes sure that quadrature_point_history is
  // empty and of the correct size
  unsigned int size = 0;
  for (auto cell : dof_handler.active_cell_iterators()) {
    triangulation.clear_user_data();
    if (check_material_id(cell->material_id()) == 1) {
      size += quadrature_formula_bulk.size();
    } else if (check_material_id(cell->material_id()) == 2) {
      size += quadrature_formula_inter.size();
    } else {
      cexc::not_imp_error exc;
      BOOST_THROW_EXCEPTION(exc);
    }
  }
  std::vector<PointHistory<dim>> tmp;
  tmp.swap(quadrature_point_history);
  quadrature_point_history.resize(size);

  // Link quadrature_point_history vector to quadrature points
  unsigned int history_index = 0;
  for (auto cell : dof_handler.active_cell_iterators()) {
    cell->set_user_pointer(&quadrature_point_history[history_index]);
    if (check_material_id(cell->material_id()) == 1) {
      history_index += quadrature_formula_bulk.size();
    } else if (check_material_id(cell->material_id()) == 2) {
      history_index += quadrature_formula_inter.size();
    } else {
      cexc::not_imp_error exc;
      BOOST_THROW_EXCEPTION(exc);
    }
  }
  Assert(history_index == quadrature_point_history.size(), ExcInternalError());
}

// This section demonstrates how to update the
// quadrature_point_history after the solution of the global system
// has been obtained  it has no real use in the current form of this
// code
template <int dim, int spacedim, bool num>
void TopLevel<dim, spacedim, num>::update_quadrature_point_history() {

  FEValues<dim, spacedim> fe_values(fe_bulk, quadrature_formula_bulk,
                                    update_values | update_gradients);
  std::vector<std::vector<Tensor<1, dim>>> displacement_increment_grads(
      quadrature_formula_bulk.size(), std::vector<Tensor<1, dim>>(dim));
  // loop over cells
  for (auto cell : dof_handler.active_cell_iterators()) {
    // get the local_quadrature_points_history of each cell
    PointHistory<dim> *local_quadrature_points_history =
        reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
    Assert(local_quadrature_points_history >= &quadrature_point_history.front(),
           ExcInternalError());
    Assert(local_quadrature_points_history < &quadrature_point_history.back(),
           ExcInternalError());

    // this section calculates the updated internal variables
    if (check_material_id(cell->material_id()) == 1) {
      for (unsigned int q = 0; q < quadrature_formula_bulk.size(); ++q) {
        unsigned int dofs_per_cell = fe_bulk.dofs_per_cell;
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
        Vector<double> Ue(dofs_per_cell);
        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
          Ue[i] = solution(local_dof_indices[i]);
        }
        SymmetricTensor<2, dim> stress =
            bulk.calc_stress(fe_bulk, cell, quadrature_formula_bulk, Ue, q);
        local_quadrature_points_history[q].bulk.old_stress = stress;
      }
    } else if (check_material_id(cell->material_id()) == 2) {
      for (unsigned int q = 0; q < quadrature_formula_inter.size(); ++q) {
        if (local_quadrature_points_history[q].inter.pen) {
          // local_quadrature_points_history[q].inter.alpha *= 2;
          local_quadrature_points_history[q].inter.pen = false;
        } else {
          local_quadrature_points_history[q].inter.alpha =
              local_quadrature_points_history[q].inter.alpha0;
        }
      }
    } else {
      cexc::not_imp_error exc;
      BOOST_THROW_EXCEPTION(exc);
    } // end cell case
  }   // end cell loop
}

template <int dim, int spacedim, bool num>
void TopLevel<dim, spacedim, num>::do_initial_timestep() {
  Timer t;
  std::cout << "doing initial timestep" << std::endl;
  make_grid();
  std::cout << "grid made" << std::endl;
  setup_system();
  std::cout << "initial step completed in " << t.elapsed() << " seconds "
            << std::endl;
  output_results();
}

template <int dim, int spacedim, bool num>
void TopLevel<dim, spacedim, num>::do_timestep() {
  Timer t_outer, t_inner;

  // increase the time
  time.increment();
  std::cout << "Timestep No. " << time.get_timestep() << " time "
            << time.get_current() << std::endl;
  // initialise norm of residual vector to one, so the loop is entered
  double rsn = 1.;
  double rsn_old = rsn;
  double step_length = 1e8;
  double old_step_length;
  unsigned int iter = 0;
  damping_parameter = 1;
  // create constraints
  create_constraints();
  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  constraints.condense(dsp);
  sparsity_pattern.copy_from(dsp);
  // Newton-Raphson loop
  auto res = residual;
  std::cout << "  Entering Newton Raphson Solver" << std::endl;
  while (rsn > tol_newton) {

    // if (iter > 1) {
    //   old_step_length = step_length;
    //   step_length = solution_update.l2_norm();
    //   if (step_length > old_step_length) {
    //     damping_parameter /= 2;
    //   } else {
    //     damping_parameter = 1;
    //   }
    // }
    t_inner.reset();

    // Assemble the system of equations
    assemble_system();
    residual = system_rhs;
    // std::cout << "    System assembled in " << t_inner.elapsed() << "
    // seconds" << std::endl;
    t_inner.reset();
    if (iter == 0) {
      constraints.condense(system_matrix, residual);
    } else {
      constraints.condense(system_matrix);
      constraints.set_zero(residual);
    }

    rsn = residual.l2_norm();
    std::cout << "    Iteration = " << iter << ",  current residual = " << rsn
              << std::endl;
    if (rsn > tol_newton) {
      t_inner.reset();
      solve(iter);
      // std::cout << "    Linear System solved in " << t_inner.elapsed() << "
      // seconds" << std::endl;
      iter++;
    }
    if (iter > max_iter_newton) {
      cexc::convergence_error exc;
      BOOST_THROW_EXCEPTION(exc);
    }
    std::cout << "    Newton Raphson Iteration Completed, Elapsed Time = "
              << t_outer.elapsed() << std::endl;
  }
  std::cout << "Step completed in " << t_outer.elapsed() << " seconds"
            << std::endl;
  update_quadrature_point_history();
  output_results();
}

template <int dim, int spacedim, bool num>
void TopLevel<dim, spacedim, num>::run() {
  do_initial_timestep();

  while (time.get_current() + 1e-8 < time.get_end()) {
    do_timestep();
  }
  //  output_results();
  // std::cout << "output results" << std::endl;
}
// FE Functions End ----------------------------------------------------

} // namespace fem

#endif