/**
 * @file DamInter.h
 * @author Till Budde (tilljanis.budde@tu-dortmund.de)
 * @brief
 * @version 0.1
 * @date 2021-09-15
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "bodyforce.h"
#include "pointhistory.h"
#include <algorithm>
#include <boost/exception/diagnostic_information.hpp>
#include <deal.II/base/point.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <tuple>
#include <vector>

/**
 * @brief Namespace of damage interface
 *
 */
namespace DamInter {

using namespace dealii;
using namespace fem;

/**
 * @brief Fiber damage model
 *
 * @tparam dim number of dimensions
 * @tparam spacedim  number of spatial dimensions
 */
template <int dim, int spacedim> class Material {

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
   * @brief Calculate Element contributions (Ke and RHS)
   *
   * @param fe  The finite Element (System)
   * @param cell  Pointer to the cell / DofCellAccessor
   * @param quadrature_formula  The quadrature formula object
   * @param Ue  The nodal displacements
   * @return std::tuple<FullMatrix<double>,Vector<double>> Return the \ref
   * cell_matrix as first and the \ref cell_rhs as second argument
   *
   * This function loops over all quadrature points.
   * At each point the jump / traction is calculated before a loop over
   * the dofs of the cell assembles the elemental stiffness matrix and
   * right hand side vector.
   */
  std::tuple<FullMatrix<double>, Vector<double>> calc_cell_contrib(
      FESystem<dim, spacedim> &fe,
      typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
      QGauss<dim - 1> quadrature_formula, Vector<double> Ue,
      std::vector<PointHistory<dim>> &quadrature_point_history);

  /**
   * @brief Calculate the stress Tensor
   *
   * @param SymmetricTensor<2,dim> Strain Tensor
   * @return SymmetricTensor<2,dim> Stress Tensor
   */
  SymmetricTensor<2, dim> calc_stress(SymmetricTensor<2, dim> strain);

private:
  /**
   * amplitude of displacement jump at onset of damage
   */
  double kappa_0 = 0.005;
  double stiffness = 210; /*!< the stiffness of the virgin spring */
  double Q_0 = 100;       /*!< strength of the interface */
  double Q_f = 50;        /*!< fracture energy */
};

template <int dim, int spacedim>
std::tuple<FullMatrix<double>, Vector<double>>
Material<dim, spacedim>::calc_cell_contrib(
    FESystem<dim, spacedim> &fe,
    typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    QGauss<dim - 1> quadrature_formula, Vector<double> Ue,
    std::vector<PointHistory<dim>> &quadrature_point_history) {
  // the constructor of feface_values needs fe & quadrature_formula,
  // the other arguments determine which values are calculated by
  // fe_values(reinit)
  // FEFaceValues could be constructed outside the material box and
  // given as an argument for better computational efficency
  // [ ] TODO ^
  FEFaceValues<dim, spacedim> fe_values(fe, quadrature_formula,
                                        update_values | update_JxW_values |
                                            update_quadrature_points);
  // get the number of dofs and quadrature points
  const unsigned int n_q_points = quadrature_formula.size();
  const unsigned int n_faces = cell->n_faces();
  // The FeFaceSystem class is discontinues at the connections between
  // different faces. Meaning that the nodes at the edges have one set of
  // degrees of freedom for each face it is a part of.
  // Since the model expects only Degrees of Freedom at the upper (+) and
  // lower (-) side, the correct numer of dofs has to be calculated like this
  const unsigned int dofs_per_cell = fe.dofs_per_cell / n_faces * 2;
  // allocate memory for the cell_matrix and cell_rhs, as well as finte,fvole
  FullMatrix<double> cell_matrix(dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);
  // allocate memory for the damage and damage_history variables
  double damage, kappa;

  // The DofHandler::active_cell_iterator points to an object of the bulk
  // material type. Using this to call fe_values.reinit() leads to an
  // exception throw, since the interface material is of a different type
  // and has different number of dofs and different shape functions.
  // By casting the pointer into the more general
  // Triangulation::active_cell_iterator we can circumnavigate this problem,
  // since the Triangulation has no knowledge of dofs anyways.
  typename Triangulation<dim, spacedim>::active_cell_iterator t = cell;

  // get the local_quadrature_points_history of each cell
  PointHistory<dim> *local_quadrature_points_history =
      reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
  Assert(local_quadrature_points_history >= &quadrature_point_history.front(),
         ExcInternalError());
  Assert(local_quadrature_points_history < &quadrature_point_history.back(),
         ExcInternalError());

  // Calculate the stiffness matrix
  auto identity = Physics::Elasticity::StandardTensors<dim>::I;
  auto C = stiffness * identity;
  // allocate memory for the jump and traction
  Tensor<1, dim, double> jump;
  Tensor<1, dim, double> traction;

  // loop over quadrature points
  for (unsigned int q_point = 0; q_point != n_q_points; ++q_point) {

    // reset jump to zero
    jump = 0;

    /**
     * It should be noted, that in general the FeFaceValues class loops over
     * its dofs by going [x_1,x_2,y_1,_y2] for the first face and then goes on
     * for every other face.
     * This has two problems:
     * 	- 1. the number of dofs is higher than the number of dofs of the
     * actual interface element.
     * 	- 2. This is inconsistent with the order in which fevalues iterates over
     * its dofs, going [x_1,y_1,x_2,y_2 ...]. This is the order in which
     * the displacement vector \ref Ue as well as the \ref cell_rhs and the
     * \ref cell_matrix are ordered.
     *
     * In order to solve this problem the program loops over faces and then over
     * the dofs per face and calculates the coressponding vector_index and
     * shape function by use of some of dealii functions that allow the
     * calculation of nodes, vector components and cell_indices.
     *
     */

    // loop over the upper (n_faces -1) and lower face (n_faces-2)
    for (unsigned int face : {n_faces - 2, n_faces - 1}) {

      // calculate the values of the master element
      fe_values.reinit(t, face);
      unsigned int dofs_per_face = fe.n_dofs_per_face(face, 0);
      // loop over the dofs of the face
      for (unsigned int face_dof = 0; face_dof != dofs_per_face; ++face_dof) {
        // get the direction / vector component of the jump (x,y,z)
        auto component_i =
            fe.face_system_to_component_index(face_dof, face).first;
        // get the associated node
        auto node_i = fe.face_system_to_component_index(face_dof, face).second;
        auto vector_index =
            dofs_per_face * (face - (n_faces - 2)) + auto vector_index =
                dofs_per_face * (face - (n_faces - 2)) +
                node_i * fe.n_components() + component_i;
        // calculate jump
        jump[component_i] +=
            Ue[vector_index] *
            fe_values.shape_value(
                fe.face_to_cell_index(
                    face_dof, face, cell->face_orientation(face), false, false),
                q_point) *
            ((face == n_faces - 2) ? -1 : 1);
      }
    }

    // get damage_history variable of last timestep
    kappa = local_quadrature_points_history[q_point].inter.kappa;
    // check is current jump is greater than last timestep
    kappa = std::max(kappa, jump.norm());
    // calculate damage
    if (kappa < kappa_0) {
      damage = 0;
    } else {
      damage = 1 - kappa_0 / kappa * exp(-(kappa - kappa_0) * Q_0 / Q_f);
    }

    // calculate traction
    traction = (1 - damage) * C * jump;

    // first loop over faces
    for (unsigned int face_i : {n_faces - 2, n_faces - 1}) {

      // loop over dofs of first face
      unsigned int dofs_per_face = fe.n_dofs_per_face(face_i, 0);
      for (unsigned int face_dof_i = 0; face_dof_i != dofs_per_face;
           ++face_dof_i) {
        // reinit values of the master element
        fe_values.reinit(t, face_i);
        // get the direction / vector component (x,y,z)
        auto component_i =
            fe.face_system_to_component_index(face_dof_i, face_i).first;
        // get the associated node
        auto node_i =
            fe.face_system_to_component_index(face_dof_i, face_i).second;
        // get the index of the displacement vector
        auto index_i = dofs_per_face * (face_i - (n_faces - 2)) +
                       node_i * fe.n_components() + component_i;
        // get the value of the shape function
        double N_i = fe_values.shape_value(
            fe.face_to_cell_index(face_dof_i, face_i,
                                  cell->face_orientation(face_i), false, false),
            q_point);
        // get the index of the cell_rhs vector / first index of cell_matrix
        auto vector_index = dofs_per_face * (face_i - (n_faces - 2)) +
                            node_i * fe.n_components() + component_i;

        // calculate the internal force fector
        auto fint_i = fe_values.shape_value(
                          fe.face_to_cell_index(face_dof_i, face_i,
                                                cell->face_orientation(face_i),
                                                false, false),
                          q_point) *
                      traction[component_i] * fe_values.JxW(q_point) *
                      ((face_i == n_faces - 2) ? -1 : +1);

        // add contributions to cell_rhs
        cell_rhs(vector_index) += fint_i;
        AssertIsFinite(cell_rhs(vector_index));

        // second loop over faces
        for (unsigned int face_j : {n_faces - 2, n_faces - 1}) {
          // reinit values of master element to second face
          fe_values.reinit(t, face_j);
          // loop over dofs of face
          for (unsigned int face_dof_j = 0; face_dof_j != dofs_per_face;
               ++face_dof_j) {
            // get the direction / vector component (x,y,z)
            auto component_j =
                fe.face_system_to_component_index(face_dof_j, face_j).first;
            // get the associated node
            auto node_j =
                fe.face_system_to_component_index(face_dof_j, face_j).second;
            // get the second index of the cell_matrix
            auto index_j = dofs_per_face * (face_j - (n_faces - 2)) +
                           node_j * fe.n_components() + component_j;
            // get the value of the shape function
            double N_j = fe_values.shape_value(
                fe.face_to_cell_index(face_dof_j, face_j,
                                      cell->face_orientation(face_j), false,
                                      false),
                q_point);

            // add contribution to stiffness matrix
            cell_matrix(index_i, index_j) +=
                (1 - damage) * N_i * C[component_i][component_j] * N_j *
                fe_values.JxW(q_point) * ((face_i == n_faces - 2) ? 1 : -1) *
                ((face_j == n_faces - 2) ? 1 : -1);
            AssertIsFinite(cell_matrix(index_i, index_j));
          } // close face_dof_j
        }   // close face_j
      }     // close_face_dof_i
    }       // close_face_i
  }         // quadrature point end

  return std::make_tuple(cell_matrix, cell_rhs);
}

} // namespace DamInter