/**
 * @file LinElaInter.h
 * @author Till Janis Budde
 * @brief  Generalized Hyperelastic interface model
 * @version 0.1
 * @date 2021-05-19
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "LagrangePolynomial.h"
#include "bodyforce.h"
#include "pointhistory.h"
#include <boost/exception/diagnostic_information.hpp>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <tuple>
#include <vector>

/**
 * @brief Namespace of Generalized Hyperelastic interface model
 *
 */
namespace GenElaInter {

using namespace dealii;
using namespace fem;

/**
 * @brief Generalized Hyperelastic interface model
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
      QGauss<dim - 1> quadrature_formula, Vector<double> Ue);

  /**
   * @brief Calculate the stress Tensor
   *
   * @param SymmetricTensor<2,dim> Strain Tensor
   * @return SymmetricTensor<2,dim> Stress Tensor
   */
  SymmetricTensor<2, dim> calc_stress(SymmetricTensor<2, dim> strain);

  void calc_derivatives(Tensor<1, dim> u, Tensor<1, dim> g1, Tensor<1, dim> g2,
                        Tensor<1, dim> g3);

  std::tuple<Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>,
             Tensor<2, dim>, Tensor<2, dim>, Tensor<2, dim>, Tensor<2, dim>,
             Tensor<2, dim>, Tensor<2, dim>, Tensor<2, dim>, Tensor<2, dim>,
             Tensor<2, dim>, Tensor<2, dim>>
  calc_derivatives(Tensor<1, dim> u, Tensor<1, dim> g1, Tensor<1, dim> g2);

private:
  double c_n = 2e2;
  double c_t = 10e2;
  double fa_I1 = 1, fa_I2 = 1, fa_I3 = 1;
  double stiffness = 210; /*!< The stiffness of the spring */
};

template <int dim, int spacedim>
std::tuple<FullMatrix<double>, Vector<double>>
Material<dim, spacedim>::calc_cell_contrib(
    FESystem<dim, spacedim> &fe,
    typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    QGauss<dim - 1> quadrature_formula, Vector<double> Ue) {

  // the constructor of feface_values needs fe & quadrature_formula,
  // the other arguments determine which values are calculated by
  // fe_values(reinit)
  // FEFaceValues could be constructed outside the material box and
  // given as an argument for better computational efficency
  // [ ] TODO ^
  FEFaceValues<dim, spacedim> fe_values(fe, quadrature_formula,
                                        update_values | update_JxW_values |
                                            update_normal_vectors |
                                            update_quadrature_points);
  LagrangePolynomial<dim - 1> gamma(quadrature_formula.size() / (dim - 1));
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

  // The DofHandler::active_cell_iterator points to an object of the bulk
  // material type. Using this to call fe_values.reinit() leads to an
  // exception throw, since the interface material is of a different type
  // and has different number of dofs and different shape functions.
  // By casting the pointer into the more general
  // Triangulation::active_cell_iterator we can circumnavigate this problem,
  // since the Triangulation has no knowledge of dofs anyways.
  typename Triangulation<dim, spacedim>::active_cell_iterator t = cell;

  // Calculate the stiffness matrix
  auto identity = Physics::Elasticity::StandardTensors<dim>::I;
  auto C = stiffness * identity;
  // allocate memory for the jump and traction
  Tensor<1, dim, double> jump;
  Tensor<1, dim, double> traction;
  Tensor<1, dim, double> normal_vector;
  Tensor<1, dim, double> grad;
  Tensor<1, dim, double> g1;
  Tensor<1, dim, double> g2;

  // loop over quadrature points
  for (unsigned int q_point = 0; q_point != n_q_points; ++q_point) {

    // reset jump to zero
    jump = 0;
    normal_vector = 0;
    g1 = 0;
    g2 = 0;

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
        // get the index of the displacement vector
        auto vector_index = dofs_per_face * (face - (n_faces - 2)) +
                            node_i * fe.n_components() + component_i;
        // calculate jump
        jump[component_i] +=
            Ue[vector_index] *
            fe_values.shape_value(
                fe.face_to_cell_index(
                    face_dof, face, cell->face_orientation(face), false, false),
                q_point) *
            ((face == n_faces - 2) ? -1 : 1);
        // g1 = gamma.get_first_derivative(node_i, q_point);
        g1[component_i] +=
            (0.5 * gamma.get_first_derivative(node_i, q_point)[0] *
             (cell->vertex(node_i)[component_i] + Ue[vector_index]));
        if (dim == 3) {
          g2[component_i] +=
              (0.5 * gamma.get_first_derivative(node_i, q_point)[1] *
               (cell->vertex(node_i)[component_i] + Ue[vector_index]));
        }
      } // end loop face_dofs

      normal_vector +=
          fe_values.normal_vector(q_point) * ((face == n_faces - 2) ? -1 : 1);
    } // end loop faces

    /*--------------------------------------------------------------------------
    Calculate tangential and normal base vectors (g) of the interface
    --------------------------------------------------------------------------*/
    normal_vector /= normal_vector.norm();
    // calculate traction
    traction = C * jump;

    Tensor<1, dim> Psi_u, Psi_g1, Psi_g2, Psi_g3;
    Tensor<2, dim> Psi_u_u, Psi_u_g1, Psi_u_g2, Psi_u_g3, Psi_g1_g1, Psi_g1_g2,
        Psi_g1_g3, Psi_g2_g2, Psi_g2_g3, Psi_g3_g3;
    { // this scope is only so my editor can fold this region
      Psi_u = 0;
      Psi_g1 = 0;
      Psi_g2 = 0;
      Psi_g3 = 0;
      Psi_u_u = 0;
      Psi_u_g1 = 0;
      Psi_u_g2 = 0;
      Psi_u_g3 = 0;
      Psi_g1_g1 = 0;
      Psi_g1_g2 = 0;
      Psi_g1_g3 = 0;
      Psi_g2_g2 = 0;
      Psi_g2_g3 = 0;
      Psi_g3_g3 = 0;
    }
    if (dim == 2) {
      std::tie(Psi_u, Psi_g1, Psi_g2, Psi_g3, Psi_u_u, Psi_u_g1, Psi_u_g2,
               Psi_u_g3, Psi_g1_g1, Psi_g1_g2, Psi_g1_g3, Psi_g2_g2, Psi_g2_g3,
               Psi_g3_g3) = calc_derivatives(jump, g1, normal_vector);
    } else if (dim == 3) {
      calc_derivatives(jump, g1, g2, normal_vector);
    } else {
    }

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

        Tensor<1, dim - 1> gamma_i =
            gamma.get_first_derivative(node_i, q_point);
        /*----------------------------------------------------------------------
        Calculate derivative of Psi w.r.t [[u]] and g
        ----------------------------------------------------------------------*/

        // calculate the internal force fector
        double fint_i =
            (N_i * ((face_i == n_faces - 2) ? -1 : 1) * Psi_u[component_i] +
             0.5 * (gamma_i[0] * Psi_g1[component_i] +
                    ((dim == 3) ? gamma_i[1] : 0) * Psi_g2[component_i])) *
            fe_values.JxW(q_point);

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

            Tensor<1, dim - 1> gamma_j =
                gamma.get_first_derivative(node_j, q_point);

            // add contribution to stiffness matrix
            cell_matrix(index_i, index_j) +=
                (N_i * ((face_i == n_faces - 2) ? -1 : 1) *
                     (Psi_u_u[component_i][component_j] * N_j *
                          ((face_j == n_faces - 2) ? -1 : 1) +
                      // Psi_u_g1 transponiert - deckt sich nicht
                      // mit herleitung - passt zu numerischer
                      // tangente
                      Psi_u_g1[component_j][component_i] * 0.5 * gamma_j[0] +
                      Psi_u_g2[component_i][component_j] * 0.5 *
                          ((dim == 3) ? gamma_j[1] : 0)) +
                 0.5 * gamma_i[0] *
                     (Psi_g1_g1[component_i][component_j] * gamma_j[0] +
                      Psi_g1_g2[component_i][component_j] *
                          ((dim == 3) ? gamma_j[1] : 0) +
                      // Psi_u_g1 nicht transponiert - deckt sich
                      // nicht mit herleitung - passt zu numerischer
                      // tangente
                      Psi_u_g1[component_i][component_j] * N_j *
                          ((face_j == n_faces - 2) ? -1 : 1)) +
                 0.5 * ((dim == 3) ? gamma_i[1] : 0) *
                     (transpose(Psi_g1_g2)[component_i][component_j] *
                          gamma_j[0] +
                      Psi_g2_g2[component_i][component_j] *
                          ((dim == 3) ? gamma_j[1] : 0) +
                      transpose(Psi_u_g2)[component_i][component_j] * N_j *
                          ((face_j == n_faces - 2) ? -1 : 1))) *
                fe_values.JxW(q_point);
            AssertIsFinite(cell_matrix(index_i, index_j));
          } // close face_dof_j
        }   // close face_j
      }     // close_face_dof_i
    }       // close_face_i
  }         // quadrature point end

  return std::make_tuple(cell_matrix, cell_rhs);
}

template <int rank_1, int rank_2, int dim, typename Number,
          typename OtherNumber>
constexpr inline DEAL_II_ALWAYS_INLINE typename Tensor<
    rank_1 + rank_2 - 2, dim,
    typename ProductType<Number, OtherNumber>::type>::tensor_type
single_contract(const Tensor<rank_1, dim, Number> &src1,
                const Tensor<rank_2, dim, OtherNumber> &src2) {

  return contract<src1.rank - 1, 0>(src1, src2);
}

template <int dim, int spacedim>
void Material<dim, spacedim>::calc_derivatives(Tensor<1, dim> u,
                                               Tensor<1, dim> g1_cov,
                                               Tensor<1, dim> g2_cov,
                                               Tensor<1, dim> g3_cov) {

  Tensor<1, dim, double> g1_con, g2_con, g3_con, g1xg3, g2xg3;

  //   g1xg3 = cross_product_3d(g1_cov, g3_cov);
  //   g2xg3 = cross_product_3d(g2_cov, g3_cov);

  //   g1_con = (g2xg3) / (scalar_product(g1_cov, g2xg3));
  //   g2_con = (g1xg3) / (scalar_product(g2_cov, g1xg3));
  //   g3_con = g3_cov;

  //   double I1 = scalar_product(u, u);
}

template <int dim, int spacedim>
std::tuple<Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>,
           Tensor<2, dim>, Tensor<2, dim>, Tensor<2, dim>, Tensor<2, dim>,
           Tensor<2, dim>, Tensor<2, dim>, Tensor<2, dim>, Tensor<2, dim>,
           Tensor<2, dim>, Tensor<2, dim>>
Material<dim, spacedim>::calc_derivatives(Tensor<1, dim> u,
                                          Tensor<1, dim> g1_cov,
                                          Tensor<1, dim> g3_cov) {

  // This function assumes a two dimensional space
  assert(dim == 2);

  Tensor<1, dim, double> g1_con, g3_con;

  g1_con = g1_cov / (g1_cov.norm() * g1_cov.norm());
  g3_con = g3_cov;

  double g11 = single_contract(g1_con, g1_con);
  //   double g33 = single_contract(g3_con, g3_con);

  //   double I1 = single_contract(u, u);
  //   double I2 = single_contract(u, g1_cov);

  Tensor<2, dim> u_dy_g1_con = outer_product(u, g1_con);
  Tensor<2, dim> u_dy_g1_cov = outer_product(u, g1_cov);
  Tensor<2, dim> u_dy_g3_con = outer_product(u, g3_con);
  Tensor<2, dim> u_dy_g3_cov = outer_product(u, g3_cov);

  Tensor<2, dim> g1_con_dy_g1_con = outer_product(g1_con, g1_con);
  Tensor<2, dim> g3_con_dy_g3_con = outer_product(g3_con, g3_con);

  double u_sc_g1_con = single_contract(u, g1_con);
  double u_sc_g3_con = single_contract(u, g3_con);
  double u_sc_g1_cov = single_contract(u, g1_cov);
  //   double u_sc_g3_cov = single_contract(u, g3_cov);

  /*----------------------------------------------------------------------------
  Psi = 1/2 ( c_n *(I1-I2) + c_t(I2))
  ----------------------------------------------------------------------------*/
  Tensor<1, dim> Psi_u, Psi_g1, Psi_g2, Psi_g3;
  double Psi_I1, Psi_I2;

  Psi_I1 = 0.5 * c_n;
  Psi_I2 = 0.5 * (-c_n + c_t);

  Tensor<1, dim> I1_u = 2 * u;
  Tensor<1, dim> I2_u = u_sc_g1_con * g1_cov + u_sc_g1_cov * g1_con;

  Tensor<1, dim> I1_g1;
  I1_g1 = 0;
  Tensor<1, dim> I2_g1 =
      u_sc_g1_con * u +
      u_sc_g1_cov * (-u_sc_g1_con * g1_con + g11 * u_sc_g3_con * g3_con);

  Psi_u = fa_I1 * Psi_I1 * I1_u + fa_I2 * Psi_I2 * I2_u;
  Psi_g1 = fa_I1 * Psi_I1 * I1_g1 + fa_I2 * Psi_I2 * I2_g1;
  Psi_g2 = 0;
  Psi_g3 = 0;

  auto identity = Physics::Elasticity::StandardTensors<dim>::I;

  Tensor<2, dim> Psi_u_u, Psi_u_g1, Psi_u_g2, Psi_u_g3, I1_u_u, I2_u_u, I2_u_g1,
      I2_g1_u, I2_g1_g1, Psi_g1_g1, Psi_g1_g2, Psi_g1_g3, Psi_g2_g2, Psi_g2_g3,
      Psi_g3_g3;

  I1_u_u = 2 * identity;
  I2_u_u = outer_product(g1_con, g1_cov) + outer_product(g1_cov, g1_con);

  I2_u_g1 = u_dy_g1_con + u_sc_g1_con * identity +
            outer_product(-u_sc_g1_con * g1_con + g11 * u_sc_g3_con * g3_con,
                          g1_cov) +
            u_sc_g1_cov * (-g1_con_dy_g1_con + g11 * g3_con_dy_g3_con);
  I2_g1_u = transpose(I2_u_g1);

  I2_g1_g1 =
      outer_product(u, (-u_sc_g1_con * g1_con + g11 * u_sc_g3_con * g3_con)) +
      outer_product((-u_sc_g1_con * g1_con + g11 * u_sc_g3_con * g3_con), u) +
      u_sc_g1_con * (outer_product(u_sc_g1_con * g1_con, g1_con) -
                     g11 * outer_product(u_sc_g3_con * g1_con, g3_con) +
                     u_sc_g1_con * (g1_con_dy_g1_con - g11 * g3_con_dy_g3_con) -
                     2 * (g11 * outer_product(u_sc_g3_con * g3_con, g1_con)) -
                     g11 * outer_product(u_sc_g1_con * g3_con, g3_con) -
                     g11 * outer_product(u_sc_g3_con * g1_con, g3_con));

  Psi_u_u = Psi_I1 * I1_u_u + Psi_I2 * I2_u_u;
  Psi_u_g1 = Psi_I2 * I2_u_g1;
  Psi_u_g2 = 0;
  Psi_u_g3 = 0;

  Psi_g1_g1 = Psi_I2 * I2_g1_g1;
  Psi_g1_g2 = 0;
  Psi_g1_g3 = 0;
  Psi_g2_g2 = 0;
  Psi_g2_g3 = 0;
  Psi_g3_g3 = 0;

  return {Psi_u,     Psi_g1,    Psi_g2,    Psi_g3,    Psi_u_u,
          Psi_u_g1,  Psi_u_g2,  Psi_u_g3,  Psi_g1_g1, Psi_g1_g2,
          Psi_g1_g3, Psi_g2_g2, Psi_g2_g3, Psi_g3_g3};
}

} // namespace GenElaInter
