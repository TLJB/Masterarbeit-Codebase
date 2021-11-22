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
#ifndef NeoHooke
#define NeoHooke

#include "bodyforce.h"
#include "pointhistory.h"
#include <deal.II/base/tensor.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/numerics/vector_tools.h>
#include <tuple>

/**
 * @brief Namespace for small-strain elasticity
 *
 */
namespace NeoHooke {

using namespace dealii;
using namespace fem;

/**
 * @brief Small Strain Linear Elasticity
 *
 * @tparam dim number of dimensions
 * @tparam spacedim number of spatial dimensions
 */
template <int dim, int spacedim> class Material {
public:
  /**
   * @brief Construct a new Material object
   *
   * Initialise the 4th order stress-strain tensor
   *
   */
  Material();
  /**
   * @brief Destroy the Material object
   *
   */
  ~Material();

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
   * At each point the strain / stress is calculated before a loop over
   * the dofs of the cell assembles the elemental stiffness matrix and
   * right hand side vector.
   */
  std::tuple<FullMatrix<double>, Vector<double>> calc_cell_contrib(
      FESystem<dim, spacedim> &fe,
      typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
      QGauss<dim> quadrature_formula, Vector<double> Ue);

  /**
   * @brief
   *
   * @param fe  The finite Element (System)
   * @param cell  Pointer to the cell / DofCellAccessor
   * @param quadrature_formula  The quadrature formula object
   * @param Ue  The nodal displacements
   * @return SymmetricTensor<2,dim> The stress tensor
   */
  SymmetricTensor<2, dim>
  calc_stress(FESystem<dim, spacedim> &fe,
              typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
              QGauss<dim> quadrature_formula, Vector<double> Ue,
              unsigned int q_point);

private:
  SymmetricTensor<4, dim>
      E_tensor;  /*!< The isotropic elastic stress strain tensor */
  double lambda; /*!< The first Lame parameter */
  double mu;     /*!< The second Lame parameter */
  double youngs_modulus = 210; /*!< The youngs modulus */
  double poisson_ratio = 0.33; /*!< The poisson ratio */
};

// Material Functions --------------------------------------------------

template <int dim, int spacedim> Material<dim, spacedim>::Material() {
  // calculate the first Lame parameter
  lambda = (youngs_modulus * poisson_ratio) /
           ((1 + poisson_ratio) * (1 - 2 * poisson_ratio));
  // calculate the second Lame parameter
  mu = youngs_modulus / (2 * (1 + poisson_ratio));
}

template <int dim, int spacedim> Material<dim, spacedim>::~Material() {}

template <int dim, int spacedim>
std::tuple<FullMatrix<double>, Vector<double>>
Material<dim, spacedim>::calc_cell_contrib(
    FESystem<dim, spacedim> &fe,
    typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    QGauss<dim> quadrature_formula, Vector<double> Ue) {
  // the constructor of fe_values needs fe & quadrature_formula,
  // the other arguments determine which values are calculated by
  // fe_values(reinit)
  // FEValues could be constructed outside the material box and
  // given as an argument for better computational efficency
  // [ ] TODO ^
  FEValues<dim, spacedim> fe_values(fe, quadrature_formula,
                                    update_values | update_gradients |
                                        update_quadrature_points |
                                        update_JxW_values);
  // get the number of dofs and quadrature points
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();
  // allocate memory for the cell_matrix and cell_rhs, as well as finte,fvole
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell), finte(dofs_per_cell),
      fvole(dofs_per_cell);

  // calculate the values of the master element
  fe_values.reinit(cell);
  // Calulate the bodyforce values
  BodyForce<dim> body_force;
  std::vector<Vector<double>> body_force_values(n_q_points,
                                                Vector<double>(dim));
  body_force.vector_value_list(fe_values.get_quadrature_points(),
                               body_force_values);

  // loop over all quadrature_points
  Tensor<2, dim> F;
  for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {

    // calculate the strain tensor
    F = 0;
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int component_i = fe.system_to_component_index(i).first;
      for (unsigned int j = 0; j < dim; ++j) {

        F[component_i][j] += Ue[i] * (fe_values.shape_grad(i, q_point))[j];
      }
    }
    auto identity = Physics::Elasticity::StandardTensors<dim>::I;
    F = F + identity;
    Tensor<2, dim> Finv = invert(F);
    Tensor<2, dim> Finvt = transpose(Finv);
    double J = determinant(F);

    // calculate the stress
    Tensor<2, dim> stress = (lambda * log(J) - mu) * Finvt + mu * F;

    // calculate Elasticity Tensor
    for (unsigned int i = 0; i != dim; ++i) {
      for (unsigned int j = 0; j != dim; ++j) {
        for (unsigned int k = 0; k != dim; ++k) {
          for (unsigned int l = 0; l != dim; ++l) {
            E_tensor[i][j][k][l] =
                mu * identity[i][k] * identity[j][l] +
                lambda * Finvt[i][j] * Finvt[k][l] -
                (lambda * log(J) - mu) * Finvt[i][l] * Finvt[j][k];
          }
        }
      }
    }

    // reset the values for finte/ fvole for each quadrature point
    finte = 0;
    fvole = 0;
    // loop over dofs - note that i=0 is the first dof of the
    // first node and i=1 is the second dof of the first node
    // deal needs this format to rearrange the cell_matrix into
    // the global stiffness matrix
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int component_i = fe.system_to_component_index(i).first;
      for (unsigned int j = 0; j < dofs_per_cell; ++j) {
        const unsigned int component_j = fe.system_to_component_index(j).first;
        // assemble the cell matrix
        cell_matrix(i, j) += (fe_values.shape_grad(i, q_point) * E_tensor *
                              fe_values.shape_grad(j, q_point) *
                              fe_values.JxW(q_point))[component_i][component_j];
        AssertIsFinite(cell_matrix(i, j));
      } // close loop over second node

      // calculate the body force
      fvole(i) += -body_force_values[q_point](component_i) *
                  fe_values.shape_value(i, q_point) * fe_values.JxW(q_point);
      // calculate the internal force fector
      finte(i) += (fe_values.shape_grad(i, q_point) * stress *
                   fe_values.JxW(q_point))[component_i];
    } // close node over first node

    // add contributions to cell_rhs
    cell_rhs.add(1, finte, -1, fvole);
    for (auto rhs_val : cell_rhs)
      AssertIsFinite(rhs_val);
  }

  return std::make_tuple(cell_matrix, cell_rhs);
}

template <int dim, int spacedim>
SymmetricTensor<2, dim> Material<dim, spacedim>::calc_stress(
    FESystem<dim, spacedim> &fe,
    typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    QGauss<dim> quadrature_formula, Vector<double> Ue, unsigned int q_point) {

  // See calc_cell_matrix
  FEValues<dim, spacedim> fe_values(
      fe, quadrature_formula,
      update_values | update_gradients | update_quadrature_points |
          update_inverse_jacobians | update_JxW_values);
  fe_values.reinit(cell);
  const unsigned int n_q_points = quadrature_formula.size();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  Vector<double> cell_rhs(dofs_per_cell), finte(dofs_per_cell),
      fvole(dofs_per_cell);

  // this section gets the bodyforce values
  BodyForce<dim> body_force;
  std::vector<Vector<double>> body_force_values(n_q_points,
                                                Vector<double>(dim));
  body_force.vector_value_list(fe_values.get_quadrature_points(),
                               body_force_values);

  // // assemble rhs by looping over dofs and q_points
  SymmetricTensor<2, dim> strain;
  strain = 0;
  for (unsigned int i = 0; i < dofs_per_cell; ++i) {
    const unsigned int component_i = fe.system_to_component_index(i).first;
    for (unsigned int j = 0; j < dim; ++j) {

      strain[component_i][j] += Ue[i] * (fe_values.shape_grad(i, q_point))[j];
    }
  }
  strain = 0.5 * (strain + transpose(strain));

  SymmetricTensor<2, dim> stress = E_tensor * strain;
  return stress;
}

// Material Functions End ----------------------------------------------

} // namespace NeoHooke

#endif