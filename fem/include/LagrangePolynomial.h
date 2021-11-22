/**
 * @file LagrangePolynomial.h
 * @author Till Budde (tilljanis.budde@tu-dortmund.de)
 * @brief Implementation of LagrangePolynomials
 * @version 0.1
 * @date 2021-10-05
 *
 * @copyright Copyright (c) 2021
 *
 * This class is explicitly meant to substitute
 * FEFaceValues<dim,spacedim>::shape_grad(int,int) as this function is not
 * implemented. Consequently its only use is to return the derivative of the
 * shape function N with respect to the generalized coordinates \xi.
 *
 * It should be noted that the generalized coordinates of the master element
 * in dealii run from from 0 to 1 (\xi \in [0,1]) and not from the usual
 * -1 to 1.
 *
 */

#ifndef LAGRANGEPOLYNOMIAL_H
#define LAGRANGEPOLYNOMIAL_H

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include "CustomExceptions.h"
#include <boost/throw_exception.hpp>

#include <vector>

namespace fem {

using namespace dealii;

/**
 * @brief Implement the LagrangePolynomials and its derivative
 *
 * @tparam dim
 *
 * This class is explicitly meant to substitute
 * FEFaceValues<dim,spacedim>::shape_grad(int,int) as this function is not
 * implemented. Consequently its only use is to return the derivative of the
 * shape function N with respect to the generalized coordinates \xi.
 *
 * It should be noted that the generalized coordinates of the master element
 * in dealii run from from 0 to 1 (\xi \in [0,1]) and not from the usual
 * -1 to 1.
 */
template <int dim> class LagrangePolynomial {

public:
  /**
   * @brief Construct a new Lagrange Polynomial object
   *
   * The constructor defaults to 2 quadrature points
   *
   */
  LagrangePolynomial() : quadrature(2) { std::cout << init(1); };
  /**
   * @brief Construct a new Lagrange Polynomial object
   *
   * @param q The number of quadrature points
   */
  LagrangePolynomial(int q) : quadrature(q) { init(1); };
  /**
   * @brief Construct a new Lagrange Polynomial object
   * 
   * @param q the number of quadrature points
   * @param p the degree of the lagrange polynomial
   */
  LagrangePolynomial(int q, int p) : quadrature(q) { init(p); };
  /**
   * @brief Destroy the Lagrange Polynomial object
   *
   */
  ~LagrangePolynomial() { delete tpp; }
  /**
   * @brief Get the first derivative of the Lagrange Polynomial
   *
   * @param i The (number of the) shape function
   * @param q The quadrature point at which the function is evaluated
   * @return Tensor<1, dim> The derivative with respect to generalized
   * coordinates
   */
  Tensor<1, dim> get_first_derivative(int i, int q);
  /**
   * @brief Get the value object
   *
   * @param i The (number of the) shape function
   * @param q The quadrature point at which the function is evaluated
   * @return double The function value
   *
   * Provided for the sake of completeness. In general
   * fe(face)values::shape_value should be used.
   */
  double get_value(int i, int q);

private:
  /**
   * @brief Assemble the dim-dimensional LagrangePolynomial from 1-dimensional
   * Polynomials
   *
   * Called from the constructor.
   *
   */
  void init(int poly_degree);
  /**
   * @brief Internal Gauss-Quadrature-object.
   *
   * Used to get evaluation points
   *
   */
  QGauss<dim> quadrature;
  /**
   * @brief Pointer to the TensorProduct of Lagrange Functions
   *
   * Usage of pointer is required as the arguments of it's constructor are
   * generated after the object is created.
   *
   */
  TensorProductPolynomials<dim, Polynomials::Polynomial<double>> *tpp;
};

template <int dim> void LagrangePolynomial<dim>::init(int poly_degree) {

  std::vector<Point<1>> support_points;
  if (poly_degree == 1) {
    support_points = std::vector<Point<1>>{Point<1>(0), Point<1>(1)};
  } else if (poly_degree == 2) {
    support_points = std::vector<Point<1>>{Point<1>(0), Point<1>(0.5), Point<1>(1)};
  } else {
    delete tpp;
    cexc::not_imp_error exc;
    BOOST_THROW_EXCEPTION(exc);
  }
  std::vector<Polynomials::Polynomial<double>> pol;
  for (unsigned int i=0; i<support_points.size(); ++i) {
    pol.push_back(Polynomials::Polynomial<double>(support_points, i));
  }
  tpp = new TensorProductPolynomials<dim, Polynomials::Polynomial<double>>(pol);
}

template <int dim>
Tensor<1, dim> LagrangePolynomial<dim>::get_first_derivative(int i, int q) {

  Point<dim> q_point = quadrature.point(q);
  return tpp->compute_1st_derivative(i, q_point);
};
template <int dim> double LagrangePolynomial<dim>::get_value(int i, int q) {

  Point<dim> q_point = quadrature.point(q);
  return tpp->compute_value(i, q_point);
};
} // namespace fem

#endif