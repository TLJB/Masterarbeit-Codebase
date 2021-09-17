/**
 * @file bodyforce.h
 * @author Till Budde (tilljanis.budde@tu-dortmund.de)
 * @brief Implements class to apply body forces
 * @version 0.1
 * @date 2021-06-28
 *
 * @copyright Copyright (c) 2021
 *
 */
#ifndef BODYFORCE_H
#define BODYFORCE_H

#include "CustomExceptions.h"
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

namespace fem {

using namespace dealii;

/**
 * @brief Class to calculate body forces
 *
 * @tparam dim number of dimensions
 */
template <int dim> class BodyForce : public Function<dim> {
public:
  /**
   * @brief Construct a new Body Force object
   *
   */
  BodyForce();
  /**
   * @brief calculates the body force at a  point
   *
   * @param p  The point at which the function is evaluated
   * @param values Output: The value of the function at that point
   *
   * This function sets the vertical component of the value vector to the
   * body force
   * The size of the Point \ref p and the \ref values Vector have to be equal.
   */
  virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;
  /**
   * @brief calculates the body force for vector of points
   *
   * @param points  The vector of points at which the body force is
   * evaluated
   * @param value_list Output: The values to which the body force is set
   *
   * This function loops over all points and sets the entries of value_list
   * to the value of vector_value.
   *
   * The size of the \ref points vector and the \ref value_list have to be
   * equal.
   */
  virtual void vector_value_list(const std::vector<Point<dim>> &points,
                                 std::vector<Vector<double>> &value_list) const;

  /**
   * @brief gravitational acceleration
   *
   */
  double g;
  /**
   * @brief density
   *
   */
  double rho;
};

// Bodyforce Functions -------------------------------------------------

template <int dim> BodyForce<dim>::BodyForce() : Function<dim>(dim) {
  g = 0;
  rho = 7850;
}
template <int dim>
inline void BodyForce<dim>::vector_value(const Point<dim> & /*p*/,
                                         Vector<double> &values) const {
  // Check wether the dimension of both vectors is equal
  Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
  values = 0;
  // set the horizontal component to the current_displacement
  if (dim > 1) {
    values(1) = -rho * g;
  } else {
    cexc::exception_base exc;
    BOOST_THROW_EXCEPTION(exc);
  }
}

template <int dim>
void BodyForce<dim>::vector_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>> &value_list) const {

  const unsigned int n_points = points.size();
  // Check that for each point, at which the function should be evaluated,
  // an entry for the function value exists
  Assert(value_list.size() == n_points,
         ExcDimensionMismatch(value_list.size(), n_points));

  // loop over all points
  for (unsigned int p = 0; p < n_points; ++p)
    BodyForce<dim>::vector_value(points[p], value_list[p]);
}

// Bodyforce Functions End ---------------------------------------------

} // namespace fem

#endif