/**
 * @file boundaryvalues.h
 * @author Till Budde (tilljanis.budde@tu-dortmund.de)
 * @brief Implements Boundary Values for the FEM
 * @version 0.1
 * @date 2021-06-28
 *
 * @copyright Copyright (c) 2021
 *
 */
#ifndef BOUNDARVALUES_H
#define BOUNDARVALUES_H

#include "CustomExceptions.h"
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

namespace fem {

using namespace dealii;

/**
 * @brief Class to calculate boundary values
 *
 * @tparam dim number of dimensions
 *
 * In order to assemble the AffineConstraints object the
 * VectorTools::interpolate_boundary_values() function requires an object of
 * type Function that implements the vector_value_list function to calculate
 * boundary values at each quadrature point on the face
 */
template <int dim> class BoundaryValues : public Function<dim> {
public:
  /**
   * @brief Construct a new Boundary Values object
   *
   * @param int The current timestep
   * @param int The total number of timesteps
   * @param double The total displacement
   *
   * Initializes the current displacement on the faces under the assumption of
   * a linear loadcurve
   */
  BoundaryValues(unsigned int, unsigned int, double);
  /**
   * @brief Calculate the Function value at a point
   *
   * @param p  The point at which the function is evaluated
   * @param values Output: The value of the function at that point
   *
   * This function sets the vertical component of the value vector to the
   * \ref current_displacement.
   * The size of the Point \ref p and the \ref values Vector have to be equal.
   */
  virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;

  /**
   * @brief Call vector_value to assign each quadrature point its boundary value
   *
   * @param points  The vector of points at which the boundary condition is
   * evaluated
   * @param value_list Output: The values to which the boundary condition is set
   *
   * This function loops over all points and sets the entries of value_list
   * to the value of vector_value.
   *
   * The size of the \ref points vector and the \ref value_list have to be
   * equal.
   */
  virtual void vector_value_list(const std::vector<Point<dim>> &points,
                                 std::vector<Vector<double>> &value_list) const;

private:
  double total_displacement;
  double displacement_step;
  double current_displacement;
};

// Boundaryvalues Functions --------------------------------------------

template <int dim>
BoundaryValues<dim>::BoundaryValues(unsigned int timestep,
                                    unsigned int no_timesteps,
                                    double total_displacement)
    : Function<dim>(dim) {
  displacement_step = total_displacement / no_timesteps;
  current_displacement = displacement_step * timestep;
}

template <int dim>
void BoundaryValues<dim>::vector_value(const Point<dim> & /*p*/,
                                       Vector<double> &values) const {
  // Check wether the dimension of both vectors is equal
  Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
  values = 0;
  // set the horizontal component to the current_displacement
  if (dim > 1) {
    values(1) = displacement_step;
  } else {
    cexc::exception_base exc;
    BOOST_THROW_EXCEPTION(exc);
  }
}

template <int dim>
void BoundaryValues<dim>::vector_value_list(
    const std::vector<Point<dim>> &points,
    std::vector<Vector<double>> &value_list) const {

  const unsigned int n_points = points.size();
  // Check that for each point, at which the function should be evaluated,
  // an entry for the function value exists
  Assert(value_list.size() == n_points,
         ExcDimensionMismatch(value_list.size(), n_points));

  // loop over all points
  for (unsigned int p = 0; p < n_points; ++p) {
    BoundaryValues<dim>::vector_value(points[p], value_list[p]);
  }
}

// Boundaryvalues Functions End ----------------------------------------

} // namespace fem

#endif