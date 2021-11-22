/**
 * @file pointhistory.h
 * @author Till Budde
 * @brief Implements structure to save state dependent variables at quadrature
 * points
 * @version 0.1
 * @date 2021-06-28
 *
 * @copyright Copyright (c) 2021
 *
 */
#ifndef POINHISTORY_H
#define POINHISTORY_H

#include <deal.II/base/symmetric_tensor.h>

namespace fem {

using namespace dealii;

/**
 * @brief State dependent variables of bulk material
 *
 * @tparam dim The dimension
 *
 * Members are state dependent variables and other values (such as stress)
 * That can be saved at each quadrature point
 *
 */
template <int dim> struct PointHistoryBulk {
  SymmetricTensor<2, dim> strain_pl;
  double alpha;
  SymmetricTensor<2, dim> old_stress;
};

/**
 * @brief State dependent variables of interface material
 *
 * @tparam dim The dimension
 *
 * Members are state dependent variables and other values (such as traction)
 * That can be saved at each quadrature point
 *
 */
template <int dim> struct PointHistoryInter {
  double kappa = 0.005;
  double alpha = 1e9;
  double alpha0 = alpha;
  bool pen = false;
};
/**
 * @brief Contains structs to save State dependent variable for interface and
 * 				bulk material
 *
 * @tparam dim  number of dimensions
 *
 * State dependent variables of interface or bulk can be accessed through
 * members \ref bulk and \ref inter
 *
 */
template <int dim> struct PointHistory {
  /**
   * @brief State dependent variables of bulk material
   *
   */
  PointHistoryBulk<dim> bulk;
  /**
   * @brief State dependent variables of interface material
   *
   */
  PointHistoryInter<dim> inter;
};

} // namespace fem

#endif