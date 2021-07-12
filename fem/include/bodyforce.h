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

#include <deal.II/lac/vector.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>

namespace fem {

  using namespace dealii;

	/**
	 * @brief Class to apply body forces
	 * 
	 * @tparam dim number of dimensions
	 */
	template <int dim>
	class BodyForce :  public Function<dim>
	{
		public:
		/**
		 * @brief Construct a new Body Force object
		 * 
		 */
		BodyForce ();
		/**
		 * @brief calculates the body force at a  point
		 * 
		 * @param p Point
		 * @param values returns Bodyforce values at point 
		 */
		virtual
		void
		vector_value (const Point<dim> &p,
						Vector<double>   &values) const;
		/**
		 * @brief calculates the body force for vector of points
		 * 
		 * @param points List of points
		 * @param value_list returns Bodyforces at points
		 */
		virtual
		void
		vector_value_list (const std::vector<Point<dim> > &points,
							std::vector<Vector<double> >   &value_list) const;
	};


// Bodyforce Functions -------------------------------------------------
 
	template <int dim>
	BodyForce<dim>::BodyForce ()
	:
	Function<dim> (dim)
	{}
	template <int dim>
	inline
	void
	BodyForce<dim>::vector_value (const Point<dim> &/*p*/,
									Vector<double>   &values) const
	{
		Assert (values.size() == dim,
				ExcDimensionMismatch (values.size(), dim));
		const double g   = 0;	
		const double rho = 7850;
		values = 0;
		values(dim-1) = -rho * g;
	}
	
	// vector_value_list ditributes the bodyforces by calling 
	// vector_value at each point	
	template <int dim>
	void
	BodyForce<dim>::vector_value_list (const std::vector<Point<dim> > &points,
										std::vector<Vector<double> >   &value_list) const
	{
		const unsigned int n_points = points.size();
		Assert (value_list.size() == n_points,
				ExcDimensionMismatch (value_list.size(), n_points));
		for (unsigned int p=0; p<n_points; ++p)
		BodyForce<dim>::vector_value (points[p],
										value_list[p]);
	}

// Bodyforce Functions End ---------------------------------------------
}


#endif