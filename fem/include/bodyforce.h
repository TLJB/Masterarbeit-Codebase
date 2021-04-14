#ifndef BODYFORCE_H
#define BODYFORCE_H

#include <deal.II/lac/vector.h>
#include <deal.II/base/point.h>
#include <deal.II/base/function.h>

namespace fem {

    using namespace dealii;

	// The BodyForce class handles bodyforces		
	template <int dim>
	class BodyForce :  public Function<dim>
	{
		public:
		BodyForce ();
		virtual
		void
		vector_value (const Point<dim> &p,
						Vector<double>   &values) const;
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
	// body_force.vector_value calculates  the bodyforces	
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
	// vector_value at each quadrature_point	
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