#ifndef BOUNDARVALUES_H
#define BOUNDARVALUES_H

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>


namespace fem {

    using namespace dealii;

	template<int dim>
	class BoundaryValues : public Function<dim>
	{	
		public:
		BoundaryValues(unsigned int, unsigned int, double);
		virtual void vector_value (const Point<dim> &p,
									Vector<double>  &values) const;
							
		virtual void vector_value_list (const std::vector<Point<dim> > &points,
										std::vector<Vector<double> > &value_list) const;
		private:
		double total_displacement;
		double displacement_step;
		double current_displacement;
	};  


// Boundaryvalues Functions --------------------------------------------

	template <int dim>
	BoundaryValues<dim>::BoundaryValues(unsigned int timestep, unsigned int no_timesteps, double total_displacement)
	:
	Function<dim> (dim)
	{
	displacement_step = total_displacement/no_timesteps;
	current_displacement = displacement_step*timestep;
	}  

	// The Boundary value assigned is an incremental displacement in 
	// vertical direction
	template<int dim>
	void BoundaryValues<dim>:: vector_value (const Point<dim> &/*p*/,
											Vector<double>   &values) const
	{
		Assert (values.size() == dim,
				ExcDimensionMismatch (values.size(), dim));
		values = 0;
		values(dim-1) = current_displacement;
	}

	// vector_value_list calls vector_value at each quadrature_point
	template <int dim>
	void BoundaryValues<dim>::vector_value_list (const std::vector<Point<dim> > &points,
											 std::vector<Vector<double> >   &value_list) const
	{
		const unsigned int n_points = points.size();
	
		Assert (value_list.size() == n_points,
				ExcDimensionMismatch (value_list.size(), n_points));
	
		for (unsigned int p=0; p < n_points; ++p)
		{
		BoundaryValues<dim>::vector_value (points[p], value_list[p]);
		}
	}

// Boundaryvalues Functions End ----------------------------------------

}

#endif