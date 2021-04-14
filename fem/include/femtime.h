#ifndef FEMTIME_H
#define FEMTIME_H

namespace fem {
    
	// Time manages all time related variables (and separates them from 
	// outside influences)
	class Time
	{
		public:
		Time ()
		:
		timestep(0),
		time_current(0.0),
		time_end(10),
		dt(0.5)
		{
			double tmp =  time_end/dt;
			no_timesteps = (unsigned int) tmp;
		}
		virtual ~Time()
		{}
		double get_current() const
		{
			return time_current;
		}
		double get_end() const
		{
			return time_end;
		}	  
		double get_dt() const
		{
			return dt;
		}
		unsigned int get_timestep() const
		{
			return timestep;
		}
		unsigned int get_no_timesteps() const
		{
			return no_timesteps;
		}
		void increment()
		{
			time_current += dt;
			++timestep;
		}
		private:
		unsigned int timestep;
		double time_current;
		double time_end;
		double dt;
		unsigned int no_timesteps;
	};	
}

#endif