/**
 * @file femtime.h
 * @author Till Budde (tilljanis.budde@tu-dortmund.de)
 * @brief Implements a "Time" object
 * @version 0.1
 * @date 2021-06-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef FEMTIME_H
#define FEMTIME_H

namespace fem {
    
	/**
	 * @brief Class to manage all "time" related issues. 
	 * 
	 */
	class Time
	{
		public:
		/**
		 * @brief Construct a new Time object
		 * 
		 */
		Time ()
		:
		timestep(0),
		time_current(0.0),
		time_end(1.),
		dt(.25)
		{
			double tmp =  time_end/dt;
			no_timesteps = (unsigned int) tmp;
		}
		/**
		 * @brief Destroy the Time object
		 * 
		 */
		virtual ~Time()
		{}
		/**
		 * @brief Get the current time
		 * 
		 * @return double current time
		 */
		double get_current() const
		{
			return time_current;
		}
		/**
		 * @brief Get the maximum time
		 * 
		 * @return double end time
		 */
		double get_end() const
		{
			return time_end;
		}	  
		/**
		 * @brief Get the size of the timestep
		 * 
		 * @return double dt
		 */
		double get_dt() const
		{
			return dt;
		}
		/**
		 * @brief Get the number of the current timestep
		 * 
		 * @return unsigned int timestep
		 */
		unsigned int get_timestep() const
		{
			return timestep;
		}
		/**
		 * @brief Get the number of timesteps
		 * 
		 * @return unsigned int number of timesteps
		 */
		unsigned int get_no_timesteps() const
		{
			return no_timesteps;
		}
		/**
		 * @brief increase time
		 * 
		 */
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