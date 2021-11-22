/**
 * @file timer.h
 * @author ??? Ripped from stackoverflow
 * @brief Implement Class to simplify measuring time
 * @version 0.1
 * @date 2021-11-22
 * 
 * 
 */
#ifndef TIMER_H
#define TIMER_H

#include <chrono> // for std::chrono functions


/**
 * @brief Class to simplify measuring time
 * 
 */
class Timer
{
private:
	// Type aliases to make accessing nested type easier
	using clock_type = std::chrono::steady_clock;
	using second_type = std::chrono::duration<double, std::ratio<1> >;
	
	std::chrono::time_point<clock_type> m_beg;

public:
	/**
	 * @brief Construct a new Timer object
	 * 
	 */
	Timer() : m_beg { clock_type::now() }
	{
	}
	
	/**
	 * @brief Reset the clock
	 * 
	 */
	void reset()
	{
		m_beg = clock_type::now();
	}
	
	/**
	 * @brief get elapsed time since last reset or construction
	 * 
	 * @return double elapsed time in seconds
	 */
	double elapsed() const
	{
		return std::chrono::duration_cast<second_type>(clock_type::now() - m_beg).count();
	}
};
#endif