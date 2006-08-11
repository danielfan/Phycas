#ifndef PHO_STOPWATCH_H
#define PHO_STOPWATCH_H

#if defined(POL_NEW_STOPWATCH)
#include "phycas/misc/stopwatch_new.hpp"
#else

#include <ctime>

/*----------------------------------------------------------------------------------------------------------------------
|	Implements a very basic stopwatch. Call the Start member function to start or restart the stopwatch, and ElapsedTime
|	to report the number of seconds since Start was last called. Note that ElapsedTime does not stop the watch.
*/
class StopWatch
	{
	public:
				StopWatch();
				~StopWatch();
		
		void	Start();
		void	Stop();

		bool	IsRunning() const;
		double	ElapsedTime() const;
		std::clock_t	ElapsedTicks() const;
		std::clock_t	TicksPerSecond() const;

	private:
		bool	running;
		std::clock_t started;
		std::clock_t	stopped;
	};

inline StopWatch::StopWatch()
	{
	running	= false;
	started	= 0;
	stopped	= 0;
	}

inline StopWatch::~StopWatch()
	{
	}

inline bool	StopWatch::IsRunning() const
	{
	return running;
	}

inline void StopWatch::Start()
	{
	if (running)
		Stop();
	running = true;
	started = std::clock();
	}

inline void StopWatch::Stop()
	{
	if (running) 
		stopped = std::clock();
	running = false;
	}
		
inline double StopWatch::ElapsedTime() const
	{
	if (running)
		return ((double)(std::clock() - started)/(double)CLOCKS_PER_SEC);
	else
		return ((double)(stopped - started)/(double)CLOCKS_PER_SEC);
	}

inline std::clock_t StopWatch::ElapsedTicks() const
	{
	if (running)
		return (std::clock() - started);
	else
		return (stopped - started);
	}

inline std::clock_t StopWatch::TicksPerSecond() const
	{
	return CLOCKS_PER_SEC;
	}

#endif

#endif	// #if defined(POL_NEW_STOPWATCH)

