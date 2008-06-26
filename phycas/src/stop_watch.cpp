/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <cmath>
#include <limits>
#include "phycas/src/stop_watch.hpp"    // includes <ctime>

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets ...
*/
StopWatch::StopWatch() 
	{
    reset();
    max_ticks = std::numeric_limits<std::clock_t>::max();
    min_ticks = std::numeric_limits<std::clock_t>::min();
    ticks_per_sec = (double)CLOCKS_PER_SEC;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor clears ...
*/
StopWatch::~StopWatch() 
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|   
*/
void StopWatch::reset()
	{
    running = false;
    elapsed_seconds = 0.0;
    start_ticks = (clock_t)0;
    stop_ticks = (clock_t)0;
    start_time = (time_t)0;
    stop_time = (time_t)0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   
*/
void StopWatch::start()
	{
    start_time = time(NULL);
    stop_time = start_time;
    start_ticks = clock();
    stop_ticks = start_ticks;
    elapsed_seconds = 0.0;
    running = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Stops the StopWatch object and updates `elapsed_seconds'.
*/
void StopWatch::stop()
	{
    stop_time = time(NULL);
    stop_ticks = clock();
    normalize();
    running = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns number of seconds that have elapsed since start() was called. Does not call stop().
*/
double StopWatch::elapsedSeconds()
	{
    normalize();
    return elapsed_seconds;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Updates `elapsed_seconds' data member by computing the difference between a clock() call made now and the one made
|   when timing was started. Overflow in the result from clock() is taken into account.
*/
void StopWatch::normalize()
	{
    if (running)
        {
        stop_ticks = clock();
        if (stop_ticks < 0 && start_ticks > 0)
            {
            double a = (double)(max_ticks - start_ticks);
            double b = (double)(stop_ticks - min_ticks);
            elapsed_seconds += (a + b + 1)/ticks_per_sec;
            }
        else
            {
            elapsed_seconds += (double)(stop_ticks - start_ticks)/ticks_per_sec;
            }
        
        // sanity check
        stop_time = time(NULL);
        double crude_elapsed_seconds = (double)(stop_time - start_time);
        if (crude_elapsed_seconds > 0.0)
            {
            double diff = std::fabs(elapsed_seconds - crude_elapsed_seconds);
            //PHYCAS_ASSERT(diff < 1.0);   // allow up to 1 sec. difference between clock() and time(NULL) estimates
            }

        start_ticks = stop_ticks;    
        }
	}

long StopWatch::stopTicks()
	{
    return stop_ticks;
	}

std::string StopWatch::tobinary(long x)
	{
    std::string s;
    long unity = 1L;
    unsigned nbytes = (unsigned)sizeof(long);
    unsigned nbits = 8*nbytes;
    for (unsigned i = 0; i < nbits; ++i)
        {
        if (x & unity)
            s.push_back('1');
        else
            s.push_back('0');
        x >>= 1;
        }
    std::string revs;
    for (std::string::reverse_iterator it = s.rbegin(); it != s.rend(); ++it)
        {
        revs.push_back(*it);
        }
    return revs;
    }

#if 0
void StopWatch::doofus()
	{
#if 1
    std::string s;
    s = tobinary(LONG_MIN);
    std::cerr << "LONG_MIN + 0L = " << (LONG_MIN) << " = " << s << std::endl;
    s = tobinary(LONG_MIN + 1L);
    std::cerr << "LONG_MIN + 1L = " << (LONG_MIN + 1L) << " = " << s << std::endl;
    s = tobinary(LONG_MIN + 2L);
    std::cerr << "LONG_MIN + 2L = " << (LONG_MIN + 2L) << " = " << s << std::endl;
    s = tobinary(-1L);
    std::cerr << "-1L = " << (-1L) << " = " << s << std::endl;
    double sz;
    sz = (double)(std::numeric_limits<std::clock_t>::max());
    std::cerr << "std::numeric_limits<clock_t>::max() = " << sz << std::endl;
    sz = (double)(std::numeric_limits<std::clock_t>::min());
    std::cerr << "std::numeric_limits<clock_t>::min() = " << sz << std::endl;
#endif

#if 1
    long clocks_per_sec = 1000000;
    long x = 0L;
    double z = 0.0;
    for (unsigned i = 0; i < 2148; ++i)
        {
        long prev_x = x;
        x += clocks_per_sec;
        if (x < 0 && prev_x > 0)
            {
            long a = LONG_MAX - prev_x;
            long b = x - LONG_MIN;
            z += (double)(a + b + 1);
            }
        else
            {
            z += (double)(x - prev_x);
            }
        std::cerr << i << '\t' << prev_x << '\t' << x << '\t' << z << std::endl;
        }
#endif
	}
#endif

}	// namespace phycas
