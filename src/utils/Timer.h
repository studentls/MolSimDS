//------------------------------------------------------------------------------------------------
// (C) 2012 by F.Dietz & L.Spiegelberg
// include License here...

//At the moment you must not distribute or use this code at all,
//cause license details are not clear yet
//If you want to receive updates feel free to subscribe
//------------------------------------------------------------------------------------------------
// File Timer.h
// implements a high accurate timer class on multiple platforms
//------------------------------------------------------------------------------------------------

#ifndef TIMER_HEADER_
#define TIMER_HEADER_

#include "Base.h"

//windows
#ifdef WINDOWS
#include <Windows.h>
#endif

//linux
#ifdef LINUX
#include <sys/time.h>
#endif

// mac os x
#ifdef MACINTOSH
#include <sys/time.h>
#endif

namespace utils
{
	/// defines a timer class, used to accurately measure time with high accurate timers
	class Timer
	{
	private:
			
#ifdef WINDOWS
		// last saved timestamp
		LARGE_INTEGER 	lLastTime;
#endif

#ifdef LINUX
		// last saved timestamp
		timespec	lLastTime;
#endif
#ifdef MACINTOSH
        // last timestamp
        timeval lLastTime;
#endif

	public:

		Timer()
		{
            reset();
		}

		/// reset timer
		void	reset()
		{
#ifdef WINDOWS
		QueryPerformanceCounter(&lLastTime);
#endif
#ifdef LINUX
		//use CLOCK_PROCESS_CPUTIME_ID to measure CPU Performance, not IO
		clock_gettime(CLOCK_MONOTONIC, &lLastTime);
#endif
#ifdef MACINTOSH
            gettimeofday(&lLastTime, NULL);
#endif
		}

		/// returns elapsed time since last call of constructor or reset() in seconds
		double	getElapsedTime(const bool reset = false)
		{
#ifdef WINDOWS
			LARGE_INTEGER lTime;
			LARGE_INTEGER lFrequency;

			QueryPerformanceCounter(&lTime);
			QueryPerformanceFrequency(&lFrequency);

			lTime.QuadPart -= lLastTime.QuadPart;

			double elapsed = (double)lTime.QuadPart / (double)lFrequency.QuadPart;

			if(reset)this->reset();

			return  elapsed;
#endif

#ifdef LINUX
			//see http://stackoverflow.com/questions/538609/high-resolution-timer-with-c-and-linux
			timespec lTime;
			timespec lDiff;

			clock_gettime(CLOCK_MONOTONIC, &lTime);
			
			lDiff.tv_sec = lTime.tv_sec - lLastTime.tv_sec;
			lDiff.tv_nsec = lTime.tv_nsec - lLastTime.tv_nsec;

			//is there a carry?
			if(lDiff.tv_nsec < 0)
			{
				lDiff.tv_sec -= 1;
				lDiff.tv_nsec += 1000000000; //10^9
			}
			
			double elapsed = (double)lDiff.tv_sec + (double)lDiff.tv_nsec / 1000000000.0;
			
			if(reset)this->reset();
			
			return elapsed;
#endif
#ifdef MACINTOSH
            timeval lTime, lDiff;
            gettimeofday(&lTime, NULL);
            
            lDiff.tv_sec = lTime.tv_sec - lLastTime.tv_sec;
			lDiff.tv_usec = lTime.tv_usec - lLastTime.tv_usec;
            
			//is there a carry?
			if(lDiff.tv_usec < 0)
			{
				lDiff.tv_sec -= 1;
				lDiff.tv_usec += 1000000000; //10^9
			}
			
			double elapsed = (double)lDiff.tv_sec + (double)lDiff.tv_usec / 1000000000.0;
            
            if(reset)this->reset();
            
            return elapsed;
#endif
			return 0.0;
		}

	};
}



#endif
