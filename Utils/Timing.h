#ifndef _TIMING_H
#define _TIMING_H

#include <sys/time.h>

/**
 * This method returns the current timestamp.
 *
 * @return the current time in second(s)
 */
double dtime() 
{
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime, (struct timezone*)0);
    // tseconds = (double)(mytime.tv_sec + mytime.tv_usec*1.0e-6);
    tseconds = (double)(mytime.tv_sec*1.0e3 + mytime.tv_usec*1.0e-3);
    return (tseconds);
}

#endif
