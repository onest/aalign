/*
* (c) 2015 Virginia Polytechnic Institute & State University (Virginia Tech)   
*                                                                              
*   This program is free software: you can redistribute it and/or modify       
*   it under the terms of the GNU General Public License as published by       
*   the Free Software Foundation, version 2.1                                  
*                                                                              
*   This program is distributed in the hope that it will be useful,            
*   but WITHOUT ANY WARRANTY; without even the implied warranty of             
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              
*   GNU General Public License, version 2.1, for more details.                 
*                                                                              
*   You should have received a copy of the GNU General Public License          
*                                                                              
*  Acknowledgement                                                             
*   This research was supported in part by NSF-BIGDATA program via             
*   IIS-1247693                                                                
*                                                                              
*/



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
