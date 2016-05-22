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



#ifndef _OPT_STR_ALIGN_DISTANCE_MIC_H
#define _OPT_STR_ALIGN_DISTANCE_MIC_H

int sequential_OSAD_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_OSAD_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_OSAD_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_OSAD_MIC(const char* queryStr, int m, const char* subjectStr, int n);

#endif
