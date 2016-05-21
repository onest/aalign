/*
* (c) 2015 Virginia Polytechnic Institute & State University (Virginia Tech)   
*                                                                              
*   This program is free software: you can redistribute it and/or modify       
*   it under the terms of the GNU General Public License as published by       
*   the Free Software Foundation, either version 3 of the License, or          
*   (at your option) any later version.                                        
*                                                                              
*   This program is distributed in the hope that it will be useful,            
*   but WITHOUT ANY WARRANTY; without even the implied warranty of             
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              
*   GNU General Public License for more details.                               
*                                                                              
*   You should have received a copy of the GNU General Public License          
*                                                                              
*  Acknowledgement                                                             
*   This research was supported in part by NSF-BIGDATA program via             
*   IIS-1247693                                                                
*                                                                              
*/


#ifndef _NEEDLEMAN_WUNSCH_CPU_H
#define _NEEDLEMAN_WUNSCH_CPU_H

int sequential_NW_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_NW_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_NW_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_NW_CPU(const char* queryStr, int m, const char* subjectStr, int n);

#endif
