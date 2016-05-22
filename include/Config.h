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


#ifndef _CONFIG_H
#define _CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#define LD_INDEL 1
#define LD_SUBSTITUTION 2 
#define OSAD_INDEL 1
#define OSAD_SUBSTITUTION 1 
#define OSAD_TRANSPOSITION 1 
#define GAP_OPEN_PENALTY  -12
// #define GAP_OPEN_PENALTY  -10
#define GAP_EXT_PENALTY   -2
// #define GAP_EXT_PENALTY   -1
#define GAP_PENALTY    -4
#define MATCH_SCORE     5 
#define MISMATCH_SCORE -3 
#define AVX3_WIDTH 16
#define AVX2_WIDTH 8 
#define MAX_AMINO_ACIDS 23 
#define SUBS_MATRIX_SIZE (MAX_AMINO_ACIDS + 1)
#define LD_MATRIX_SIZE 26 
#define OSAD_MATRIX_SIZE 26 

#define USE_BLOSUM
// A convenient way to access BLOSUM62 matrix
#define BLOSUM62(i, j) blosum62[(i)*(32)+(j)]
static const int8_t blosum62[32*32] = {
 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  4, -3,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4, -3,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0, -2,  4, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1, -3,  4, -1,  0,  0,  0,  0,  0,  0,  0,  0,
 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0, -3,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3,  3, -3, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4,  3, -3, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0, -3,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3,  2, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3,  0, -3, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -3, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,
 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -2, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -1, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,
 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3,  2, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-2, -1,  4,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4, -3,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -3,  2,  0, -3, -2, -1, -2, -1,  2, -3,  3, -3, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -2, -2, -2,  0, -3,  4, -1,  0,  0,  0,  0,  0,  0,  0,  0,
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
};

int blosum62char2int(char base);
char blosum62int2char(int index);
void striped_vec_scan_max(int *, int *, int , int , int);
void striped_vec_scan_min(int *, int *, int , int , int);
std::string read_file(char *);
std::vector<std::string> read_db_file(char *, char *&, int *&, int &, std::vector<size_t>&);

#endif
