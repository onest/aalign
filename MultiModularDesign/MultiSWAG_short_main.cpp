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



#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include <climits>
#include "Config.h"
#include "Timing.h"
using namespace std;
void align_db(const char *queryStr, int queryNum, const char *dbStrs, int *dbOffset, int dbNum, int16_t *output, int numThread);
int main(int argc, char** argv)
{
    int opt;
    char *queryFile;
    char *dbFile;
    int numThreads = 1;
    while((opt = getopt(argc, argv, "q:d:t:")) != -1) 
    {
        switch(opt) 
        {
            case 'q':
                queryFile = optarg; // contains only one sequence
                break;
            case 'd':
                dbFile = optarg; // contains a database of sequences
                break;
            case 't':
                numThreads = atoi(optarg); // contains a database of sequences
                break;
            default:
                std::cerr << "Usage" << argv[0] << " [-q query file (-)] [-d database file (-)] [-t thread num (1)]" << std::endl;
                exit(EXIT_FAILURE);
        }
    }

    string queryStr;
    char *dbStrs;
    int *dbOffset;
    int dbNum;
    queryStr = read_file(queryFile);
    vector<size_t> ind;
    vector<string> name = read_db_file(dbFile, dbStrs, dbOffset, dbNum, ind);
    cout << "loading finished" << endl;

    int16_t *output = (int16_t *)malloc(sizeof(int16_t)*dbNum);
    // cout << "query:" << queryStr << endl;
    // cout << "dbNum:" << dbNum << endl;
    // for(int i = 0; i < dbNum; i++)
    // {
        // string db(dbStrs+dbOffset[i], dbOffset[i+1]-dbOffset[i]);
        // cout << db << endl;
    // }
    double tstart, tstop, ttime;
    tstart = dtime();
    align_db(queryStr.c_str(), queryStr.length(), dbStrs, dbOffset, dbNum, output, numThreads);
    tstop  = dtime();
    ttime  = tstop - tstart;
    double gcups = ((double)dbOffset[dbNum]) * queryStr.length();
    gcups /= 1000000 * (ttime);
    printf("Runtime: %f ms %f GCUPS\n", ttime, gcups);

    for(int i = 0; i < 10; i++)
    {
        int16_t *it = max_element(output, output+dbNum);
        if(*it != SHRT_MIN)
        {
            cout << "score " << *it << ": " << name[ind[it - output]] << endl;
            output[it-output] = SHRT_MIN;
        }
    }
    // for(int i =0; i < dbNum; i++)
    // {
        // cout << output[i] << endl;
    // }
    free(output);
    free(dbStrs);
    free(dbOffset);
}


