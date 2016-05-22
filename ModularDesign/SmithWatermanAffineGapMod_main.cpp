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
#include <string>
#include <immintrin.h> 
#include <getopt.h>
#include "Config.h"
#include "Timing.h"
using namespace std;

int mod_striped_iterate_SWAG(const char* queryStr, int m, const char* subjectStr, int n);
int mod_striped_scan_SWAG(const char* queryStr, int m, const char* subjectStr, int n);
int mod_striped_merged_SWAG(const char* queryStr, int m, const char* subjectStr, int n);
int seq_SWAG(const char* queryStr, int m, const char* subjectStr, int n);
int seq_scan_SWAG(const char* queryStr, int m, const char* subjectStr, int n);

int main(int argc, char **argv)
{
    string queryStr = //"GACTTAC";
        "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAA"
        "PRVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKT"
        "CPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRN"
        "TFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACAGR"
        "DRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALEL"
        "KDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD";
    string subjectStr = //"CGTGAATTCAT";
        "MTAMEESQSDISLELPLSQETFSGLWKLLPPEDILPSPHCMDDLLLPQDVEEFFEGPSEA"
        "LRVSGAPAAQDPVTETPGPAAPAPATPWPLSSFVPSQKTYQGNYGFHLGFLQSGTAKSVM"
        "CTYSPPLNKLFCQLAKTCPVQLWVSATPPAGSRVRAMAIYKKSQHMTEVVRRCPHHERCS"
        "DGDGLAPPQHLIRVEGNLYPEYLEDRQTFRHSVVVPYEPPEAGSEYTTIHYKYMCNSSCM"
        "GGMNRRPILTIITLEDSSGNLLGRDSFEVRVCACPGRDRRTEEENFRKKEVLCPELPPGS"
        "AKRALPTCTSASPPQKKKPLDGEYFTLKIRGRKRFEMFRELNEALELKDAHATEESGDSR"
        "AHSSYLKTKKGQSTSRHKKTMVKKVGPDSD";

    int opt;
    char* _queryFile;
    char* _dbFile;
    while((opt = getopt(argc, argv, "q:d:")) != -1) 
    {
        switch(opt) 
        {
            case 'q':
                _queryFile = optarg; 
                break;
            case 'd':
                _dbFile= optarg; 
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }
    queryStr = read_file(_queryFile);
    subjectStr = read_file(_dbFile);

    double tstart, tstop, ttime;
    tstart = dtime();
    int ret_seq = seq_SWAG(queryStr.c_str(), queryStr.length(), subjectStr.c_str(), subjectStr.length());
    tstop  = dtime();
    ttime  = tstop - tstart;
    printf("time\tseq_orig\t%lf\t%d\n", ttime, ret_seq);

    tstart = dtime();
    int ret_sca = seq_scan_SWAG(queryStr.c_str(), queryStr.length(), subjectStr.c_str(), subjectStr.length());
    tstop  = dtime();
    ttime  = tstop - tstart;
    printf("time\tseq_scan\t%lf\t%d\n", ttime, ret_sca);

    tstart = dtime();
    int ret_vec_mod = mod_striped_scan_SWAG(queryStr.c_str(), queryStr.length(), subjectStr.c_str(), subjectStr.length());
    tstop  = dtime();
    ttime  = tstop - tstart;
    printf("time\tvec_scan\t%lf\t%d\n", ttime, ret_vec_mod);
    // printf("time\tvec_scan(avg)\t%lf\t%d\n", ttime/subjectStr.length(), ret_vec_mod);
    tstart = dtime();
    int ret_vec_mod1 = mod_striped_iterate_SWAG(queryStr.c_str(), queryStr.length(), subjectStr.c_str(), subjectStr.length());
    tstop  = dtime();
    ttime  = tstop - tstart;
    printf("time\tvec_iter\t%lf\t%d\n", ttime, ret_vec_mod1);
    // printf("time\tvec_iter(avg)\t%lf\t%d\n", ttime/subjectStr.length(), ret_vec_mod1);

    tstart = dtime();
    int ret_vec_mod2 = mod_striped_merged_SWAG(queryStr.c_str(), queryStr.length(), subjectStr.c_str(), subjectStr.length());
    tstop  = dtime();
    ttime  = tstop - tstart;
    printf("time\tvec_merg\t%lf\t%d\n", ttime, ret_vec_mod2);
}
