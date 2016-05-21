#include <iostream>
#include <string>
#include <immintrin.h> 
#include <getopt.h>
#include "Config.h"
#include "Timing.h"
using namespace std;

int seq_SW(const char* queryStr, int m, const char* subjectStr, int n);
int seq_scan_SW(const char* queryStr, int m, const char* subjectStr, int n);
int mod_striped_iterate_SW(const char* queryStr, int m, const char* subjectStr, int n);
int mod_striped_scan_SW(const char* queryStr, int m, const char* subjectStr, int n);
int mod_striped_merged_SW(const char* queryStr, int m, const char* subjectStr, int n);

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
    int ret_seq = seq_SW(queryStr.c_str(), queryStr.length(), subjectStr.c_str(), subjectStr.length());
    tstop  = dtime();
    ttime  = tstop - tstart;
    printf("time\tseq_orig\t%lf\t%d\n", ttime, ret_seq);

    tstart = dtime();
    int ret_sca = seq_scan_SW(queryStr.c_str(), queryStr.length(), subjectStr.c_str(), subjectStr.length());
    tstop  = dtime();
    ttime  = tstop - tstart;
    printf("time\tseq_scan\t%lf\t%d\n", ttime, ret_sca);

    tstart = dtime();
    int ret_vec_mod = mod_striped_scan_SW(queryStr.c_str(), queryStr.length(), subjectStr.c_str(), subjectStr.length());
    tstop  = dtime();
    ttime  = tstop - tstart;
    printf("time\tvec_scan\t%lf\t%d\n", ttime, ret_vec_mod);
    // printf("time\tvec_scan(avg)\t%lf\t%d\n", ttime/subjectStr.length(), ret_vec_mod);
    tstart = dtime();
    int ret_vec_mod1 = mod_striped_iterate_SW(queryStr.c_str(), queryStr.length(), subjectStr.c_str(), subjectStr.length());
    tstop  = dtime();
    ttime  = tstop - tstart;
    printf("time\tvec_iter\t%lf\t%d\n", ttime, ret_vec_mod1);
    // printf("time\tvec_iter(avg)\t%lf\t%d\n", ttime/subjectStr.length(), ret_vec_mod1);
    tstart = dtime();
    int ret_vec_mod2 = mod_striped_merged_SW(queryStr.c_str(), queryStr.length(), subjectStr.c_str(), subjectStr.length());
    tstop  = dtime();
    ttime  = tstop - tstart;
    printf("time\tvec_merg\t%lf\t%d\n", ttime, ret_vec_mod2);
}
