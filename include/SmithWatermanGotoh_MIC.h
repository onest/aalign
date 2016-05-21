#ifndef _SMITH_WATERMAN_GOTOH_MIC_H
#define _SMITH_WATERMAN_GOTOH_MIC_H

int sequential_SWAGotoh_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_SWAGotoh_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_SWAGotoh_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_SWAGotoh_MIC(const char* queryStr, int m, const char* subjectStr, int n);

#endif
