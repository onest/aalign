#ifndef _SMITH_WATERMAN_MIC_H
#define _SMITH_WATERMAN_MIC_H

int sequential_SWA_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_SWA_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_SWA_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_SWA_MIC(const char* queryStr, int m, const char* subjectStr, int n);

#endif
