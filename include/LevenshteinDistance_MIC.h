#ifndef _LEVENSHTEIN_DISTANCE_MIC_H
#define _LEVENSHTEIN_DISTANCE_MIC_H

int sequential_LD_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_LD_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_LD_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_LD_MIC(const char* queryStr, int m, const char* subjectStr, int n);

#endif
