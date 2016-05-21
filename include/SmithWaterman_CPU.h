#ifndef _SMITH_WATERMAN_CPU_H
#define _SMITH_WATERMAN_CPU_H

int sequential_SWA_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_SWA_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_SWA_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_SWA_CPU(const char* queryStr, int m, const char* subjectStr, int n);

#endif
