#ifndef _SMITH_WATERMAN_GOTOH_CPU_H
#define _SMITH_WATERMAN_GOTOH_CPU_H

int sequential_SWAGotoh_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_SWAGotoh_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_SWAGotoh_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_SWAGotoh_CPU(const char* queryStr, int m, const char* subjectStr, int n);

#endif
