#ifndef _LEVENSHTEIN_DISTANCE_CPU_H
#define _LEVENSHTEIN_DISTANCE_CPU_H

int sequential_LD_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_LD_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_LD_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_LD_CPU(const char* queryStr, int m, const char* subjectStr, int n);

#endif
