#ifndef _OPT_STR_ALIGN_DISTANCE_CPU_H
#define _OPT_STR_ALIGN_DISTANCE_CPU_H

int sequential_OSAD_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_OSAD_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_OSAD_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_OSAD_CPU(const char* queryStr, int m, const char* subjectStr, int n);

#endif
