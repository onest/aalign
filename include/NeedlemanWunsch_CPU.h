#ifndef _NEEDLEMAN_WUNSCH_CPU_H
#define _NEEDLEMAN_WUNSCH_CPU_H

int sequential_NW_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_NW_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_NW_CPU(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_NW_CPU(const char* queryStr, int m, const char* subjectStr, int n);

#endif
