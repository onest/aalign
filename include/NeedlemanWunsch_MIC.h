#ifndef _NEEDLEMAN_WUNSCH_MIC_H
#define _NEEDLEMAN_WUNSCH_MIC_H

int sequential_NW_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_NW_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_NW_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_NW_MIC(const char* queryStr, int m, const char* subjectStr, int n);

#endif
