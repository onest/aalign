#ifndef _OPT_STR_ALIGN_DISTANCE_MIC_H
#define _OPT_STR_ALIGN_DISTANCE_MIC_H

int sequential_OSAD_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int sequential_scan_OSAD_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_OSAD_MIC(const char* queryStr, int m, const char* subjectStr, int n);
int vectorized_scan_OSAD_MIC(const char* queryStr, int m, const char* subjectStr, int n);

#endif
