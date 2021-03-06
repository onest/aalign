/*
* (c) 2015 Virginia Polytechnic Institute & State University (Virginia Tech)   
*                                                                              
*   This program is free software: you can redistribute it and/or modify       
*   it under the terms of the GNU General Public License as published by       
*   the Free Software Foundation, version 2.1                                  
*                                                                              
*   This program is distributed in the hope that it will be useful,            
*   but WITHOUT ANY WARRANTY; without even the implied warranty of             
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              
*   GNU General Public License, version 2.1, for more details.                 
*                                                                              
*   You should have received a copy of the GNU General Public License          
*                                                                              
*  Acknowledgement                                                             
*   This research was supported in part by NSF-BIGDATA program via             
*   IIS-1247693                                                                
*                                                                              
*/


#include <iostream>
#include <climits>
#include "Config.h"
#include "Modules.h"
using namespace std;
// #define __DEBUG_SEQ
// #define __DEBUG_ITER
// #define __DEBUG_SCAN
// #define __DEBUG_MERG
int seq_scan_SW(const char* queryStr, int m, const char* subjectStr, int n)
{
    int8_t* profTDia = (int8_t*)malloc(sizeof(int8_t)*32*m);
    fill_n(profTDia, 32*m, 0);
    for(int i = 0; i < SUBS_MATRIX_SIZE; i++)
    {
        for(int j = 0; j < m; j++)
        {
            profTDia[i*m+j] = (int8_t)(BLOSUM62(i, blosum62char2int(queryStr[j])));
        }
    }

    int score = INT_MIN;
    int *SWBuf1 = (int *)malloc(sizeof(int)*(m+1));
    int *SWBuf2 = (int *)malloc(sizeof(int)*(m+1));
    int *SWBuf2_scan = (int *)malloc(sizeof(int)*(m+1));
    fill_n(SWBuf1, m+1, 0); // Init 

    int bufMax[3];
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m+1; j++)
        {
            if(j == 0)
            {
                SWBuf2[j] = 0; // Init
            }
            else
            {
                bufMax[0] = SWBuf1[j  ] + GAP_PENALTY;
                bufMax[1] = 0;
#ifdef USE_BLOSUM
                // bufMax[1] = NWBuf1[j-1] + BLOSUM62(blosum62char2int(queryStr[j-1]), 
                        // blosum62char2int(subjectStr[i]));
                bufMax[2] = SWBuf1[j-1] + profTDia[blosum62char2int(subjectStr[i])*m+j-1];
#else
                bufMax[2] = SWBuf1[j-1] + (queryStr[j-1] == subjectStr[i]?MATCH_SCORE:MISMATCH_SCORE); 
#endif
                SWBuf2[j] = *max_element(bufMax, bufMax+3);
                score = max(score, SWBuf2[j]);
            }
        }
        for(int j = 0; j < m+1; j++)
        {
            if(j == 0)
            {
                SWBuf2_scan[j] = 0; // Init
            }
            else
            {
                SWBuf2_scan[j] = max(SWBuf2_scan[j-1] + GAP_PENALTY, SWBuf2[j-1] + GAP_PENALTY);
            }
        }
        for(int j = 0; j < m+1; j++)
        {
            SWBuf2[j] = max(SWBuf2[j], SWBuf2_scan[j]);
            // if(j>0)
                // cout << NWBuf2[j] << " ";
            score = max(score, SWBuf2[j]);
        }
        // cout << endl;
        swap(SWBuf1, SWBuf2);
    }

    free(SWBuf1);
    free(SWBuf2);
    free(SWBuf2_scan);
    return score;
}
int seq_SW(const char* queryStr, int m, const char* subjectStr, int n)
{
    // Kaixi new profile
    int8_t* profTDia = (int8_t*)malloc(sizeof(int8_t)*32*m);
    fill_n(profTDia, 32*m, 0);
    for(int i = 0; i < SUBS_MATRIX_SIZE; i++)
    {
        for(int j = 0; j < m; j++)
        {
            profTDia[i*m+j] = (int8_t)(BLOSUM62(i, blosum62char2int(queryStr[j])));
        }
    }    
    int score = INT_MIN;
    int *SWABuf1 = (int *)malloc(sizeof(int)*(m+1));
    int *SWABuf2 = (int *)malloc(sizeof(int)*(m+1));
    fill_n(SWABuf1, m+1, 0); // Init 

    int bufMax[4];
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m+1; j++)
        {
            if(j == 0)
            {
                SWABuf2[j] = 0; // Init
            }
            else
            {
                bufMax[0] = 0;
                bufMax[1] = SWABuf2[j-1] + GAP_PENALTY;
                bufMax[2] = SWABuf1[j  ] + GAP_PENALTY;
#ifdef USE_BLOSUM
                // bufMax[3] = SWABuf1[j-1] + BLOSUM62(blosum62char2int(queryStr[j-1]), 
                        // blosum62char2int(subjectStr[i]));
                bufMax[3] = SWABuf1[j-1] + profTDia[blosum62char2int(subjectStr[i])*m+j-1];
#else
                bufMax[3] = SWABuf1[j-1] + (queryStr[j-1] == subjectStr[i]?MATCH_SCORE:MISMATCH_SCORE); 
#endif
                SWABuf2[j] = *max_element(bufMax, bufMax+4);
#ifdef __DEBUG_SEQ
                cout << SWABuf2[j] << " ";
#endif
                score = max(score, SWABuf2[j]);
            }
        }
#ifdef __DEBUG_SEQ
        cout << endl;
#endif
        swap(SWABuf1, SWABuf2);
    }
    free(SWABuf1);
    free(SWABuf2);
    return score;
}

int mod_striped_iterate_SW(const char* queryStr, int m, const char* subjectStr, int n)
{
    // Step 1: Stripe the query sequence
    int w = vec_length;
    int k = (m + w - 1) / w;
    int m_striped = k * w;
    char *stripedQueryStr = (char *)malloc(sizeof(char) * (m_striped + 1));
    for(int i = 0; i < w; i++)
    {
        for(int j = 0; j < k; j++)
        {
            stripedQueryStr[j * w + i] = (i*k+j>=m)?'-':queryStr[i*k+j];
        }
    }
    stripedQueryStr[m_striped]='\0';

    // assume the maximum number of different characters is 32
    int8_t* profTDia = (int8_t*)_mm_malloc(sizeof(int8_t)*32*m_striped, vec_length);
    fill_n(profTDia, 32*m_striped, 0);
    for(int i = 0; i < SUBS_MATRIX_SIZE; i++)
    {
        for(int j = 0; j < m_striped; j++)
        {
            profTDia[i*m_striped+j] = (int8_t)(BLOSUM62(i, blosum62char2int(stripedQueryStr[j])));
        }
    }

    int score = 1;
    int *arrT1 = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    int *arrT2 = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length); 
    for(int j = 0; j < k; j++) // Init
    {
        for(int i = 0; i < w; i++)
        {
            arrT1[j * w + i] = 0;
        }
    }
    // fill_n(arrE, m_striped, 0); 

    vec vTDia, vTLeft, vTUp, vT;
    // vec vE, vF;
    vec vMax = broadcast(INT_MIN);
    vec vGapTLeft = broadcast(GAP_PENALTY); 
    vec vGapTUp = broadcast(GAP_PENALTY); 
    vec vConst = broadcast(0);
    // vec vGapE = broadcast(GAP_EXT_PENALTY);
    // vec vGapF = broadcast(GAP_EXT_PENALTY);
    for(int i = 0; i < n; i++)
    {
        int si = blosum62char2int(subjectStr[i]);
        vTDia = right_shift_x_fill(arrT1+(k-1)*vec_length, 1, 0);
        vTUp = set_vector(m_striped, 0, GAP_PENALTY);
        vTUp = add_vector(vTUp, vGapTUp);
        // vF = set_vector(m_striped, 0, GAP_EXT_PENALTY);
        // vF = add_vector(vF, vGapF);
        // vF = max_vector(vF, vTUp);
        for(int j = 0; j < k; j++)
        {
            vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
            vTLeft = add_array(arrT1+j*vec_length, vGapTLeft);
            // vE = add_array(arrE+j*vec_length, vGapE);
            // vE = max_vector(vE, vTLeft);
            // store_vector(arrE+j*vec_length, vE);
            vT = max_vector(vTDia, vTLeft, vTUp, vConst);
            store_vector(arrT2+j*vec_length, vT);
            vMax = max_vector(vMax, vT);
            vTDia = load_vector(arrT1+j*vec_length);
            vTUp = vT;
            vTUp = add_vector(vTUp, vGapTUp);
            // vF = add_vector(vF, vGapF);
            // vF = max_vector(vF, vTUp);
        }
        vTUp = right_shift_x_fill(vTUp, 1, GAP_PENALTY);
        int j = 0;
        vT = load_vector(arrT2+j*vec_length); 
        // float ratio = 0.f;
        // int count = 0;
        // while(count < vec_length-1)
        // while(count < 15)
        while(influence_test_max(vTUp, add_vector(vT, sub_vector(vGapTUp, vGapTUp))))
        {
            vT = max_vector(vTUp, vT);
            store_vector(arrT2+j*vec_length, vT);
            vMax = max_vector(vMax, vT);
            vTUp = add_vector(vTUp, vGapTUp);
            if(++j >= k)
            {
                vTUp = right_shift_x_fill(vTUp, 1, GAP_PENALTY);
                j = 0;
                // ratio += 1.f;
                // count = 0;
            }
            vT = load_vector(arrT2+j*vec_length);
            // count++;
        }
        // ratio += (float)count / k;
        // cout << ratio << endl;

#ifdef __DEBUG_ITER
        for(int ii = 0; ii < w; ii++)
        {
            for(int jj = 0; jj < k; jj++)
            {
                if(ii*k+jj < m)
                    cout << arrT2[jj * w + ii] << " "; 
            }
        }
        cout << endl;
#endif
        std::swap(arrT1, arrT2);
    }
    score = reduce_max(vMax);
    _mm_free(arrT1);
    _mm_free(arrT2);
    // _mm_free(arrE);
    _mm_free(profTDia);
    free(stripedQueryStr);

    return score;
}

int mod_striped_scan_SW(const char* queryStr, int m, const char* subjectStr, int n)
{
    // Step 1: Stripe the query sequence
    int w = vec_length;
    int k = (m + w - 1) / w;
    int m_striped = k * w;
    char *stripedQueryStr = (char *)malloc(sizeof(char) * (m_striped + 1));
    for(int i = 0; i < w; i++)
    {
        for(int j = 0; j < k; j++)
        {
            stripedQueryStr[j * w + i] = (i*k+j>=m)?'-':queryStr[i*k+j];
        }
    }
    stripedQueryStr[m_striped]='\0';

    // assume the maximum number of different characters is 32
    int8_t* profTDia = (int8_t*)_mm_malloc(sizeof(int8_t)*32*m_striped, vec_length);
    fill_n(profTDia, 32*m_striped, 0);
    for(int i = 0; i < SUBS_MATRIX_SIZE; i++)
    {
        for(int j = 0; j < m_striped; j++)
        {
            profTDia[i*m_striped+j] = (int8_t)(BLOSUM62(i, blosum62char2int(stripedQueryStr[j])));
        }
    }

    int score = 1;
    int *arrT1 = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    int *arrT2 = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length); 
    int *arrScan = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    for(int j = 0; j < k; j++) // Init
    {
        for(int i = 0; i < w; i++)
        {
            arrT1[j * w + i] = 0;
        }
    }

    vec vTDia, vTLeft, vTUp, vT;
    // vec vE;
    vec vMax = broadcast(INT_MIN);
    vec vGapTLeft = broadcast(GAP_PENALTY); 
    // vec vGapE = broadcast(GAP_EXT_PENALTY);
    vec vConst = broadcast(0);
    for(int i = 0; i < n; i++)
    {
        int si = blosum62char2int(subjectStr[i]);
        vTDia = right_shift_x_fill(arrT1+(k-1)*vec_length, 1, 0);
        for(int j = 0; j < k; j++)
        {
            vTLeft = add_array(arrT1+j*vec_length, vGapTLeft);
            // vE = add_array(arrE+j*vec_length, vGapE);
            // vE = max_vector(vE, vTLeft);
            // store_vector(arrE+j*vec_length, vE);
            vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
            vT = max_vector(vTDia, vTLeft, vConst);
            store_vector(arrT2+j*vec_length, vT);
            vTDia = load_vector(arrT1+j*vec_length);
        }
        max_scan_striped_array(arrT2, arrScan, m_striped, 0, GAP_PENALTY, GAP_PENALTY);
        for(int j = 0; j < k; j++)
        {
            vTUp = load_vector(arrScan+j*vec_length);
            vT = load_vector(arrT2+j*vec_length);
            vT = max_vector(vT, vTUp);
            vMax = max_vector(vMax, vT);
            store_vector(arrT2+j*vec_length, vT);
        }
#ifdef __DEBUG_SCAN
        for(int ii = 0; ii < w; ii++)
        {
            for(int jj = 0; jj < k; jj++)
            {
                if(ii*k+jj < m)
                    cout << arrT2[jj * w + ii] << " "; 
            }
        }
        cout << endl;
#endif
        std::swap(arrT1, arrT2);
    } 

    score = reduce_max(vMax);
    _mm_free(arrT1);
    _mm_free(arrT2);
    _mm_free(arrScan);
    _mm_free(profTDia);
    free(stripedQueryStr);

    return score;
}

int mod_striped_merged_SW(const char* queryStr, int m, const char* subjectStr, int n)
{
    int probe_stride = m / 20;
    // Step 1: Stripe the query sequence
    int w = vec_length;
    int k = (m + w - 1) / w;
    int m_striped = k * w;
    char *stripedQueryStr = (char *)malloc(sizeof(char) * (m_striped + 1));
    for(int i = 0; i < w; i++)
    {
        for(int j = 0; j < k; j++)
        {
            stripedQueryStr[j * w + i] = (i*k+j>=m)?'-':queryStr[i*k+j];
        }
    }
    stripedQueryStr[m_striped]='\0';

    // assume the maximum number of different characters is 32
    int8_t* profTDia = (int8_t*)_mm_malloc(sizeof(int8_t)*32*m_striped, vec_length);
    fill_n(profTDia, 32*m_striped, 0);
    for(int i = 0; i < SUBS_MATRIX_SIZE; i++)
    {
        for(int j = 0; j < m_striped; j++)
        {
            profTDia[i*m_striped+j] = (int8_t)(BLOSUM62(i, blosum62char2int(stripedQueryStr[j])));
        }
    }

    int score = 1;
    int *arrT1 = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    int *arrT2 = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length); 
    int *arrScan = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    for(int j = 0; j < k; j++) // Init
    {
        for(int i = 0; i < w; i++)
        {
            arrT1[j * w + i] = 0;
        }
    }
    // fill_n(arrE, m_striped, 0); 

    vec vTDia, vTLeft, vTUp, vT;
    // vec vE, vF;
    vec vMax = broadcast(INT_MIN);
    vec vGapTLeft = broadcast(GAP_PENALTY); 
    vec vGapTUp = broadcast(GAP_PENALTY); 
    vec vConst = broadcast(0);
    // vec vGapE = broadcast(GAP_EXT_PENALTY);
    // vec vGapF = broadcast(GAP_EXT_PENALTY);
    bool flag_iter = true;
    bool probe_flag_iter = false;
    float ratio = 0.f;
    int acc = 0;
    for(int i = 0; i < n; i++)
    {
        if(flag_iter || probe_flag_iter)
        {
            // cout << "0\n";
            ratio = 0.f;
            int si = blosum62char2int(subjectStr[i]);
            vTDia = right_shift_x_fill(arrT1+(k-1)*vec_length, 1, 0);
            vTUp = set_vector(m_striped, 0, GAP_PENALTY);
            vTUp = add_vector(vTUp, vGapTUp);
            // vF = set_vector(m_striped, 0, GAP_EXT_PENALTY);
            // vF = add_vector(vF, vGapF);
            // vF = max_vector(vF, vTUp);
            for(int j = 0; j < k; j++)
            {
                vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
                vTLeft = add_array(arrT1+j*vec_length, vGapTLeft);
                // vE = add_array(arrE+j*vec_length, vGapE);
                // vE = max_vector(vE, vTLeft);
                // store_vector(arrE+j*vec_length, vE);
                vT = max_vector(vTDia, vTLeft, vTUp, vConst);
                store_vector(arrT2+j*vec_length, vT);
                vMax = max_vector(vMax, vT);
                vTDia = load_vector(arrT1+j*vec_length);
                vTUp = vT;
                vTUp = add_vector(vTUp, vGapTUp);
                // vF = add_vector(vF, vGapF);
                // vF = max_vector(vF, vTUp);
            }
            vTUp = right_shift_x_fill(vTUp, 1, GAP_PENALTY);
            int j = 0;
            vT = load_vector(arrT2+j*vec_length); 
            int count = 0;
            // while(count < vec_length-1)
            // while(count < 15)
            while(influence_test_max(vTUp, add_vector(vT, sub_vector(vGapTUp, vGapTUp))))
            {
                vT = max_vector(vTUp, vT);
                store_vector(arrT2+j*vec_length, vT);
                vMax = max_vector(vMax, vT);
                vTUp = add_vector(vTUp, vGapTUp);
                if(++j >= k)
                {
                    vTUp = right_shift_x_fill(vTUp, 1, GAP_PENALTY);
                    j = 0;
                    count = 0;
                    ratio += 1.f;
                }
                vT = load_vector(arrT2+j*vec_length);
                count++;
            }
            ratio += (float)count / k;
            // cout << ratio << endl;
#ifdef __DEBUG_MERG
            for(int ii = 0; ii < w; ii++)
            {
                for(int jj = 0; jj < k; jj++)
                {
                    if(ii*k+jj < m)
                        cout << arrT2[jj * w + ii] << " "; 
                }
            }
            cout << endl;
#endif
            std::swap(arrT1, arrT2);
        }else
        {
            // cout << "1\n";
            acc ++;
            int si = blosum62char2int(subjectStr[i]);
            vTDia = right_shift_x_fill(arrT1+(k-1)*vec_length, 1, 0);
            for(int j = 0; j < k; j++)
            {
                vTLeft = add_array(arrT1+j*vec_length, vGapTLeft);
                // vE = add_array(arrE+j*vec_length, vGapE);
                // vE = max_vector(vE, vTLeft);
                // store_vector(arrE+j*vec_length, vE);
                vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
                vT = max_vector(vTDia, vTLeft, vConst);
                store_vector(arrT2+j*vec_length, vT);
                vTDia = load_vector(arrT1+j*vec_length);
            }
            max_scan_striped_array(arrT2, arrScan, m_striped, 0, GAP_PENALTY, GAP_PENALTY);
            for(int j = 0; j < k; j++)
            {
                vTUp = load_vector(arrScan+j*vec_length);
                vT = load_vector(arrT2+j*vec_length);
                vT = max_vector(vT, vTUp);
                vMax = max_vector(vMax, vT);
                store_vector(arrT2+j*vec_length, vT);
            }
#ifdef __DEBUG_MERG
            for(int ii = 0; ii < w; ii++)
            {
                for(int jj = 0; jj < k; jj++)
                {
                    if(ii*k+jj < m)
                        cout << arrT2[jj * w + ii] << " "; 
                }
            }
            cout << endl;
#endif
            std::swap(arrT1, arrT2);
        
        }

        if(probe_flag_iter)
        {
            probe_flag_iter = false;
            if(ratio <= 3.0f)
                flag_iter = true;
            else
                acc = 0;
        }
        if(flag_iter && ratio > 3.0f)
        {
            flag_iter = false;
            acc = 0;
        } else if(!flag_iter && acc == probe_stride)
        {
            probe_flag_iter = true;
        }
    }
    score = reduce_max(vMax);
    _mm_free(arrT1);
    _mm_free(arrT2);
    // _mm_free(arrE);
    _mm_free(arrScan);
    _mm_free(profTDia);
    free(stripedQueryStr);

    return score;
}


