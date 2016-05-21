/*
* (c) 2015 Virginia Polytechnic Institute & State University (Virginia Tech)   
*                                                                              
*   This program is free software: you can redistribute it and/or modify       
*   it under the terms of the GNU General Public License as published by       
*   the Free Software Foundation, either version 3 of the License, or          
*   (at your option) any later version.                                        
*                                                                              
*   This program is distributed in the hope that it will be useful,            
*   but WITHOUT ANY WARRANTY; without even the implied warranty of             
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              
*   GNU General Public License for more details.                               
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

int seq_scan_SWAG(const char* queryStr, int m, const char* subjectStr, int n)
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
    int *SWAH1 = (int *)malloc(sizeof(int)*(m+1));
    int *SWAH2 = (int *)malloc(sizeof(int)*(m+1));
    int *SWAE1 = (int *)malloc(sizeof(int)*(m+1));
    // int *SWAE2 = (int *)malloc(sizeof(int)*(m+1));
    int *SWAF  = (int *)malloc(sizeof(int)*(m+1));
    fill_n(SWAH1, m+1, 0); // Init
    fill_n(SWAE1, m+1, 0); // Init

    int bufMax[3];
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m+1; j++)
        {
            if(j == 0)
            {
                // SWAE2[j] = 0; // Init
                SWAE1[j] = 0; // Init
                SWAH2[j] = 0; // Init
                // SWAF[j]  = 0; // Init
            }
            else
            {
                // SWAE2[j] = max(SWAE1[j]+GAP_EXT_PENALTY, SWAH1[j]+GAP_OPEN_PENALTY);
                SWAE1[j] = max(SWAE1[j]+GAP_EXT_PENALTY, SWAH1[j]+GAP_OPEN_PENALTY);
                // SWAF[j] = max(SWAF[j-1]+GAP_EXT_PENALTY, SWAH2[j-1]+GAP_OPEN_PENALTY);
                bufMax[0] = 0;
                // bufMax[1] = SWAE2[j];
                bufMax[1] = SWAE1[j];
                // bufMax[2] = SWAF[j];
#ifdef USE_BLOSUM
                // bufMax[2] = SWAH1[j-1] + BLOSUM62(blosum62char2int(queryStr[j-1]), 
                        // blosum62char2int(subjectStr[i]));;
                bufMax[2] = SWAH1[j-1] + profTDia[blosum62char2int(subjectStr[i])*m+j-1];
#else
                bufMax[2] = SWAH1[j-1] + (queryStr[j-1] == subjectStr[i]?MATCH_SCORE:MISMATCH_SCORE); 
#endif
                SWAH2[j] = *max_element(bufMax, bufMax+3);
                // score = max(score, SWAH2[j]);
            }
        }
        for(int j = 0; j < m+1; j++)
        {
            if(j == 0)
            {
                SWAF[j] = 0; // Init
            }
            else
            {
                SWAF[j] = max(SWAF[j-1] + GAP_EXT_PENALTY, SWAH2[j-1] + GAP_EXT_PENALTY);
            }
        }
        for(int j = 0; j < m+1; j++)
        {
            SWAH2[j] = max(SWAH2[j], SWAF[j]+GAP_OPEN_PENALTY-GAP_EXT_PENALTY);
            // if(j>0)
                // cout << SWAH2[j] << " ";
            score = max(score, SWAH2[j]);
        }
        // cout << endl;

        swap(SWAH1, SWAH2);
        // swap(SWAE1, SWAE2);
    }
    return score;
}

int seq_SWAG(const char* queryStr, int m, const char* subjectStr, int n)
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
    int *SWAH1 = (int *)malloc(sizeof(int)*(m+1));
    int *SWAH2 = (int *)malloc(sizeof(int)*(m+1));
    int *SWAE1 = (int *)malloc(sizeof(int)*(m+1));
    // int *SWAE2 = (int *)malloc(sizeof(int)*(m+1));
    // int *SWAF  = (int *)malloc(sizeof(int)*(m+1));
    // Kaixi new 
    int F_prev, F_curr;

    fill_n(SWAH1, m+1, 0); // Init
    fill_n(SWAE1, m+1, 0); // Init

    int bufMax[4];
    for(int i = 0; i < n; i++)
    {
        // int up = 0, left = 0, dia = 0, con = 0; 
#pragma vector always
        for(int j = 0; j < m+1; j++)
        {
            if(j == 0)
            {
                // SWAE2[j] = 0; // Init
                SWAE1[j] = 0; // Init
                SWAH2[j] = 0; // Init
                // SWAF[j]  = 0; // Init
                F_curr  = 0; // Init
            }
            else
            {
                // SWAE2[j] = max(SWAE1[j]+GAP_EXT_PENALTY, SWAH1[j]+GAP_OPEN_PENALTY);
                SWAE1[j] = max(SWAE1[j]+GAP_EXT_PENALTY, SWAH1[j]+GAP_OPEN_PENALTY);
                // SWAF[j] = max(SWAF[j-1]+GAP_EXT_PENALTY, SWAH2[j-1]+GAP_OPEN_PENALTY);
                F_curr = max(F_prev+GAP_EXT_PENALTY, SWAH2[j-1]+GAP_OPEN_PENALTY);
                bufMax[0] = 0;
                // bufMax[1] = SWAE2[j];
                bufMax[1] = SWAE1[j];
                // bufMax[2] = SWAF[j];
                bufMax[2] = F_curr;
#ifdef USE_BLOSUM
                // bufMax[3] = SWAH1[j-1] + BLOSUM62(blosum62char2int(queryStr[j-1]), 
                        // blosum62char2int(subjectStr[i]));;
                bufMax[3] = SWAH1[j-1] + profTDia[blosum62char2int(subjectStr[i])*m+j-1];
#else
                bufMax[3] = SWAH1[j-1] + (queryStr[j-1] == subjectStr[i]?MATCH_SCORE:MISMATCH_SCORE); 
#endif
                // if(max_element(bufMax, bufMax+4)==bufMax+0)
                    // con++;
                // else if(max_element(bufMax, bufMax+4)==bufMax+1)
                    // left++;
                // else if(max_element(bufMax, bufMax+4)==bufMax+2)
                    // up++;
                // else if(max_element(bufMax, bufMax+4)==bufMax+3)
                    // dia++;

                SWAH2[j] = *max_element(bufMax, bufMax+4);
#ifdef __DEBUG_SEQ
                cout << SWAH2[j] << " ";
#endif
                score = max(score, SWAH2[j]);
            }
            swap(F_prev, F_curr);
        }
        // cout << (float)up/m << endl;
#ifdef __DEBUG_SEQ
        cout << endl;
#endif
        swap(SWAH1, SWAH2);
        // swap(SWAE1, SWAE2);
    }
    
    return score;
}

int mod_striped_iterate_SWAG(const char* queryStr, int m, const char* subjectStr, int n)
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
    int *arrE  = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    fill_n(arrT1, m_striped, 0); 
    fill_n(arrE, m_striped, 0); 

    vec vTDia, vTLeft, vTUp, vT;
    vec vE, vF;
    vec vMax = broadcast(INT_MIN);
    vec vGapTLeft = broadcast(GAP_OPEN_PENALTY); 
    vec vGapTUp = broadcast(GAP_OPEN_PENALTY); 
    vec vConst = broadcast(0);
    vec vGapE = broadcast(GAP_EXT_PENALTY);
    vec vGapF = broadcast(GAP_EXT_PENALTY);
    for(int i = 0; i < n; i++)
    {
        int si = blosum62char2int(subjectStr[i]);
        vTDia = right_shift_x_fill(arrT1+(k-1)*vec_length, 1, 0);
        vTUp = set_vector(m_striped, 0, GAP_OPEN_PENALTY);
        vTUp = add_vector(vTUp, vGapTUp);
        vF = set_vector(m_striped, 0, GAP_EXT_PENALTY);
        vF = add_vector(vF, vGapF);
        vF = max_vector(vF, vTUp);
        for(int j = 0; j < k; j++)
        {
            vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
            vTLeft = add_array(arrT1+j*vec_length, vGapTLeft);
            vE = add_array(arrE+j*vec_length, vGapE);
            vE = max_vector(vE, vTLeft);
            store_vector(arrE+j*vec_length, vE);
            vT = max_vector(vTDia, vE, vF, vConst);
            store_vector(arrT2+j*vec_length, vT);
            vMax = max_vector(vMax, vT);
            vTDia = load_vector(arrT1+j*vec_length);
            vTUp = vT;
            vTUp = add_vector(vTUp, vGapTUp);
            vF = add_vector(vF, vGapF);
            vF = max_vector(vF, vTUp);
        }
        vF = right_shift_x_fill(vF, 1, GAP_OPEN_PENALTY);
        int j = 0;
        vT = load_vector(arrT2+j*vec_length); 
        int count = 0;
        // float ratio = 0.f;
        // while(count < vec_length-1)
        // while(count < 1)
        while(influence_test_max(vF, add_vector(vT, sub_vector(vGapTUp, vGapF))))
        {
            vT = max_vector(vF, vT);
            store_vector(arrT2+j*vec_length, vT);
            vMax = max_vector(vMax, vT);
            vF = add_vector(vF, vGapF);
            if(++j >= k)
            {
                vF = right_shift_x_fill(vF, 1, GAP_OPEN_PENALTY);
                j = 0;
                // count = 0;
                // ratio += 1.f;
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
    _mm_free(arrE);
    _mm_free(profTDia);
    free(stripedQueryStr);

    return score;
}

int mod_striped_scan_SWAG(const char* queryStr, int m, const char* subjectStr, int n)
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
    int *arrE  = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    int *arrScan = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    fill_n(arrT1, m_striped, 0); 
    fill_n(arrE, m_striped, 0); 

    vec vTDia, vTLeft, vTUp, vT;
    vec vE;
    vec vMax = broadcast(INT_MIN);
    vec vGapTLeft = broadcast(GAP_OPEN_PENALTY); 
    vec vGapE = broadcast(GAP_EXT_PENALTY);
    vec vConst = broadcast(0);
    for(int i = 0; i < n; i++)
    {
        int si = blosum62char2int(subjectStr[i]);
        vTDia = right_shift_x_fill(arrT1+(k-1)*vec_length, 1, 0);
        for(int j = 0; j < k; j++)
        {
            vTLeft = add_array(arrT1+j*vec_length, vGapTLeft);
            vE = add_array(arrE+j*vec_length, vGapE);
            vE = max_vector(vE, vTLeft);
            store_vector(arrE+j*vec_length, vE);
            vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
            vT = max_vector(vTDia, vE, vConst);
            // vMax = max_vector(vMax, vT); // test
            store_vector(arrT2+j*vec_length, vT);
            vTDia = load_vector(arrT1+j*vec_length);
        }
        max_scan_striped_array(arrT2, arrScan, m_striped, 0, GAP_EXT_PENALTY, GAP_OPEN_PENALTY);
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
    _mm_free(arrE);
    _mm_free(arrScan);
    _mm_free(profTDia);
    free(stripedQueryStr);

    return score;
}

int mod_striped_merged_SWAG(const char* queryStr, int m, const char* subjectStr, int n)
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
    int *arrE  = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    int *arrScan = (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    fill_n(arrT1, m_striped, 0); 
    fill_n(arrE, m_striped, 0); 

    vec vTDia, vTLeft, vTUp, vT;
    vec vE, vF;
    vec vMax = broadcast(INT_MIN);
    vec vGapTLeft = broadcast(GAP_OPEN_PENALTY); 
    vec vGapTUp = broadcast(GAP_OPEN_PENALTY); 
    vec vConst = broadcast(0);
    vec vGapE = broadcast(GAP_EXT_PENALTY);
    vec vGapF = broadcast(GAP_EXT_PENALTY);
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
            vTUp = set_vector(m_striped, 0, GAP_OPEN_PENALTY);
            vTUp = add_vector(vTUp, vGapTUp);
            vF = set_vector(m_striped, 0, GAP_EXT_PENALTY);
            vF = add_vector(vF, vGapF);
            vF = max_vector(vF, vTUp);
            for(int j = 0; j < k; j++)
            {
                vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
                vTLeft = add_array(arrT1+j*vec_length, vGapTLeft);
                vE = add_array(arrE+j*vec_length, vGapE);
                vE = max_vector(vE, vTLeft);
                store_vector(arrE+j*vec_length, vE);
                vT = max_vector(vTDia, vE, vF, vConst);
                store_vector(arrT2+j*vec_length, vT);
                vMax = max_vector(vMax, vT);
                vTDia = load_vector(arrT1+j*vec_length);
                vTUp = vT;
                vTUp = add_vector(vTUp, vGapTUp);
                vF = add_vector(vF, vGapF);
                vF = max_vector(vF, vTUp);
            }
            vF = right_shift_x_fill(vF, 1, GAP_OPEN_PENALTY);
            int j = 0;
            vT = load_vector(arrT2+j*vec_length); 
            int count = 0;
            // while(count < vec_length-1)
            // while(count < 1)
            while(influence_test_max(vF, add_vector(vT, sub_vector(vGapTUp, vGapF))))
            {
                vT = max_vector(vF, vT);
                store_vector(arrT2+j*vec_length, vT);
                vMax = max_vector(vMax, vT);
                vF = add_vector(vF, vGapF);
                if(++j >= k)
                {
                    vF = right_shift_x_fill(vF, 1, GAP_OPEN_PENALTY);
                    j = 0;
                    count = 0;
                    ratio += 1.f;
                }
                vT = load_vector(arrT2+j*vec_length);
                count++;
            }
            ratio += (float)count / k;

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
                vE = add_array(arrE+j*vec_length, vGapE);
                vE = max_vector(vE, vTLeft);
                store_vector(arrE+j*vec_length, vE);
                vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
                vT = max_vector(vTDia, vE, vConst);
                // vMax = max_vector(vMax, vT); // test
                store_vector(arrT2+j*vec_length, vT);
                vTDia = load_vector(arrT1+j*vec_length);
            }
            max_scan_striped_array(arrT2, arrScan, m_striped, 0, GAP_EXT_PENALTY, GAP_OPEN_PENALTY);
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
            if(ratio <= 1.0f)
                flag_iter = true;
            else
                acc = 0;
        }
        if(flag_iter && ratio > 1.0f)
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
    _mm_free(arrE);
    _mm_free(arrScan);
    _mm_free(profTDia);
    free(stripedQueryStr);

    return score;
}
