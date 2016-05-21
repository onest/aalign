#include <iostream>
#include <climits>
#include "Config.h"
#include "Modules.h"
using namespace std;
// #define __DEBUG_SEQ
// #define __DEBUG_ITER
// #define __DEBUG_SCAN
// #define __DEBUG_MERG

int seq_scan_NWAG(const char* queryStr, int m, const char* subjectStr, int n)
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
    int *NWAH1 = (int *)malloc(sizeof(int)*(m+1));
    int *NWAH2 = (int *)malloc(sizeof(int)*(m+1));
    int *NWAE1 = (int *)malloc(sizeof(int)*(m+1));
    // int *SWAE2 = (int *)malloc(sizeof(int)*(m+1));
    int *NWAF  = (int *)malloc(sizeof(int)*(m+1));
    // fill_n(NWAH1, m+1, 0); // Init
    // fill_n(NWAE1, m+1, 0); // Init
    for(int i = 0; i < m+1; i++) 
        NWAH1[i] = INT_MIN/2;
    for(int i = 0; i < m+1; i++) 
        NWAE1[i] = GAP_OPEN_PENALTY + GAP_EXT_PENALTY * (i-1);
    NWAH1[0] = 0;
    NWAE1[0] = INT_MIN/2;

    int bufMax[2];
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m+1; j++)
        {
            if(j == 0)
            {
                // SWAE2[j] = 0; // Init
                NWAE1[j] = INT_MIN/2; // Init
                NWAH2[j] = INT_MIN/2; // Init
                // SWAF[j]  = 0; // Init
            }
            else
            {
                // SWAE2[j] = max(SWAE1[j]+GAP_EXT_PENALTY, SWAH1[j]+GAP_OPEN_PENALTY);
                NWAE1[j] = max(NWAE1[j]+GAP_EXT_PENALTY, NWAH1[j]+GAP_OPEN_PENALTY);
                // SWAF[j] = max(SWAF[j-1]+GAP_EXT_PENALTY, SWAH2[j-1]+GAP_OPEN_PENALTY);
                // bufMax[0] = 0;
                // bufMax[1] = SWAE2[j];
                bufMax[0] = NWAE1[j];
                // bufMax[2] = SWAF[j];
#ifdef USE_BLOSUM
                // bufMax[2] = SWAH1[j-1] + BLOSUM62(blosum62char2int(queryStr[j-1]), 
                        // blosum62char2int(subjectStr[i]));;
                bufMax[1] = NWAH1[j-1] + profTDia[blosum62char2int(subjectStr[i])*m+j-1];
#else
                bufMax[1] = NWAH1[j-1] + (queryStr[j-1] == subjectStr[i]?MATCH_SCORE:MISMATCH_SCORE); 
#endif
                NWAH2[j] = *max_element(bufMax, bufMax+2);
                // score = max(score, SWAH2[j]);
            }
        }
        for(int j = 0; j < m+1; j++)
        {
            if(j == 0)
            {
                NWAF[j] = 0; // Init
            }
            else
            {
                NWAF[j] = max(NWAF[j-1] + GAP_EXT_PENALTY, NWAH2[j-1] + GAP_EXT_PENALTY);
            }
        }
        for(int j = 0; j < m+1; j++)
        {
            NWAH2[j] = max(NWAH2[j], NWAF[j]+GAP_OPEN_PENALTY-GAP_EXT_PENALTY);
            // if(j>0)
                // cout << SWAH2[j] << " ";
            // score = max(score, NWAH2[j]);
        }
        // cout << endl;

        swap(NWAH1, NWAH2);
        // swap(SWAE1, SWAE2);
    }
    score = NWAH1[m];
    return score;

}

int seq_NWAG(const char* queryStr, int m, const char* subjectStr, int n)
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
    int *NWBufH1 = (int *)malloc(sizeof(int)*(m+1));
    int *NWBufH2 = (int *)malloc(sizeof(int)*(m+1));
    int *NWBufE1 = (int *)malloc(sizeof(int)*(m+1));
    // int *NWBufE2 = (int *)malloc(sizeof(int)*(m+1));
    // int *NWBufF  = (int *)malloc(sizeof(int)*(m+1));

    for(int i = 0; i < m+1; i++) 
        NWBufH1[i] = INT_MIN/2;
    for(int i = 0; i < m+1; i++) 
        NWBufE1[i] = GAP_OPEN_PENALTY + GAP_EXT_PENALTY * (i-1);
    NWBufH1[0] = 0;
    NWBufE1[0] = INT_MIN/2;
    int F_prev, F_curr;

    int bufMax[3];
    for(int i = 0; i < n; i++)
    {
        // int up = 0, left = 0, dia = 0; 
        for(int j = 0; j < m+1; j++)
        {
            if(j == 0)
            {
                NWBufH2[j] = INT_MIN/2; 
                // NWBufE2[i] = INT_MIN;
                // NWBufF[j]  = GAP_OPEN_PENALTY + GAP_EXT_PENALTY * i;
                F_curr  = GAP_OPEN_PENALTY + GAP_EXT_PENALTY * i; // Init
            }
            else
            {
                // NWBufE2[j] = max(NWBufE1[j]+GAP_EXT_PENALTY, NWBufH1[j]+GAP_OPEN_PENALTY);
                NWBufE1[j] = max(NWBufE1[j]+GAP_EXT_PENALTY, NWBufH1[j]+GAP_OPEN_PENALTY);
                // NWBufF[j]  = max(NWBufF[j-1]+GAP_EXT_PENALTY, NWBufH2[j-1]+GAP_OPEN_PENALTY);
                F_curr = max(F_prev+GAP_EXT_PENALTY, NWBufH2[j-1]+GAP_OPEN_PENALTY);
                // bufMax[0] = NWBufE2[j];
                bufMax[0] = NWBufE1[j];
                // bufMax[1] = NWBufF[j];
                bufMax[1] = F_curr;
#ifdef USE_BLOSUM
                // bufMax[2] = NWBufH1[j-1] + BLOSUM62(blosum62char2int(queryStr[j-1]), 
                        // blosum62char2int(subjectStr[i]));
                bufMax[2] = NWBufH1[j-1] + profTDia[blosum62char2int(subjectStr[i])*m+j-1];
#else
                bufMax[2] = NWBufH1[j-1] + (queryStr[j-1] == subjectStr[i]?MATCH_SCORE:MISMATCH_SCORE); 
#endif
                // if(max_element(bufMax, bufMax+3)==bufMax+0)
                    // up++;
                // else if(max_element(bufMax, bufMax+3)==bufMax+1)
                    // left++;
                // else if(max_element(bufMax, bufMax+3)==bufMax+2)
                    // dia++;

                NWBufH2[j] = *max_element(bufMax, bufMax+3);
                // score = max(score, NWBufH2[j]);
#ifdef __DEBUG_SEQ
                cout << NWBufH2[j] << " ";
#endif
            }
            swap(F_prev, F_curr);
        }
#ifdef __DEBUG_SEQ
        cout << endl;
#endif
        // cout << (float)left/m << endl;
        swap(NWBufH1, NWBufH2);
        // swap(NWBufE1, NWBufE2);
    }
    score = NWBufH1[m];
    free(NWBufH1);
    free(NWBufH2);
    free(NWBufE1);
    // free(NWBufE2);
    // free(NWBufF);
    return score;
}

int mod_striped_iterate_NWAG(const char* queryStr, int m, const char* subjectStr, int n)
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
    for(int j = 0; j < k; j++) // Init
    {
        for(int i = 0; i < w; i++)
        {
            arrE[j * w + i] = GAP_OPEN_PENALTY + GAP_EXT_PENALTY * (i * k + j + 1 - 1);
        }
    }
    fill_n(arrT1, m_striped, INT_MIN/2); 

    vec vTDia, vTLeft, vTUp, vT;
    vec vE, vF;
    // vec vMax = broadcast(INT_MIN);
    vec vGapTLeft = broadcast(GAP_OPEN_PENALTY); 
    vec vGapTUp = broadcast(GAP_OPEN_PENALTY); 
    // vec vConst = broadcast(0);
    vec vGapE = broadcast(GAP_EXT_PENALTY);
    vec vGapF = broadcast(GAP_EXT_PENALTY);
    for(int i = 0; i < n; i++)
    {
        int si = blosum62char2int(subjectStr[i]);
        vTDia = right_shift_x_fill(arrT1+(k-1)*vec_length, 1, i==0?0:INT_MIN/2);
        vTUp = set_vector(m_striped, INT_MIN/2, GAP_OPEN_PENALTY);
        vTUp = add_vector(vTUp, vGapTUp);
        vF = set_vector(m_striped, GAP_OPEN_PENALTY+GAP_EXT_PENALTY*(i), GAP_EXT_PENALTY);
        vF = add_vector(vF, vGapF);
        vF = max_vector(vF, vTUp);
        for(int j = 0; j < k; j++)
        {
            vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
            vTLeft = add_array(arrT1+j*vec_length, vGapTLeft);
            vE = add_array(arrE+j*vec_length, vGapE);
            vE = max_vector(vE, vTLeft);
            store_vector(arrE+j*vec_length, vE);
            vT = max_vector(vTDia, vE, vF);
            store_vector(arrT2+j*vec_length, vT);
            // vMax = max_vector(vMax, vT);
            vTDia = load_vector(arrT1+j*vec_length);
            vTUp = vT;
            vTUp = add_vector(vTUp, vGapTUp);
            vF = add_vector(vF, vGapF);
            vF = max_vector(vF, vTUp);
        }
        vF = right_shift_x_fill(vF, 1, INT_MIN/2);
        int j = 0;
        vT = load_vector(arrT2+j*vec_length); 
        // float ratio = 0.f;
        // int count = 0;
        // while(count < vec_length-1)
        // while(count < 15)
        while(influence_test_max(vF, add_vector(vT, sub_vector(vGapTUp, vGapF))))
        {
            vT = max_vector(vF, vT);
            store_vector(arrT2+j*vec_length, vT);
            // vMax = max_vector(vMax, vT);
            vF = add_vector(vF, vGapF);
            if(++j >= k)
            {
                vF = right_shift_x_fill(vF, 1, INT_MIN/2);
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
    // score = reduce_max(vMax);
    score = arrT1[((m-1)%k)*w+(m-1)/k];
    _mm_free(arrT1);
    _mm_free(arrT2);
    _mm_free(arrE);
    _mm_free(profTDia);
    free(stripedQueryStr);

    return score;
}

int mod_striped_scan_NWAG(const char* queryStr, int m, const char* subjectStr, int n)
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
    for(int j = 0; j < k; j++) // Init
    {
        for(int i = 0; i < w; i++)
        {
            arrE[j * w + i] = GAP_OPEN_PENALTY + GAP_EXT_PENALTY * (i * k + j + 1 - 1);
        }
    }
    fill_n(arrT1, m_striped, INT_MIN/2); 

    vec vTDia, vTLeft, vTUp, vT;
    vec vE;
    // vec vMax = broadcast(INT_MIN);
    vec vGapTLeft = broadcast(GAP_OPEN_PENALTY); 
    vec vGapE = broadcast(GAP_EXT_PENALTY);
    // vec vConst = broadcast(0);
    for(int i = 0; i < n; i++)
    {
        int si = blosum62char2int(subjectStr[i]);
        vTDia = right_shift_x_fill(arrT1+(k-1)*vec_length, 1, i==0?0:INT_MIN/2);
        for(int j = 0; j < k; j++)
        {
            vTLeft = add_array(arrT1+j*vec_length, vGapTLeft);
            vE = add_array(arrE+j*vec_length, vGapE);
            vE = max_vector(vE, vTLeft);
            store_vector(arrE+j*vec_length, vE);
            vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
            vT = max_vector(vTDia, vE);
            store_vector(arrT2+j*vec_length, vT);
            vTDia = load_vector(arrT1+j*vec_length);
        }
        max_scan_striped_array(arrT2, arrScan, m_striped, GAP_OPEN_PENALTY+GAP_EXT_PENALTY*(i), GAP_EXT_PENALTY, GAP_OPEN_PENALTY);
        for(int j = 0; j < k; j++)
        {
            vTUp = load_vector(arrScan+j*vec_length);
            vT = load_vector(arrT2+j*vec_length);
            vT = max_vector(vT, vTUp);
            // vMax = max_vector(vMax, vT);
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

    // score = reduce_max(vMax);
    score = arrT1[((m-1)%k)*w+(m-1)/k];
    _mm_free(arrT1);
    _mm_free(arrT2);
    _mm_free(arrE);
    _mm_free(arrScan);
    _mm_free(profTDia);
    free(stripedQueryStr);

    return score;
}

int mod_striped_merged_NWAG(const char* queryStr, int m, const char* subjectStr, int n)
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
    for(int j = 0; j < k; j++) // Init
    {
        for(int i = 0; i < w; i++)
        {
            arrE[j * w + i] = GAP_OPEN_PENALTY + GAP_EXT_PENALTY * (i * k + j + 1 - 1);
        }
    }
    fill_n(arrT1, m_striped, INT_MIN/2); 

    vec vTDia, vTLeft, vTUp, vT;
    vec vE, vF;
    // vec vMax = broadcast(INT_MIN);
    vec vGapTLeft = broadcast(GAP_OPEN_PENALTY); 
    vec vGapTUp = broadcast(GAP_OPEN_PENALTY); 
    // vec vConst = broadcast(0);
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
            vTDia = right_shift_x_fill(arrT1+(k-1)*vec_length, 1, i==0?0:INT_MIN/2);
            vTUp = set_vector(m_striped, INT_MIN/2, GAP_OPEN_PENALTY);
            vTUp = add_vector(vTUp, vGapTUp);
            vF = set_vector(m_striped, GAP_OPEN_PENALTY+GAP_EXT_PENALTY*(i), GAP_EXT_PENALTY);
            vF = add_vector(vF, vGapF);
            vF = max_vector(vF, vTUp);
            for(int j = 0; j < k; j++)
            {
                vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
                vTLeft = add_array(arrT1+j*vec_length, vGapTLeft);
                vE = add_array(arrE+j*vec_length, vGapE);
                vE = max_vector(vE, vTLeft);
                store_vector(arrE+j*vec_length, vE);
                vT = max_vector(vTDia, vE, vF);
                store_vector(arrT2+j*vec_length, vT);
                // vMax = max_vector(vMax, vT);
                vTDia = load_vector(arrT1+j*vec_length);
                vTUp = vT;
                vTUp = add_vector(vTUp, vGapTUp);
                vF = add_vector(vF, vGapF);
                vF = max_vector(vF, vTUp);
            }
            vF = right_shift_x_fill(vF, 1, INT_MIN/2);
            int j = 0;
            vT = load_vector(arrT2+j*vec_length); 
            int count = 0;
            // while(count < vec_length-1)
            // while(count < 15)
            while(influence_test_max(vF, add_vector(vT, sub_vector(vGapTUp, vGapF))))
            {
                vT = max_vector(vF, vT);
                store_vector(arrT2+j*vec_length, vT);
                // vMax = max_vector(vMax, vT);
                vF = add_vector(vF, vGapF);
                if(++j >= k)
                {
                    vF = right_shift_x_fill(vF, 1, INT_MIN/2);
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
            vTDia = right_shift_x_fill(arrT1+(k-1)*vec_length, 1, i==0?0:INT_MIN/2);
            for(int j = 0; j < k; j++)
            {
                vTLeft = add_array(arrT1+j*vec_length, vGapTLeft);
                vE = add_array(arrE+j*vec_length, vGapE);
                vE = max_vector(vE, vTLeft);
                store_vector(arrE+j*vec_length, vE);
                vTDia = add_profile(profTDia+si*m_striped+j*vec_length, vTDia);
                vT = max_vector(vTDia, vE);
                store_vector(arrT2+j*vec_length, vT);
                vTDia = load_vector(arrT1+j*vec_length);
            }
            max_scan_striped_array(arrT2, arrScan, m_striped, GAP_OPEN_PENALTY+GAP_EXT_PENALTY*(i), GAP_EXT_PENALTY, GAP_OPEN_PENALTY);
            for(int j = 0; j < k; j++)
            {
                vTUp = load_vector(arrScan+j*vec_length);
                vT = load_vector(arrT2+j*vec_length);
                vT = max_vector(vT, vTUp);
                // vMax = max_vector(vMax, vT);
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
    // score = reduce_max(vMax);
    score = arrT1[((m-1)%k)*w+(m-1)/k];
    _mm_free(arrT1);
    _mm_free(arrT2);
    _mm_free(arrE);
    _mm_free(arrScan);
    _mm_free(profTDia);
    free(stripedQueryStr);

    return score;
}


