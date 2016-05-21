#include <iostream>
#include <climits>
#include "Config.h"
#include "Modules_int.h"
#include <omp.h> 
using namespace std;

// int mod_striped_iterate_SWAG(const char* queryStr, int m, const char* subjectStr, int n)
int mod_striped_iterate_SWAG(const char *stripedQueryStr, int m_striped, const char *subjectStr, int n, int8_t *profTDia, int m, int *buf)
{
    // Step 1: Stripe the query sequence
    int w = vec_length;
    int k = (m + w - 1) / w;
    /*
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
    */

    int score = 1;
    int *arrT1   = buf+0*m_striped; // (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    int *arrT2   = buf+1*m_striped; // (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length); 
    int *arrE    = buf+2*m_striped; // (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
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
    // _mm_free(arrT1);
    // _mm_free(arrT2);
    // _mm_free(arrE);
    // _mm_free(profTDia);
    // free(stripedQueryStr);

    return score;
}


// int mod_striped_merged_SWAG(const char* queryStr, int m, const char* subjectStr, int n)
int mod_striped_merged_SWAG(const char *stripedQueryStr, int m_striped, const char *subjectStr, int n, int8_t *profTDia, int m, int *buf)
{
    int probe_stride = m / 20;
    // Step 1: Stripe the query sequence
    int w = vec_length;
    int k = (m + w - 1) / w;
    /*
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
    */

    int score = 1;
    int *arrT1   = buf+0*m_striped; // (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    int *arrT2   = buf+1*m_striped; // (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length); 
    int *arrE    = buf+2*m_striped; // (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
    int *arrScan = buf+3*m_striped; // (int *)_mm_malloc(sizeof(int)*m_striped, 4*vec_length);
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
    // _mm_free(arrT1);
    // _mm_free(arrT2);
    // _mm_free(arrE);
    // _mm_free(arrScan);
    // _mm_free(profTDia);
    // free(stripedQueryStr);

    return score;
}

int mod_striped_merged_SWAG2(const char* queryStr, int m, const char* subjectStr, int n)
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
void align_db(const char *queryStr, int m, const char *dbStrs, int *dbOffset, int dbNum, int *output, int numThreads)
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

    int *buf = (int *)_mm_malloc(sizeof(int)*m_striped*numThreads*4, 4*vec_length);
    // cout << "size:" << m_striped*numThreads*4 << endl;

    int tid;
#pragma omp parallel for num_threads(numThreads) private(tid) default(shared) schedule(guided, 1)
    for(int i =0; i < dbNum; i++)
    {

		tid = omp_get_thread_num();
        int score = mod_striped_merged_SWAG(stripedQueryStr, m_striped, dbStrs+dbOffset[i], dbOffset[i+1]-dbOffset[i], profTDia, m, buf+tid*m_striped*4);
        // int score = mod_striped_iterate_SWAG(stripedQueryStr, m_striped, dbStrs+dbOffset[i], dbOffset[i+1]-dbOffset[i], profTDia, m, buf+tid*m_striped*4);
        // int score;
        // if(m < dbOffset[i+1]-dbOffset[i])
            // score = mod_striped_merged_SWAG2(queryStr, m, dbStrs+dbOffset[i], dbOffset[i+1]-dbOffset[i]);
        // else
            // score = mod_striped_merged_SWAG2(dbStrs+dbOffset[i], dbOffset[i+1]-dbOffset[i], queryStr, m);
        output[i] = score;
    }

    _mm_free(profTDia);
    _mm_free(buf);
    free(stripedQueryStr);
}
