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


#ifndef _MODULES_H
#define _MODULES_H

#include <iostream>
#include <algorithm>
#include <immintrin.h>
#include <cstdarg>
#include <cstdint>

#define max(a,b) ((a)>(b)?(a):(b))

#ifdef __MIC__
typedef __m512i vec;
static int vec_length = 16;
inline __m512i right_shift_x_fill(int *array, int num, int v) 
{
    __m512i v_ret;
    __m512i cv_rev = _mm512_set_epi32(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
    unsigned short mask = 0xffff;
    mask <<= num;
    __m512i cv_fil = _mm512_set1_epi32(v);
    
    v_ret = _mm512_load_epi32(array);
    v_ret = _mm512_permutevar_epi32(cv_rev, v_ret);
    // v_ret = _mm512_alignr_epi32(cv_fil, v_ret, num); 
    // v_ret = _mm512_permutevar_epi32(cv_rev, v_ret);
    v_ret = _mm512_mask_swizzle_epi32(cv_fil, mask, v_ret, _MM_SWIZ_REG_NONE);

    return v_ret;
}
inline __m512i right_shift_x_fill(__m512i v_ret, int num, int v) 
{
    __m512i cv_rev = _mm512_set_epi32(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
    unsigned short mask = 0xffff;
    mask <<= num;
    __m512i cv_fil = _mm512_set1_epi32(v);
    
    v_ret = _mm512_permutevar_epi32(cv_rev, v_ret);
    v_ret = _mm512_mask_swizzle_epi32(cv_fil, mask, v_ret, _MM_SWIZ_REG_NONE);

    return v_ret;
}
inline __m512i set_vector(int m_striped, int init, int gap)
{
    int k = m_striped / vec_length;
    __m512i v = _mm512_set_epi32(
            init + gap * k * 15,
            init + gap * k * 14,
            init + gap * k * 13,
            init + gap * k * 12,
            init + gap * k * 11,
            init + gap * k * 10,
            init + gap * k * 9 ,
            init + gap * k * 8 ,
            init + gap * k * 7 ,
            init + gap * k * 6 ,
            init + gap * k * 5 ,
            init + gap * k * 4 ,
            init + gap * k * 3 ,
            init + gap * k * 2 ,
            init + gap * k * 1 ,
            init); 
    return v;
}
inline __m512i add_vector(__m512i va, __m512i vb)
{
    return _mm512_add_epi32(va, vb);
}
inline __m512i sub_vector(__m512i va, __m512i vb)
{
    return _mm512_sub_epi32(va, vb);
}
inline __m512i max_vector(__m512i v1, __m512i v2) 
{
    return _mm512_max_epi32(v1, v2);
}
inline __m512i max_vector(__m512i v1, __m512i v2, __m512i v3) 
{
    __m512i v = _mm512_max_epi32(v1, v2);
    v = _mm512_max_epi32(v, v3);
    return v;
}
inline __m512i max_vector(__m512i v1, __m512i v2, __m512i v3, __m512i v4) 
{
    __m512i v = _mm512_max_epi32(v1, v2);
    v = _mm512_max_epi32(v, v3);
    v = _mm512_max_epi32(v, v4);
    return v;
}
inline __m512i add_profile(int8_t *profile, __m512i v) 
{
    __m512i v_prf = _mm512_extload_epi32(profile, 
            _MM_UPCONV_EPI32_SINT8, _MM_BROADCAST32_NONE, _MM_HINT_NONE); 
    v = _mm512_add_epi32(v, v_prf);
    return v;
}
inline void store_vector(int *addr, __m512i v)
{
    _mm512_store_epi32(addr, v);
}
inline __m512i load_vector(int *addr)
{
    return _mm512_load_epi32(addr);
}
inline bool influence_test_max(__m512i v_a, __m512i v_b)
{
    unsigned short mask = _mm512_cmpgt_epi32_mask(v_a, v_b);
    return mask > 0x00;
}
inline int reduce_max(__m512i v)
{
    return  _mm512_reduce_max_epi32(v);
}

inline __m512i broadcast(int val)
{
    return  _mm512_set1_epi32(val);
}

inline __m512i add_array(int *array, __m512i v) 
{
    __m512i v_arr = _mm512_load_epi32(array);
    v = _mm512_add_epi32(v_arr, v);
    return v;
}
inline void max_scan_striped_array(int *inp, int *out, int m_striped, int init, int gapB, int gapC)
{
    int k = m_striped / vec_length;
    __m512i vecGap  = _mm512_set1_epi32(gapB);
    __m512i vecGapC  = _mm512_set1_epi32(gapC);
    __m512i vecGapPureC  = _mm512_set1_epi32(gapC-gapB);
    __m512i vecScan = _mm512_load_epi32(inp);
            vecScan = _mm512_add_epi32(vecScan, vecGap);
    for(int i = 1; i < k; i++)
    {
        __m512i vecH = _mm512_load_epi32(inp+i*vec_length);
                vecH = _mm512_add_epi32(vecH, vecGap);
        vecScan = _mm512_max_epi32(vecH, _mm512_add_epi32(vecScan, vecGap));
    }

    int scan0  = init + gapB;
    int scan1  = max(scan0 +k*gapB, ((int *)&vecScan)[0 ]);
    int scan2  = max(scan1 +k*gapB, ((int *)&vecScan)[1 ]);
    int scan3  = max(scan2 +k*gapB, ((int *)&vecScan)[2 ]);
    int scan4  = max(scan3 +k*gapB, ((int *)&vecScan)[3 ]);
    int scan5  = max(scan4 +k*gapB, ((int *)&vecScan)[4 ]);
    int scan6  = max(scan5 +k*gapB, ((int *)&vecScan)[5 ]);
    int scan7  = max(scan6 +k*gapB, ((int *)&vecScan)[6 ]);
    int scan8  = max(scan7 +k*gapB, ((int *)&vecScan)[7 ]);
    int scan9  = max(scan8 +k*gapB, ((int *)&vecScan)[8 ]);
    int scan10 = max(scan9 +k*gapB, ((int *)&vecScan)[9 ]);
    int scan11 = max(scan10+k*gapB, ((int *)&vecScan)[10]);
    int scan12 = max(scan11+k*gapB, ((int *)&vecScan)[11]);
    int scan13 = max(scan12+k*gapB, ((int *)&vecScan)[12]);
    int scan14 = max(scan13+k*gapB, ((int *)&vecScan)[13]);
    int scan15 = max(scan14+k*gapB, ((int *)&vecScan)[14]);

    __m512i vecLast = _mm512_set_epi32(
                scan15,
                scan14,
                scan13,
                scan12,
                scan11,
                scan10,
                scan9 ,
                scan8 ,
                scan7 ,
                scan6 ,
                scan5 ,
                scan4 ,
                scan3 ,
                scan2 ,
                scan1 ,
                scan0 );

    if(gapC!=gapB)
    {
        _mm512_store_epi32(out, _mm512_add_epi32(vecLast, vecGapPureC));
        for(int i = 1; i < k; i++)
        {
            __m512i vecH = _mm512_load_epi32(inp+(i-1)*vec_length);
                    vecH = _mm512_add_epi32(vecH, vecGap);
            vecLast= _mm512_max_epi32(vecH, _mm512_add_epi32(vecLast, vecGap));
            _mm512_store_epi32(out+i*vec_length, _mm512_add_epi32(vecLast, vecGapPureC));
        }
    }
    else
    {
        _mm512_store_epi32(out, vecLast);
        for(int i = 1; i < k; i++)
        {
            __m512i vecH = _mm512_load_epi32(inp+(i-1)*vec_length);
                    vecH = _mm512_add_epi32(vecH, vecGap);
            vecLast= _mm512_max_epi32(vecH, _mm512_add_epi32(vecLast, vecGap));
            _mm512_store_epi32(out+i*vec_length, vecLast);
        }
    }
}
#endif
#ifdef __AVX2__
typedef __m256i vec;
// static int vec_length = 8;
// static int vec_length = 16;
static int vec_length = 32;
inline __m256i right_shift_x_fill(int8_t *array, int num, int v) 
{
    __m256i v_ret;
    __m256i cv_rev = _mm256_set_epi8(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
    v_ret = _mm256_load_si256((__m256i *)array);
    
    v_ret = _mm256_shuffle_epi8 (v_ret, cv_rev); 
    __int8 fill = _mm256_extract_epi8(v_ret, 0);
    v_ret = _mm256_insert_epi8(v_ret, fill, 16);
    v_ret = _mm256_insert_epi8(v_ret, v, 0);
    return v_ret;
}

inline __m256i right_shift_x_fill(__m256i v_ret, int num, int v) 
{
    __m256i cv_rev = _mm256_set_epi8(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
    
    v_ret = _mm256_shuffle_epi8 (v_ret, cv_rev); 
    __int8 fill = _mm256_extract_epi8(v_ret, 0);
    v_ret = _mm256_insert_epi8(v_ret, fill, 16);
    v_ret = _mm256_insert_epi8(v_ret, v, 0);
    return v_ret;
}
inline __m256i set_vector(int m_striped, int8_t init, int8_t gap)
{
    // int k = m_striped / vec_length;
    // __m256i v = _mm256_set_epi32(
            // init + gap * k * 7 ,
            // init + gap * k * 6 ,
            // init + gap * k * 5 ,
            // init + gap * k * 4 ,
            // init + gap * k * 3 ,
            // init + gap * k * 2 ,
            // init + gap * k * 1 ,
            // init); 
    int8_t k = m_striped / vec_length;
    __m256i v = _mm256_set_epi8(
            init /*+ gap * k * 31*/,
            init /*+ gap * k * 30*/,
            init /*+ gap * k * 29*/,
            init /*+ gap * k * 28*/,
            init /*+ gap * k * 27*/,
            init /*+ gap * k * 26*/,
            init /*+ gap * k * 25*/,
            init /*+ gap * k * 24*/,
            init /*+ gap * k * 23*/,
            init /*+ gap * k * 22*/,
            init /*+ gap * k * 21*/,
            init /*+ gap * k * 20*/,
            init /*+ gap * k * 19*/,
            init /*+ gap * k * 18*/,
            init /*+ gap * k * 17*/,
            init /*+ gap * k * 16*/,
            init /*+ gap * k * 15*/,
            init /*+ gap * k * 14*/,
            init /*+ gap * k * 13*/,
            init /*+ gap * k * 12*/,
            init /*+ gap * k * 11*/,
            init /*+ gap * k * 10*/,
            init /*+ gap * k * 9 */,
            init /*+ gap * k * 8 */,
            init /*+ gap * k * 7 */,
            init /*+ gap * k * 6 */,
            init /*+ gap * k * 5 */,
            init /*+ gap * k * 4 */,
            init /*+ gap * k * 3 */,
            init /*+ gap * k * 2 */,
            init /*+ gap * k * 1 */,
            init); 
    return v;
}
inline bool influence_test_max(__m256i v_a, __m256i v_b)
{
    // __m256i mask = _mm256_cmpgt_epi32(v_a, v_b);
    __m256i mask = _mm256_cmpgt_epi8(v_a, v_b);
    __m128i l = _mm256_extractf128_si256(mask, 0);
    __m128i h = _mm256_extractf128_si256(mask, 1);
    bool lflag = (bool)_mm_test_all_zeros(l, _mm_set1_epi32(0xffffffff));
    bool hflag = (bool)_mm_test_all_zeros(h, _mm_set1_epi32(0xffffffff));
    return !lflag || !hflag;
}
inline bool check_overflow(__m256i va, __m256i vb)
{
    __m256i v_zero = _mm256_set1_epi8(0);
    __m256i v_max = _mm256_set1_epi8(CHAR_MAX);
    __m256i v_min = _mm256_set1_epi8(CHAR_MIN);
    __m256i mask = _mm256_cmpgt_epi8(vb, v_zero);
    // __m256i v_pos = _mm256_blendv_epi8(v_max, vb, mask);
    __m256i v_minus = _mm256_sub_epi8(v_max, vb);
    __m256i v_pos = _mm256_blendv_epi8(v_max, v_minus, mask);
    bool res1 = influence_test_max(va, v_pos);
    // std::cout << res1 << std::endl;
    if(res1)
        return true;
    mask = _mm256_cmpgt_epi8(v_zero, vb);
    __m256i v_plus = _mm256_sub_epi8(v_min, vb);
    __m256i v_neg = _mm256_blendv_epi8(v_min, v_plus, mask);
    bool res2 = influence_test_max(v_neg, va);
    // std::cout << res2 << std::endl;
    if(res2)
        return true;

    // for(int i=0;i<32;i++)
        // std::cout << (int)((int8_t*)(&v_neg))[i] << " ";
    // std::cout << std::endl;
    return false;
}

inline __m256i add_vector(__m256i va, __m256i vb, bool *check)
{
    if(check_overflow(va, vb))
    {
        check[0] = false;
    }
    return _mm256_add_epi8(va, vb);
}
inline __m256i sub_vector(__m256i va, __m256i vb)
{
    return _mm256_sub_epi8(va, vb);
}
inline __m256i max_vector(__m256i v1, __m256i v2) 
{
    return _mm256_max_epi8(v1, v2);
}
inline __m256i max_vector(__m256i v1, __m256i v2, __m256i v3) 
{
    __m256i va = _mm256_max_epi8(v1, v2);
    __m256i v  = _mm256_max_epi8(va, v3);
    return v;
}
inline __m256i max_vector(__m256i v1, __m256i v2, __m256i v3, __m256i v4) 
{
    __m256i va = _mm256_max_epi8(v1, v2);
    __m256i vb = _mm256_max_epi8(v3, v4);
    __m256i v  = _mm256_max_epi8(va, vb);
    return v;
}
inline __m256i add_profile(int8_t *profile, __m256i v, bool *check) 
{
    // __m128i a = _mm_lddqu_si128((__m128i*)(profile));
    // __m128i b = _mm_cvtepi8_epi32(a);
    // __m256i veci = _mm256_insertf128_si256(veci, b, 0);
    // a = _mm_shuffle_epi32(a, 0xB1);
    // b = _mm_cvtepi8_epi32(a);
    // __m256i v_prf = _mm256_insertf128_si256(veci, b, 1); 
    // __m256i v_prf = _mm256_cvtepi8_epi16 (a);
    __m256i v_prf = _mm256_load_si256((__m256i *)profile);
    if(check_overflow(v, v_prf))
    {
        check[0] = false;
    }
    v = _mm256_add_epi8(v, v_prf);
    return v;
}
inline void store_vector(int8_t *addr, __m256i v)
{
    _mm256_store_si256((__m256i *)(addr), v);
}
inline __m256i load_vector(int8_t *addr)
{
    return _mm256_load_si256((__m256i *)addr);

    // __m128i a = _mm_lddqu_si128((__m128i*)addr);
    // __m128i b = _mm_cvtepi16_epi32(a);
    // __m256i v_ret = _mm256_insertf128_si256(v_ret, b, 0);
    // a = _mm_shuffle_epi32(a, 0x4E);
    // b = _mm_cvtepi16_epi32(a);
    // v_ret = _mm256_insertf128_si256(v_ret, b, 1);
    // return v_ret;
}

inline int8_t reduce_max(__m256i v)
{
    return  *std::max_element((int8_t *)&v, ((int8_t *)&v) + vec_length);
}
inline __m256i broadcast(int8_t val)
{
    return  _mm256_set1_epi8(val);
}
inline __m256i add_array(int8_t *array, __m256i v, bool *check) 
{
    __m256i v_arr = _mm256_load_si256((__m256i *)array);

    if(check_overflow(v_arr, v))
    {
        check[0] = false;
    }
    v = _mm256_add_epi8(v_arr, v);
    return v;
}
/*
inline void max_scan_striped_array(int16_t *inp, int16_t *out, int m_striped, int init, int gapB, int gapC)
{
    int k = m_striped / vec_length;
    __m256i vecGap  = _mm256_set1_epi16(gapB);
    __m256i vecGapC  = _mm256_set1_epi16(gapC);
    __m256i vecGapPureC  = _mm256_set1_epi16(gapC-gapB);
    // __m256i vecScan = _mm256_load_si256((__m256i *)inp);
    __m256i vecScan = load_vector(inp);
            vecScan = _mm256_add_epi16(vecScan, vecGap);
    for(int i = 1; i < k; i++)
    {
        // __m256i vecH = _mm256_load_si256((__m256i *)(inp+i*vec_length));
        __m256i vecH = load_vector(inp+i*vec_length);
                vecH = _mm256_add_epi16(vecH, vecGap);
        vecScan = _mm256_max_epi16(vecH, _mm256_add_epi16(vecScan, vecGap));
    }

    int16_t scan0  = init + gapB;
    int16_t scan1  = max(scan0 +k*gapB, ((int16_t *)&vecScan)[0 ]);
    int16_t scan2  = max(scan1 +k*gapB, ((int16_t *)&vecScan)[1 ]);
    int16_t scan3  = max(scan2 +k*gapB, ((int16_t *)&vecScan)[2 ]);
    int16_t scan4  = max(scan3 +k*gapB, ((int16_t *)&vecScan)[3 ]);
    int16_t scan5  = max(scan4 +k*gapB, ((int16_t *)&vecScan)[4 ]);
    int16_t scan6  = max(scan5 +k*gapB, ((int16_t *)&vecScan)[5 ]);
    int16_t scan7  = max(scan6 +k*gapB, ((int16_t *)&vecScan)[6 ]);
    int16_t scan8  = max(scan7 +k*gapB, ((int16_t *)&vecScan)[7 ]);
    int16_t scan9  = max(scan8 +k*gapB, ((int16_t *)&vecScan)[8 ]);
    int16_t scan10 = max(scan9 +k*gapB, ((int16_t *)&vecScan)[9 ]);
    int16_t scan11 = max(scan10+k*gapB, ((int16_t *)&vecScan)[10]);
    int16_t scan12 = max(scan11+k*gapB, ((int16_t *)&vecScan)[11]);
    int16_t scan13 = max(scan12+k*gapB, ((int16_t *)&vecScan)[12]);
    int16_t scan14 = max(scan13+k*gapB, ((int16_t *)&vecScan)[13]);
    int16_t scan15 = max(scan14+k*gapB, ((int16_t *)&vecScan)[14]);

    __m256i vecLast = _mm256_set_epi16(
                scan15,
                scan14,
                scan13,
                scan12,
                scan11,
                scan10,
                scan9 ,
                scan8 ,
                scan7 ,
                scan6 ,
                scan5 ,
                scan4 ,
                scan3 ,
                scan2 ,
                scan1 ,
                scan0 );

    if(gapC!=gapB)
    {
        // _mm256_store_si256((__m256i *)out, _mm256_add_epi32(vecLast, vecGapPureC));
        store_vector(out, _mm256_add_epi16(vecLast, vecGapPureC));
        for(int i = 1; i < k; i++)
        {
            // __m256i vecH = _mm256_load_si256((__m256i *)(inp+(i-1)*vec_length));
            __m256i vecH = load_vector(inp+(i-1)*vec_length);
                    vecH = _mm256_add_epi16(vecH, vecGap);
            vecLast= _mm256_max_epi16(vecH, _mm256_add_epi16(vecLast, vecGap));
            // _mm256_store_si256((__m256i *)(out+i*vec_length), _mm256_add_epi32(vecLast, vecGapPureC));
            store_vector(out+i*vec_length, _mm256_add_epi16(vecLast, vecGapPureC));
        }
    }
    else
    {
        // _mm256_store_si256((__m256i *)out, vecLast);
        store_vector(out, vecLast);
        for(int i = 1; i < k; i++)
        {
            // __m256i vecH = _mm256_load_si256((__m256i *)(inp+(i-1)*vec_length));
            __m256i vecH = load_vector(inp+(i-1)*vec_length);
                    vecH = _mm256_add_epi16(vecH, vecGap);
            vecLast= _mm256_max_epi16(vecH, _mm256_add_epi16(vecLast, vecGap));
            // _mm256_store_si256((__m256i *)(out+i*vec_length), vecLast);
            store_vector(out+i*vec_length, vecLast);
        }
    }
    
}
*/
#endif

#endif
