#include <iostream>
#include <fstream>
#include <immintrin.h>
#include "Config.h"
#include <string>
#include <algorithm>
#include <vector>

using namespace std;

int blosum62char2int(char base)
{
    switch(base)
    {
        case 'A': return 0 ;
        case 'R': return 1 ;
        case 'N': return 2 ;
        case 'D': return 3 ;
        case 'C': return 4 ;
        case 'Q': return 5 ;
        case 'E': return 6 ;
        case 'G': return 7 ;
        case 'H': return 8 ;
        case 'I': return 9 ;
        case 'L': return 10;
        case 'K': return 11;
        case 'M': return 12;
        case 'F': return 13;
        case 'P': return 14;
        case 'S': return 15;
        case 'T': return 16;
        case 'W': return 17;
        case 'Y': return 18;
        case 'V': return 19;
        case 'B': return 20;
        case 'J': return 21;
        case 'Z': return 22;
        case 'X': return 23;
        case 'U': return 24;
        case '-': return 24;
    }
    fprintf(stderr, "Error: An invalid character has been detected! ----( %c )----\n", base);
    exit(EXIT_FAILURE);
}

char blosum62int2char(int index)
{
    switch(index)
    {
        case 0 : return 'A';
        case 1 : return 'R';
        case 2 : return 'N';
        case 3 : return 'D';
        case 4 : return 'C';
        case 5 : return 'Q';
        case 6 : return 'E';
        case 7 : return 'G';
        case 8 : return 'H';
        case 9 : return 'I';
        case 10: return 'L';
        case 11: return 'K';
        case 12: return 'M';
        case 13: return 'F';
        case 14: return 'P';
        case 15: return 'S';
        case 16: return 'T';
        case 17: return 'W';
        case 18: return 'Y';
        case 19: return 'V';
        case 20: return 'B';
        case 21: return 'J';
        case 22: return 'Z';
        case 23: return 'X';
        case 24: return '-';
    }
    fprintf(stderr, "Error: An invalid index has been detected! ----( %d )----\n", index);
    exit(EXIT_FAILURE);
}

void striped_vec_scan_max(int *input, int *output, int size, int _gap, int row_i)
{
    int segLen = size/AVX3_WIDTH;
    __m512i vecGap  = _mm512_set1_epi32(_gap);
    __m512i vecScan = _mm512_load_epi32(input);
            vecScan = _mm512_add_epi32(vecScan, vecGap);
    for(int i = 1; i < segLen; i++)
    {
        __m512i vecH = _mm512_load_epi32(input+i*AVX3_WIDTH);
                vecH = _mm512_add_epi32(vecH, vecGap);
        vecScan = _mm512_max_epi32(vecH, _mm512_add_epi32(vecScan, vecGap));
    }

    int scan0  = _gap*(row_i+1)+_gap;
    int scan1  = max(scan0 +segLen*_gap, ((int *)&vecScan)[0 ]);
    int scan2  = max(scan1 +segLen*_gap, ((int *)&vecScan)[1 ]);
    int scan3  = max(scan2 +segLen*_gap, ((int *)&vecScan)[2 ]);
    int scan4  = max(scan3 +segLen*_gap, ((int *)&vecScan)[3 ]);
    int scan5  = max(scan4 +segLen*_gap, ((int *)&vecScan)[4 ]);
    int scan6  = max(scan5 +segLen*_gap, ((int *)&vecScan)[5 ]);
    int scan7  = max(scan6 +segLen*_gap, ((int *)&vecScan)[6 ]);
    int scan8  = max(scan7 +segLen*_gap, ((int *)&vecScan)[7 ]);
    int scan9  = max(scan8 +segLen*_gap, ((int *)&vecScan)[8 ]);
    int scan10 = max(scan9 +segLen*_gap, ((int *)&vecScan)[9 ]);
    int scan11 = max(scan10+segLen*_gap, ((int *)&vecScan)[10]);
    int scan12 = max(scan11+segLen*_gap, ((int *)&vecScan)[11]);
    int scan13 = max(scan12+segLen*_gap, ((int *)&vecScan)[12]);
    int scan14 = max(scan13+segLen*_gap, ((int *)&vecScan)[13]);
    int scan15 = max(scan14+segLen*_gap, ((int *)&vecScan)[14]);

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

    _mm512_store_epi32(output, vecLast);
    for(int i = 1; i < segLen; i++)
    {
        __m512i vecH = _mm512_load_epi32(input+(i-1)*AVX3_WIDTH);
                vecH = _mm512_add_epi32(vecH, vecGap);
        vecLast= _mm512_max_epi32(vecH, _mm512_add_epi32(vecLast, vecGap));
        _mm512_store_epi32(output+i*AVX3_WIDTH, vecLast);
    }
}


void striped_vec_scan_min(int *input, int *output, int size, int _gap, int row_i)
{
    int segLen = size/AVX3_WIDTH;
    __m512i vecGap  = _mm512_set1_epi32(_gap);
    __m512i vecScan = _mm512_load_epi32(input);
            vecScan = _mm512_add_epi32(vecScan, vecGap);
    for(int i = 1; i < segLen; i++)
    {
        __m512i vecH = _mm512_load_epi32(input+i*AVX3_WIDTH);
                vecH = _mm512_add_epi32(vecH, vecGap);
        vecScan = _mm512_min_epi32(vecH, _mm512_add_epi32(vecScan, vecGap));
    }

    int scan0  = _gap*(row_i+1)+_gap;
    int scan1  = min(scan0 +segLen*_gap, ((int *)&vecScan)[0 ]);
    int scan2  = min(scan1 +segLen*_gap, ((int *)&vecScan)[1 ]);
    int scan3  = min(scan2 +segLen*_gap, ((int *)&vecScan)[2 ]);
    int scan4  = min(scan3 +segLen*_gap, ((int *)&vecScan)[3 ]);
    int scan5  = min(scan4 +segLen*_gap, ((int *)&vecScan)[4 ]);
    int scan6  = min(scan5 +segLen*_gap, ((int *)&vecScan)[5 ]);
    int scan7  = min(scan6 +segLen*_gap, ((int *)&vecScan)[6 ]);
    int scan8  = min(scan7 +segLen*_gap, ((int *)&vecScan)[7 ]);
    int scan9  = min(scan8 +segLen*_gap, ((int *)&vecScan)[8 ]);
    int scan10 = min(scan9 +segLen*_gap, ((int *)&vecScan)[9 ]);
    int scan11 = min(scan10+segLen*_gap, ((int *)&vecScan)[10]);
    int scan12 = min(scan11+segLen*_gap, ((int *)&vecScan)[11]);
    int scan13 = min(scan12+segLen*_gap, ((int *)&vecScan)[12]);
    int scan14 = min(scan13+segLen*_gap, ((int *)&vecScan)[13]);
    int scan15 = min(scan14+segLen*_gap, ((int *)&vecScan)[14]);

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

    _mm512_store_epi32(output, vecLast);
    for(int i = 1; i < segLen; i++)
    {
        __m512i vecH = _mm512_load_epi32(input+(i-1)*AVX3_WIDTH);
                vecH = _mm512_add_epi32(vecH, vecGap);
        vecLast= _mm512_min_epi32(vecH, _mm512_add_epi32(vecLast, vecGap));
        _mm512_store_epi32(output+i*AVX3_WIDTH, vecLast);
    }
}

#define BUF_SIZE 1024
string read_file(char *filename)
{
    ifstream fstream_q;
    fstream_q.open(filename); 
    if(!fstream_q)
    {
        std::cerr << "Cannot open query file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    string seq_str;
    char* line_buf = new char[BUF_SIZE];
    while(!fstream_q.eof())
    {
        fstream_q.getline(line_buf, BUF_SIZE);
        if(line_buf[0] != '>')
        {
            seq_str += line_buf;
        }
    }
    delete[] line_buf;
    return seq_str;
}

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1].length() < v[i2].length();});

  return idx;
}

std::vector<std::string> read_db_file(char *filename, char *&dbStrs, int *&dbOffset, int &dbNum, std::vector<size_t> &ind)
{
    std::ifstream fstrm;
    fstrm.open(filename); 
    if(!fstrm)
    {
        std::cerr << "Cannot open db file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> pool;
    std::vector<std::string> name;
    std::string *seq_str = NULL;
    std::string *name_str = NULL;
    char* line_buf = new char[BUF_SIZE];
    int total = 0;
    while(!fstrm.eof())
    {
        fstrm.getline(line_buf, BUF_SIZE);
        if(line_buf[0] != '>')
        {
            *seq_str += line_buf;
        }else
        {
            name_str = new std::string();
            *name_str += line_buf;
            name.push_back((*name_str).substr(1, (*name_str).find_first_of(" ")-1));
            // cout << name[name.size()-1] << endl;
            if(seq_str != NULL)
            {
                pool.push_back(*seq_str);
                total += (*seq_str).length();
            }
            seq_str = new std::string();
        }
    }
    pool.push_back(*seq_str);
    total += (*seq_str).length();
    delete[] line_buf;

    struct compare 
    {
        bool operator()(const std::string& first, const std::string& second) 
        {
            return first.length() < second.length();
        }
    } cmpFunc;   
    ind = sort_indexes(pool);
    // for (auto i: ind) 
    // {
        // cout << i << endl;
    // }
    std::sort(pool.begin(), pool.end(), cmpFunc);
    // for(int i=0; i< pool.size(); i++)
        // std::cout << pool[i].length() << std::endl;
    dbStrs = (char *)malloc(sizeof(char)*total);
    dbOffset = (int *)malloc(sizeof(int)*(pool.size()+1));
    dbOffset[0] = 0;
    for(int i=0; i< pool.size(); i++)
    {
        dbOffset[i+1] = dbOffset[i] + pool[i].length();
        std::copy(pool[i].begin(), pool[i].end(), dbStrs+dbOffset[i]);
        // std::cout << dbOffset[i+1] << std::endl;
    }
    dbNum = pool.size();

    // std::cout << dbNum << std::endl;
    return name;
}

