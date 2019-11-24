// 编译使用 g++ -fopenmp -o kron_gen  filename.cpp -O3 -std=c++11
// 使用命令 ./kron_gen [点参数] [边参数] [文件名]
// 如 ./kron_gen 500 12 mykrongen.bin
// 上面将生成一个具有 500*(1024*1024) 顶点 12*(500*1024*1024) 个边的kron图数据

#include <iostream>
// #include <vector>
#include <algorithm>
#include <random>
#include <chrono>
// #include <array>
#include <cmath>
#include <fcntl.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#include <fstream>

using namespace std;

int main (int argc, char const **argv)
{
    // int const scaleFactor = atoi(argv[1]);
    int const nodeFactor = atoi(argv[1]);
    int const edgeFactor = atoi(argv[2]);
    const char* dt_name = argv[3];

    // ofstream fout(argv[3], ios::binary | ios::out | ios::trunc);

    // int64_t n_vert = (int64_t)pow(2.,scaleFactor);
    int64_t n_vert = nodeFactor*1024*1024;
    int64_t n_edge = n_vert * edgeFactor;

    int scaleFactor = (int)(log(n_vert)/log(2));
    cout << "Scale Factor is " << scaleFactor << endl;

    int64_t i, k, p;
    uint32_t *vertex = new uint32_t [2*n_edge];
    // vector< pair<uint32_t, uint32_t> > vertex(n_edge);
    #pragma omp parallel for
    for (i = 0; i < 2*n_edge; i++)
        vertex[i] = 1;
    // for (i = 0; i < n_edge; i++)
    //     vertex[i] = make_pair(1, 1);
    
    auto t0 = chrono::high_resolution_clock::now();
    
    double const a = 0.57;
    double const b = 0.19;
    double const c = 0.19;

    double a_b = a + b;
    double c_norm = c / (1-a_b);
    double a_norm = a / a_b;

    uint32_t kpow;
    
    int n_proc = omp_get_max_threads();
    int64_t block_size = n_edge/n_proc;

    for (k = 0; k < scaleFactor; k++)
    {
        kpow = (uint32_t)pow(2.,k);
        
        #pragma omp parallel for
        for (p = 0; p < n_proc; p++)
        {
            // uint32_t seed = chrono::system_clock::now().time_since_epoch().count();
            // auto tp = chrono::high_resolution_clock::now();
            auto t_seed = chrono::system_clock::now().time_since_epoch().count();

            // uint32_t seed = (chrono::duration_cast<milliseconds>(chrono::system_clock::now() - t0)) % 0xFFFF;
            default_random_engine e((t_seed.count()*(p+1)) % 0xFFFF);
            // cout << "rand seed is " << t_seed.count() + p*p << endl;
            uniform_real_distribution<double> u(0,1);

            bool h_bit, f_bit;
            int64_t z_ini = block_size*p;
            int64_t z_end = min(block_size*(p+1),n_edge);
        
            for (int64_t i = z_ini; i < z_end; i++)
            {
                h_bit = (u(e) > a_b);
                f_bit = (u(e) > (c_norm*h_bit + a_norm*(!h_bit)));
                
                vertex[2*i] += kpow*h_bit;
                vertex[2*i+1] += kpow*f_bit;
                // vertex[i].first += kpow*h_bit;
                // vertex[i].second += kpow*f_bit;
            }
        }
    }
    
    // vector <uint32_t> rand_perm (n_vert);
    uint32_t *rand_perm = new uint32_t [n_vert];
    #pragma omp parallel for
    for (i = 0; i < n_vert; i++)
        rand_perm[i] = i;

    uint32_t seed = chrono::system_clock::now().time_since_epoch().count();
    cout << seed << endl;

    auto tp = chrono::high_resolution_clock::now();
    auto t_seed = chrono::duration_cast<chrono::milliseconds>(tp.time_since_epoch());
    default_random_engine e(t_seed.count() % 0xFFFF);
    
    shuffle(rand_perm,rand_perm+n_vert,e);

    #pragma omp parallel for
    for (i = 0; i < 2*n_edge; i++)
        vertex[i] = rand_perm[vertex[i]];

    // #pragma omp parallel for
    // for (i = 0; i < n_edge; i++)
    // {
    //     vertex[i].first = rand_perm[vertex[i].first];
    //     vertex[i].second = rand_perm[vertex[i].second];
    // }
    
    uint64_t *vtx64 = (uint64_t *)vertex; 
    shuffle(vtx64, vtx64+n_edge, e);
    // shuffle(vertex.begin(),vertex.end(),e);

    FILE *fp = fopen(dt_name, "w+");
    fwrite(vtx64, sizeof(uint64_t), n_edge, fp);
    fclose(fp);

    delete [] vertex;
    delete [] rand_perm;
 
    // for (i = 0; i < n_edge; i++)
    // {
    //     fout.write((char*)&(vertex[i].first), sizeof(uint32_t));
    //     fout.write((char*)&(vertex[i].second), sizeof(uint32_t));
    // }

    // fout.close();
    auto tf = chrono::high_resolution_clock::now();
    auto t_used = chrono::duration_cast<chrono::seconds>(tf - t0);
    cout << "Data is generated successfully! Time elapses " << t_used.count() << " s." << endl;
    
    return 0;
}
