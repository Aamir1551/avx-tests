#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <glm.hpp>
#include <gtc/type_ptr.hpp>
#include <omp.h>

using namespace std::chrono;
using namespace std;

void AVX_MATMUL(float col0[4], float col1[4], float col2[4], float col3[4], float v[4]) {
    __m128 v0 = _mm_set1_ps(v[0]);
    __m128 v1 = _mm_set1_ps(v[1]);
    __m128 v2 = _mm_set1_ps(v[2]);
    __m128 v3 = _mm_set1_ps(v[3]);
    __m128 col0avx = _mm_loadu_ps(col0);
    __m128 col1avx = _mm_loadu_ps(col1);
    __m128 col2avx = _mm_loadu_ps(col2);
    __m128 col3avx = _mm_loadu_ps(col3);

    //__m128 v0 = _mm_permute_ps(v, 0);
    /*__m128 v1 = _mm_permute_ps(v, 0b01010101);
    __m128 v2 = _mm_permute_ps(v, 0b10101010);
    __m128 v3 = _mm_permute_ps(v, 0b11111111);*/
    __m128 prod0 = _mm_mul_ps(v0, col0avx);
    __m128 prod1 = _mm_mul_ps(v1, col1avx);
    __m128 s0 = _mm_add_ps(prod0, prod1);
    __m128 prod2 = _mm_mul_ps(v2, col2avx);
    __m128 prod3 = _mm_mul_ps(v3, col3avx);
    __m128 s1 = _mm_add_ps(prod2, prod3);
    __m128 out = _mm_add_ps(s0, s1);
}

void AVX_MATMUL4(__m128 &col0, __m128 &col1, __m128 &col2, __m128 &col3, __m128 &v) {

}


void GLMMatmul(glm::mat4 &m, glm::vec4 &v) {
    //glm::mat4 rr = glm::rotate(m, (float) 3.0, glm::vec3(1, 0, 0));
    glm::vec4 r = m * v;
    //glm::mat4 r = m * m;
}

void SerialMatMul4(float M[4][4], float X[4][4], float res[4][4]) {
    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            int s=0;
            for(int k=0; k<4; k++) {
                s += M[i][k] * X[k][j];
            }
            res[i][j] = s;
        }
    }
}

void SerialMatMul(float M[4][4], float v[4], float res[4]) {
    for(int i=0; i<4; i++) {
        for(int j=0; j<1; j++) {
            int s=0;
            for(int k=0; k<4; k++) {
                s += M[i][k] * v[k];
            }
            res[i] = s;
        }
    }
}

int main() {

    //float *res = (float *) &out;
    //printf("%lf %lf %lf %lf\n", res[0], res[1], res[2], res[3]);

    // TODO:
    // 4x4 matrix multiplications

    int n = 50000000;


    //int n = 10000000;

    float M[4][4];
    float X[4][4];
    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            M[i][j] = rand();
            X[i][j] = rand();
        }
    }
    float v[4];
    for(int i=0; i<4; i++) {
        v[i] = rand();
    }
    float res[4];
    float res4[4][4];
    auto start1 = high_resolution_clock::now();
//#pragma omp parallel for
    for(int i=0; i<n; i++) {
        //SerialMatMul(M, v, res);
        SerialMatMul4(M, X, res4);
    }
    auto end1 = high_resolution_clock::now();
    auto duration1 = duration_cast<microseconds>(end1 - start1);
    cout << "Time Taken For Serial: " << duration1.count() << endl;



    float m_v[16];
    for(int i=0; i<16; i++) {
        m_v[i] = rand() % 4;
    }
    glm::mat4 m = glm::make_mat4(m_v);
    glm::vec4 v_glm = glm::vec4(rand(), rand(), rand(), rand());
    auto start2 = high_resolution_clock::now();
//#pragma omp parallel for
    for(int i=0; i<n; i++) {
        GLMMatmul(m, v_glm);
    }
    auto end2 = high_resolution_clock::now();
    auto duration2 = duration_cast<microseconds>(end2 - start2);
    cout << "Time Taken For GLM:    " << duration2.count() << endl;

    /*__m128 v_avx = _mm_setr_ps(rand(), rand(), rand(), rand());
    __m128 col0 = _mm_setr_ps(rand(), rand(), rand(), rand());
    __m128 col1 = _mm_setr_ps(rand(), rand(), rand(), rand());
    __m128 col2 = _mm_setr_ps(rand(), rand(), rand(), rand());
    __m128 col3 = _mm_setr_ps(rand(), rand(), rand(), rand());*/
    float col0[4] = {(float) rand(), (float) rand(), (float) rand(), (float) rand()};
    float col1[4] = {(float) rand(), (float) rand(), (float) rand(), (float) rand()};
    float col2[4] = {(float) rand(), (float) rand(), (float) rand(), (float) rand()};
    float col3[4] = {(float) rand(), (float) rand(), (float) rand(), (float) rand()};
    float vec[4] = {(float) rand(), (float) rand(), (float) rand(), (float) rand()};
    auto start = high_resolution_clock::now();
//#pragma omp parallel for
    for(int i=0; i<n * 4; i++) {
        AVX_MATMUL(col0, col1, col2, col3, vec);
    }
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "Time Taken For AVX:    " << duration.count() << endl;

}