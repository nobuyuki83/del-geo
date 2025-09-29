#pragma once
#include <cuda/std/array>

namespace vec3 {

__device__
auto norm(const float* a) -> float
{
   const float sqnrm = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
   return sqrt(sqnrm);
}

__device__
auto add(const float* a, const float* b) -> cuda::std::array<float,3>
{
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

__device__
void add_in_place(float* a, const float* b)
{
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
}

__device__
auto sub(const float* a, const float* b) -> cuda::std::array<float,3>
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

__device__
auto cross(const float* v1, const float* v2) -> cuda::std::array<float,3>
{
    return {
        v1[1] * v2[2] - v2[1] * v1[2],
        v1[2] * v2[0] - v2[2] * v1[0],
        v1[0] * v2[1] - v2[0] * v1[1],
    };
}

__device__
auto dot(const float* a, const float* b) -> float
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

__device__
auto axpy(float a, const float* x, const float* y) -> cuda::std::array<float,3>
{
    return {
        a * x[0] + y[0],
        a * x[1] + y[1],
        a * x[2] + y[2],
    };
}

__device__
auto scale(const float* x, float a) -> cuda::std::array<float,3>
{
    return {
        a * x[0],
        a * x[1],
        a * x[2],
    };
}

__device__
void add_inplace(float* x, const float* y) {
    x[0] += y[0];
    x[1] += y[1];
    x[2] += y[2];
}

__device__
auto mult_mat3_col_major(const float* a, const float* m) -> cuda::std::array<float,3>
{
    return {
        a[0] * m[0] + a[1] * m[1] + a[2] * m[2],
        a[0] * m[3] + a[1] * m[4] + a[2] * m[5],
        a[0] * m[6] + a[1] * m[7] + a[2] * m[8],
    };
}


}