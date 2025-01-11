#pragma once
#include <cuda/std/array>
#include "vec3.h"
#include "mat3_col_major.h"

namespace mat4_col_major {

__device__
auto transform_homogeneous(
  const float* transform,
  const float* x) -> cuda::std::array<float,3>
{
    float y3 = transform[3] * x[0] + transform[7] * x[1] + transform[11] * x[2] + transform[15];
    const float y0 = transform[0] * x[0] + transform[4] * x[1] + transform[8] * x[2] + transform[12];
    const float y1 = transform[1] * x[0] + transform[5] * x[1] + transform[9] * x[2] + transform[13];
    const float y2 = transform[2] * x[0] + transform[6] * x[1] + transform[10] * x[2] + transform[14];
    return {y0/y3, y1/y3, y2/y3};
}

__device__
auto multmat(const float* a, const float* b) -> cuda::std::array<float,16>
{
    cuda::std::array<float,16> o;
    for(int i=0;i<4;++i) {
        for(int j=0;j<4;++j) {
           o[i + j * 4] = 0.;
            for(int k=0;k<4;++k) {
                o[i + j * 4] += a[i + k * 4] * b[k + j * 4];
            }
        }
    }
    return o;
}

__device__
auto jacobian_transform(const float* t, const float* p) -> cuda::std::array<float,9>
{
    const float a[9] = {
        t[0], t[1], t[2],
        t[4], t[5], t[6],
        t[8], t[9], t[10], };
    const float b[3] = {t[12], t[13], t[14]};
    const float d = t[15];
    const float c[3] = {t[3], t[7], t[11]};
    const float e = 1.f / (vec3::dot(c, p) + d);
    const float ee = e * e;
    const cuda::std::array<float,3> ap = mat3_col_major::mult_vec(a, p);
    const cuda::std::array<float,3> f = vec3::add(ap.data(), b);
    return {
        a[0] * e - f[0] * c[0] * ee,
        a[1] * e - f[1] * c[0] * ee,
        a[2] * e - f[2] * c[0] * ee,
        a[3] * e - f[0] * c[1] * ee,
        a[4] * e - f[1] * c[1] * ee,
        a[5] * e - f[2] * c[1] * ee,
        a[6] * e - f[0] * c[2] * ee,
        a[7] * e - f[1] * c[2] * ee,
        a[8] * e - f[2] * c[2] * ee,
    };
}

/// ray that goes through `pos_world: [f32;3]` that will be -z direction in the normalized device coordinate (NDC).
/// return `(ray_org: [f32;3], ray_dir: [f32;2])`
__device__
auto ray_from_transform_world2ndc(
    const float* transform_world2ndc,
    const float* pos_world,
    const float* transform_ndc2world
) -> cuda::std::pair<cuda::std::array<float,3>, cuda::std::array<float,3>> {
    const auto pos_mid_ndc = transform_homogeneous(transform_world2ndc, pos_world);
    const float pos_stt_ndc[3] = {pos_mid_ndc[0], pos_mid_ndc[1], 1.0};
    const float pos_end_ndc[3] = {pos_mid_ndc[0], pos_mid_ndc[1], -1.0};
    const auto ray_stt_world = transform_homogeneous(transform_ndc2world, pos_stt_ndc);
    const auto ray_end_world = transform_homogeneous(transform_ndc2world, pos_end_ndc);
    return cuda::std::make_pair(
        ray_stt_world,
        vec3::sub(ray_end_world.data(), ray_stt_world.data())
    );
}

}