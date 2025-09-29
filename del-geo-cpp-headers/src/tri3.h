#pragma once
#include <cuda/std/optional>
#include "vec3.h"
#include "mat3_col_major.h"

namespace tri3 {

/// normal vector of a 3D triangle (coordinates given by stack-allocated arrays)
__device__
auto normal(const float* v1, const  float* v2, const float* v3) -> cuda::std::array<float,3>
{
    return {
        (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]),
        (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]),
        (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]),
    };
}

struct DwNormal {
    cuda::std::array<float,9> d_p0;
    cuda::std::array<float,9> d_p1;
    cuda::std::array<float,9> d_p2;
};

__device__
auto dw_normal(const float* p0, const float* p1, const float* p2) -> DwNormal
{
    using V3f = cuda::std::array<float,3>;
    const auto v01 = vec3::sub(p1, p0);
    const auto v12 = vec3::sub(p2, p1);
    const auto v20 = vec3::sub(p0, p2);
    return DwNormal {
        mat3_col_major::from_vec3_to_skew_mat(v12.data()),
        mat3_col_major::from_vec3_to_skew_mat(v20.data()),
        mat3_col_major::from_vec3_to_skew_mat(v01.data()),
    };
}

__device__
auto intersection_against_ray(
    const float *p0,
    const float *p1,
    const float *p2,
    const float *ray_org,
    const float *ray_dir) -> cuda::std::optional<float>
{
    using V3f = cuda::std::array<float,3>;
    float eps = 1.0e-5;
    const V3f edge1 = vec3::sub(p1, p0);
    const V3f edge2 = vec3::sub(p2, p0);
    const V3f pvec = vec3::cross(ray_dir, edge2.data());
    const float det = vec3::dot(edge1.data(), pvec.data());
    if( det > -eps && det < eps ){
        return {};
    }
    float invdet = 1.f / det;
    const V3f tvec = vec3::sub(ray_org, p0);
    float u = invdet * vec3::dot(tvec.data(), pvec.data());
    if( u < 0.f || u > 1.f ){
        return {};
    }
    const V3f qvec = vec3::cross(tvec.data(), edge1.data());
    const float v = invdet * vec3::dot(ray_dir, qvec.data());
    if( v < 0.f || u + v > 1.f ){
        return {};
    }
    // At this stage we can compute t to find out where the intersection point is on the line.
    const float t = invdet * vec3::dot(edge2.data(), qvec.data());
    return t;
}

struct IntersectionAgainstLineBwdWrtTri {
    float t;
    float u;
    float v;
    cuda::std::array<float,3> d_p0;
    cuda::std::array<float,3> d_p1;
    cuda::std::array<float,3> d_p2;
};

/// ray triangle intersection.
/// * `dir` - any nonzero vector (not necessary to be a unit vector)
/// * `t` - ratio of `dir` vector from
/// * `u` - barycentric coordinate
/// * `v` - barycentric coordinate
///
/// `org + t * dir = (1 - u - v) * p0 + u * p1 + v * p2`
__device__
auto intersection_against_line_bwd_wrt_tri(
    const float* p0,
    const float* p1,
    const float* p2,
    const float* org,
    const float* dir,
    float d_t,
    float d_u,
    float d_v) -> cuda::std::optional<IntersectionAgainstLineBwdWrtTri>
{
    const auto e1 = vec3::sub(p0, p1);
    const auto e2 = vec3::sub(p2, p0);
    const auto n = vec3::cross(e1.data(), e2.data());
    const float n_dot_dir = vec3::dot(n.data(), dir);
    if( n_dot_dir == 0.f ){
        return {};
    }
    const float inv_det = 1.f / n_dot_dir;
    const auto c = vec3::sub(p0, org);
    const float n_dot_c = vec3::dot(n.data(), c.data());
    const float t = inv_det * n_dot_c;
    if( t < 0.f ){
        return {};
    }
    const auto r = vec3::cross(dir, c.data());
    //
    const float r_dot_e2 = vec3::dot(r.data(), e2.data());
    const float u = inv_det * r_dot_e2;
    if( u < 0.f ){
        return {};
    }
    //
    const float r_dot_e1 = vec3::dot(r.data(), e1.data());
    const float v = inv_det * r_dot_e1;
    if( v < 0.f ){
        return {};
    }
    if( u + v >= 1.f ){
        return {};
    }
    // --------------
    // below: bwd
    const float d_n_dot_c = d_t * inv_det;
    const float d_r_dot_e2 = d_u * inv_det;
    const float d_r_dot_e1 = d_v * inv_det;
    const float d_inv_det = d_t * n_dot_c + d_u * r_dot_e2 + d_v * r_dot_e1;
    //
    auto d_n = vec3::scale(c.data(), d_n_dot_c);
    auto d_c = vec3::scale(n.data(), d_n_dot_c);
    auto d_e2 = vec3::scale(r.data(), d_r_dot_e2);
    auto d_e1 = vec3::scale(r.data(), d_r_dot_e1);
    auto d_r = vec3::add( vec3::scale(e2.data(), d_r_dot_e2).data(), vec3::scale(e1.data(), d_r_dot_e1).data() );
    //
    const float d_n_dot_dir = -d_inv_det / n_dot_dir / n_dot_dir;
    vec3::add_inplace(d_n.data(), vec3::scale(dir, d_n_dot_dir).data());
    vec3::add_inplace(d_c.data(), vec3::cross(d_r.data(), dir).data());
    vec3::add_inplace(d_e2.data(), vec3::cross(d_n.data(), e1.data()).data());
    vec3::add_inplace(d_e1.data(), vec3::cross(e2.data(), d_n.data()).data());
    //
    const auto d_p0 = vec3::add( vec3::sub(d_e1.data(), d_e2.data()).data(), d_c.data());
    const auto d_p1 = vec3::scale(d_e1.data(), -1.f);
    const auto d_p2 = d_e2;
    return IntersectionAgainstLineBwdWrtTri {t, u, v, d_p0, d_p1, d_p2};
}

}