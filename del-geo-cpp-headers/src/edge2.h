#include <cuda/std/optional>
#include <cuda/std/tuple>
//
#include "tri2.h"

namespace edge2 {

///
/// * Return
///
/// Some((r0,r1)): ratio of edge
///
/// r0 * s0 + (1-r0) * e0 == r1 * s1 + (1-r1) * e1
__device__
auto intersection_edge2(
  const float* s0,
  const float* e0,
  const float* s1,
  const float* e1) -> cuda::std::optional< cuda::std::pair<float,float> >
{
    const float area1 = tri2::area(s0, e0, s1);
    const float area2 = tri2::area(s0, e0, e1);
    const float area3 = tri2::area(s1, e1, s0);
    const float area4 = tri2::area(s1, e1, e0);
    if( area1 * area2 > 0.f ){
        return cuda::std::nullopt;
    }
    if( area3 * area4 > 0.f ){
        return cuda::std::nullopt;
    }
    const float r1 = area1 / (area1 - area2);
    const float r0 = area3 / (area3 - area4);
    return cuda::std::make_pair(r0, r1);
}

/// # return
///
///  `(dls0: [T;2], dle0: [T;2], dlds1: [T;2], dlde1: [T;2])`
__device__
auto dldw_intersection_edge2(
    const float* s0,
    const float* e0,
    const float* s1,
    const float* e1,
    float dldr0,
    float dldr1
) -> cuda::std::tuple<
  cuda::std::array<float,2>,
  cuda::std::array<float,2>,
  cuda::std::array<float,2>,
  cuda::std::array<float,2>>
{
    const float a1 = tri2::area(s0, e0, s1);
    const float a2 = tri2::area(s0, e0, e1);
    const float a3 = tri2::area(s1, e1, s0);
    const float a4 = tri2::area(s1, e1, e0);
    //let r1 = a1 / (a1 - a2);
    //let r0 = a3 / (a3 - a4);
    const float dlda1 = dldr1 * (1.f / (a1 - a2) - a1 / pow(a1 - a2, 2));
    const float dlda2 = dldr1 * (a1 / pow(a1 - a2, 2));
    const float dlda3 = dldr0 * (1.f / (a3 - a4) - a3 / pow(a3 - a4, 2));
    const float dlda4 = dldr0 * (a3 / pow(a3 - a4, 2));
    const auto [dlds0_1, dlde0_1, dlds1_1] = tri2::dldw_area(s0, e0, s1, dlda1);
    const auto [dlds0_2, dlde0_2, dlde1_2] = tri2::dldw_area(s0, e0, e1, dlda2);
    const auto [dlds1_3, dlde1_3, dlds0_3] = tri2::dldw_area(s1, e1, s0, dlda3);
    const auto [dlds1_4, dlde1_4, dlde0_4] = tri2::dldw_area(s1, e1, e0, dlda4);
    return cuda::std::make_tuple(
        vec2::add_three(dlds0_1.data(), dlds0_2.data(), dlds0_3.data()),
        vec2::add_three(dlde0_1.data(), dlde0_2.data(), dlde0_4.data()),
        vec2::add_three(dlds1_1.data(), dlds1_3.data(), dlds1_4.data()),
        vec2::add_three(dlde1_2.data(), dlde1_3.data(), dlde1_4.data())
    );
}

}