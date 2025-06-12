#include <cuda/std/optional>
#include <cuda/std/tuple>
//
#include "tri2.h"

namespace edge2 {


struct RatioOfEdgeIntersection {
  float r0;
  float r1;
};

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
  const float* e1) -> cuda::std::optional< RatioOfEdgeIntersection >
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
    return RatioOfEdgeIntersection{r0, r1};
}

struct BackwardIntersectionEdge2 {
  cuda::std::array<float,2> dlds0;
  cuda::std::array<float,2> dlde0;
  cuda::std::array<float,2> dlds1;
  cuda::std::array<float,2> dlde1;
};

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
) -> BackwardIntersectionEdge2
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
    const auto s0e0s1_1 = tri2::dldw_area(s0, e0, s1, dlda1);
    const auto s0e0e1_2 = tri2::dldw_area(s0, e0, e1, dlda2);
    const auto s1e1s0_3 = tri2::dldw_area(s1, e1, s0, dlda3);
    const auto s1e1e0_4 = tri2::dldw_area(s1, e1, e0, dlda4);
    using cuda::std::get;
    return BackwardIntersectionEdge2 {
        vec2::add_three(get<0>(s0e0s1_1).data(), get<0>(s0e0e1_2).data(), get<2>(s1e1s0_3).data()),
        vec2::add_three(get<1>(s0e0s1_1).data(), get<1>(s0e0e1_2).data(), get<2>(s1e1e0_4).data()),
        vec2::add_three(get<2>(s0e0s1_1).data(), get<0>(s1e1s0_3).data(), get<0>(s1e1e0_4).data()),
        vec2::add_three(get<2>(s0e0e1_2).data(), get<1>(s1e1s0_3).data(), get<1>(s1e1e0_4).data())
    };
    /*
    const auto dlds0_1 = cuda::std::get<0>(d1);
    const auto dlde0_1 = cuda::std::get<1>(d1);
    const auto dlds1_1 = cuda::std::get<2>(d1);
    const auto dlds0_2 = cuda::std::get<0>(d2);
    const auto dlde0_2 = cuda::std::get<1>(d2);
    const auto dlde1_2 = cuda::std::get<2>(d2);
    const auto dlds1_3 = cuda::std::get<0>(d3);
    const auto dlde1_3 = cuda::std::get<1>(d3);
    const auto dlds0_3 = cuda::std::get<2>(d3);
    const auto dlds1_4 = cuda::std::get<0>(d4);
    const auto dlde1_4 = cuda::std::get<1>(d4);
    const auto dlde0_4 = cuda::std::get<2>(d4);
    return cuda::std::make_tuple(
        vec2::add_three(dlds0_1.data(), dlds0_2.data(), dlds0_3.data()),
        vec2::add_three(dlde0_1.data(), dlde0_2.data(), dlde0_4.data()),
        vec2::add_three(dlds1_1.data(), dlds1_3.data(), dlds1_4.data()),
        vec2::add_three(dlde1_2.data(), dlde1_3.data(), dlde1_4.data())
    );
    */
}

}