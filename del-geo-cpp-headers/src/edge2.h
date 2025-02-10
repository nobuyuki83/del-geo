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

}