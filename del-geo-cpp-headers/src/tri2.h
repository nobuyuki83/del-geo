#pragma once
#include <cuda/std/tuple>
//
#include "vec2.h"

namespace tri2 {

__device__
float area(
  const float* p0,
  const float* p1,
  const float* p2)
{
  return 0.5f * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
}


/// # Return
/// `(dldp0: [T;2], dldp1: [T;2], dldp2: [T;2])`
__device__
auto dldw_area(
  const float* p0,
  const float* p1,
  const float* p2,
  float dldarea) -> cuda::std::tuple<
  cuda::std::array<float,2>,
  cuda::std::array<float,2>,
  cuda::std::array<float,2>>
{
    const float dareadp0x2[3] = {p1[1] - p2[1], p2[0] - p1[0]};
    const float dareadp1x2[3] = {p2[1] - p0[1], p0[0] - p2[0]};
    const float dareadp2x2[3] = {p0[1] - p1[1], p1[0] - p0[0]};
    return cuda::std::make_tuple(
        vec2::scale(dareadp0x2, dldarea * 0.5),
        vec2::scale(dareadp1x2, dldarea * 0.5),
        vec2::scale(dareadp2x2, dldarea * 0.5)
    );
}


}