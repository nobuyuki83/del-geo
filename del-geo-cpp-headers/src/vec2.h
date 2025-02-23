#pragma once
#include <cuda/std/array>

namespace vec2 {

__device__
auto scale(const float* x, float a) -> cuda::std::array<float,2>
{
    return {
        a * x[0],
        a * x[1],
    };
}

__device__
auto add_three(
  const float* a,
  const float* b,
  const float* c) -> cuda::std::array<float,2>
{
    return {
      a[0] + b[0] + c[0],
      a[1] + b[1] + c[1]};
}

}