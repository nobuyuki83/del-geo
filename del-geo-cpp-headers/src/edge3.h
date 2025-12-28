#pragma once

#include <cuda/std/optional>
#include <cuda/std/tuple>
//

namespace edge3 {

__device__
float length(
  const float* ps,
  const float* pe) {
  float x = pe[0] - ps[0];
  float y = pe[1] - ps[1];
  float z = pe[2] - ps[2];
  return sqrt(x*x + y*y + z*z);
}

}