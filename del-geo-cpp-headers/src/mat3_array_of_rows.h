#pragma once
#include <math.h>

namespace mat3_array_of_rows {

// -------------------------------------------------------
// 3x3 column-major matrix utilities
// -------------------------------------------------------

__device__ static float mat3_det(const float m[3][3]) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[1][0] * (m[0][1] * m[2][2] - m[0][2] * m[2][1]) +
         m[2][0] * (m[0][1] * m[1][2] - m[0][2] * m[1][1]);
}

__device__ static bool mat3_inv(float inv[3][3], const float m[3][3]) {
  float d = mat3_det(m);
  if (fabsf(d) < 1e-30f) return false;
  float id = 1.f / d;
  inv[0][0] = id * (m[1][1] * m[2][2] - m[1][2] * m[2][1]);
  inv[0][1] = -id * (m[0][1] * m[2][2] - m[0][2] * m[2][1]);
  inv[0][2] = id * (m[0][1] * m[1][2] - m[0][2] * m[1][1]);
  inv[1][0] = -id * (m[1][0] * m[2][2] - m[1][2] * m[2][0]);
  inv[1][1] = id * (m[0][0] * m[2][2] - m[0][2] * m[2][0]);
  inv[1][2] = -id * (m[0][0] * m[1][2] - m[0][2] * m[1][0]);
  inv[2][0] = id * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  inv[2][1] = -id * (m[0][0] * m[2][1] - m[0][1] * m[2][0]);
  inv[2][2] = id * (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
  return true;
}

__device__ static void mat3_mulvec(float out[3], const float m[3][3],
                                   const float v[3]) {
  for (int i = 0; i < 3; ++i) {
    out[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
  }
}

// r -= (dx/dr)^{-1} * dx  where dxdr is stored as dxdr[xyz][param]
__device__ static bool newton_step(
  float r[3],
  const float dx[3],
  const float dxdr[3][3])
{
  float drdx[3][3];
  if (!mat3_inv(drdx, dxdr)) return false;
  float dr[3];
  mat3_mulvec(dr, drdx, dx);
  for (int i = 0; i < 3; ++i) {
    r[i] -= dr[i];
  }
  return true;
}

}  // namespace mat3_array_of_rows