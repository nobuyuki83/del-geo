#pragma once
#include "mat3_array_of_rows.h"

namespace hex {

// -------------------------------------------------------
// Hex (8 nodes) — Newton-Raphson
// N_i = (1/8)(1±r)(1±s)(1±t), domain r,s,t in [-1,1]
// -------------------------------------------------------
__device__ static bool parametric_coord_for_origin(
    float pco[3], const float p0[3], const float p1[3], const float p2[3],
    const float p3[3], const float p4[3], const float p5[3], const float p6[3],
    const float p7[3]) {
  pco[0] = pco[1] = pco[2] = 0.f;
  const float *pts[8] = {p0, p1, p2, p3, p4, p5, p6, p7};
  for (int it = 0; it < 20; ++it) {
    float r = pco[0], s = pco[1], t = pco[2], e = 0.125f;
    float sf[8] = {
        e * (1 - r) * (1 - s) * (1 - t), e * (1 + r) * (1 - s) * (1 - t),
        e * (1 + r) * (1 + s) * (1 - t), e * (1 - r) * (1 + s) * (1 - t),
        e * (1 - r) * (1 - s) * (1 + t), e * (1 + r) * (1 - s) * (1 + t),
        e * (1 + r) * (1 + s) * (1 + t), e * (1 - r) * (1 + s) * (1 + t)};
    float pos[3] = {0, 0, 0};
    for (int n = 0; n < 8; ++n) {
      for (int d = 0; d < 3; ++d) {
        pos[d] += sf[n] * pts[n][d];
      }
    }
    float dndr[8][3] = {
        {-e * (1 - s) * (1 - t), -e * (1 - r) * (1 - t),
         -e * (1 - r) * (1 - s)},
        {e * (1 - s) * (1 - t), -e * (1 + r) * (1 - t), -e * (1 + r) * (1 - s)},
        {e * (1 + s) * (1 - t), e * (1 + r) * (1 - t), -e * (1 + r) * (1 + s)},
        {-e * (1 + s) * (1 - t), e * (1 - r) * (1 - t), -e * (1 - r) * (1 + s)},
        {-e * (1 - s) * (1 + t), -e * (1 - r) * (1 + t), e * (1 - r) * (1 - s)},
        {e * (1 - s) * (1 + t), -e * (1 + r) * (1 + t), e * (1 + r) * (1 - s)},
        {e * (1 + s) * (1 + t), e * (1 + r) * (1 + t), e * (1 + r) * (1 + s)},
        {-e * (1 + s) * (1 + t), e * (1 - r) * (1 + t), e * (1 - r) * (1 + s)}};
    float dxdr[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int n = 0; n < 8; ++n) {
      for (int i_x = 0; i_x < 3; ++i_x) {
        for (int i_r = 0; i_r < 3; ++i_r) {
          dxdr[i_x][i_r] += pts[n][i_x] * dndr[n][i_r];
        }
      }
    }
    if (!mat3_array_of_rows::newton_step(pco, pos, dxdr)) return false;
  }
  static const float tol = 1e-5f;
  return pco[0] >= -1 - tol && pco[0] <= 1 + tol && pco[1] >= -1 - tol &&
         pco[1] <= 1 + tol && pco[2] >= -1 - tol && pco[2] <= 1 + tol;
}

}  // namespace hex
