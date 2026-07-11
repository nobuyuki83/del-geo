#pragma once
#include "mat3_array_of_rows.h"

namespace prism {

// -------------------------------------------------------
// Prism (6 nodes) — Newton-Raphson
// N0=(1-r-s)(1-t), N1=r(1-t), N2=s(1-t), N3=(1-r-s)t, N4=rt, N5=st
// domain: r>=0, s>=0, r+s<=1, t in [0,1]
// -------------------------------------------------------
__device__ static bool parametric_coord_for_origin(
    float pco[3], const float p0[3], const float p1[3], const float p2[3],
    const float p3[3], const float p4[3], const float p5[3]) {
  pco[0] = 1.f / 3;
  pco[1] = 1.f / 3;
  pco[2] = 0.5f;
  const float *pts[6] = {p0, p1, p2, p3, p4, p5};
  for (int it = 0; it < 20; ++it) {
    float r = pco[0], s = pco[1], t = pco[2], rs = 1 - r - s, tm = 1 - t;
    float sf[6] = {rs * tm, r * tm, s * tm, rs * t, r * t, s * t};
    float pos[3] = {0, 0, 0};
    for (int n = 0; n < 6; ++n) {
      for (int d = 0; d < 3; ++d) {
        pos[d] += sf[n] * pts[n][d];
      }
    }
    float dndr[6][3] = {{-tm, -tm, -rs}, {tm, 0, -r}, {0, tm, -s},
                        {-t, -t, rs},    {t, 0, r},   {0, t, s}};
    float dxdr[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int n = 0; n < 6; ++n) {
      for (int i_x = 0; i_x < 3; ++i_x) {
        for (int i_r = 0; i_r < 3; ++i_r) {
          dxdr[i_x][i_r] += pts[n][i_x] * dndr[n][i_r];
        }
      }
    }
    if (!mat3_array_of_rows::newton_step(pco, pos, dxdr)) return false;
  }
  static const float tol = 1e-5f;
  float r = pco[0], s = pco[1], t = pco[2];
  return r >= -tol && s >= -tol && r + s <= 1 + tol && t >= -tol &&
         t <= 1 + tol;
}

}  // namespace prism