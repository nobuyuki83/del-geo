#pragma once

namespace pyramid {

// -------------------------------------------------------
// Pyramid (5 nodes) — Newton-Raphson
// N0=(1-r)(1-s)(1-t), N1=r(1-s)(1-t), N2=rs(1-t), N3=(1-r)s(1-t), N4=t
// domain: r,s,t in [0,1]
// -------------------------------------------------------
__device__ static bool parametric_coord_for_origin(
    float pco[3], const float p0[3], const float p1[3], const float p2[3],
    const float p3[3], const float p4[3]) {
  pco[0] = pco[1] = pco[2] = 0.5f;
  const float *pts[5] = {p0, p1, p2, p3, p4};
  for (int it = 0; it < 20; ++it) {
    float r = pco[0], s = pco[1], t = pco[2], rm = 1 - r, sm = 1 - s,
          tm = 1 - t;
    float sf[5] = {rm * sm * tm, r * sm * tm, r * s * tm, rm * s * tm, t};
    float pos[3] = {0, 0, 0};
    for (int n = 0; n < 5; ++n) {
      for (int d = 0; d < 3; ++d) {
        pos[d] += sf[n] * pts[n][d];
      }
    }
    float dndr[5][3] = {{-sm * tm, -rm * tm, -rm * sm},
                        {sm * tm, -r * tm, -r * sm},
                        {s * tm, r * tm, -r * s},
                        {-s * tm, rm * tm, -rm * s},
                        {0, 0, 1}};
    float dxdr[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int n = 0; n < 5; ++n) {
      for (int i_x = 0; i_x < 3; ++i_x) {
        for (int i_r = 0; i_r < 3; ++i_r) {
          dxdr[i_x][i_r] += pts[n][i_x] * dndr[n][i_r];
        }
      }
    }
    if (!mat3_array_of_rows::newton_step(pco, pos, dxdr)) return false;
  }
  static const float tol = 1e-5f;
  return pco[0] >= -tol && pco[0] <= 1 + tol && pco[1] >= -tol &&
         pco[1] <= 1 + tol && pco[2] >= -tol && pco[2] <= 1 + tol;
}

}  // namespace pyramid