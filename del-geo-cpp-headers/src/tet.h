#pragma once

namespace tet {

// -------------------------------------------------------
// Tet (4 nodes): barycentric coordinates
// p0..p3 already shifted: pi = vi - query
// -------------------------------------------------------
__device__ static float tet_vol6(const float a[3], const float b[3],
                                 const float c[3], const float d[3]) {
  float ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
  float ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
  float ad[3] = {d[0] - a[0], d[1] - a[1], d[2] - a[2]};
  return ab[0] * (ac[1] * ad[2] - ac[2] * ad[1]) -
         ab[1] * (ac[0] * ad[2] - ac[2] * ad[0]) +
         ab[2] * (ac[0] * ad[1] - ac[1] * ad[0]);
}

__device__ static bool tet_parametric(float bc[3], const float p0[3],
                                      const float p1[3], const float p2[3],
                                      const float p3[3]) {
  const float z[3] = {0, 0, 0};
  float total = tet_vol6(p0, p1, p2, p3);
  if (fabsf(total) < 1e-30f) return false;
  float r0 = tet_vol6(z, p1, p2, p3) / total;
  float r1 = tet_vol6(p0, z, p2, p3) / total;
  float r2 = tet_vol6(p0, p1, z, p3) / total;
  float r3 = 1.f - r0 - r1 - r2;
  if (r0 < 0 || r1 < 0 || r2 < 0 || r3 < 0) return false;
  bc[0] = r0;
  bc[1] = r1;
  bc[2] = r2;
  return true;
}

}  // namespace tet