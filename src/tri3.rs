use num_traits::AsPrimitive;

pub fn area_<T>(p0: &[T], p1: &[T], p2: &[T]) -> T
    where T: num_traits::Float + 'static,
          f64: num_traits::AsPrimitive<T>
{
    use crate::vec3;
    assert!(p0.len() == 3 && p1.len() == 3 && p2.len() == 3);
    let v1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let v2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let na = [
        v1[1] * v2[2] - v2[1] * v1[2],
        v1[2] * v2[0] - v2[2] * v1[0],
        v1[0] * v2[1] - v2[0] * v1[1]];
    return vec3::squared_norm_(&na).sqrt() * 0.5_f64.as_();
}

pub fn normal_<T>(
    vnorm: &mut [T],
    v1: &[T],
    v2: &[T],
    v3: &[T])
    where T: std::ops::Sub<Output=T> + std::ops::Mul<Output=T> + std::ops::Sub + Copy
{
    vnorm[0] = (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]);
    vnorm[1] = (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]);
    vnorm[2] = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]);
}

pub fn unit_normal_<T>(
    n: &mut [T],
    v1: &[T],
    v2: &[T],
    v3: &[T]) -> T
    where T: std::ops::Sub<Output=T> + std::ops::Mul<Output=T> + num_traits::Float + 'static + Copy + std::ops::MulAssign,
          f32: num_traits::AsPrimitive<T>
{
    use crate::vec3;
    normal_(
        n,
        v1, v2, v3);
    let a = vec3::norm_(n) * 0.5_f32.as_();
    let invlen: T = 0.5_f32.as_() / a;
    n[0] *= invlen;
    n[1] *= invlen;
    n[2] *= invlen;
    a
}

pub fn area_and_unorm_<T>(
    v1: &[T],
    v2: &[T],
    v3: &[T]) -> (T, [T; 3])
    where T: std::ops::Sub<Output=T> + std::ops::Mul<Output=T> + num_traits::Float + 'static + Copy + std::ops::MulAssign,
          f32: num_traits::AsPrimitive<T>
{
    use crate::vec3;
    let mut n: [T; 3] = [0_f32.as_(); 3];
    normal_(
        &mut n,
        v1, v2, v3);
    let a = vec3::norm_(&n) * 0.5_f32.as_();
    let invlen: T = 0.5_f32.as_() / a;
    n[0] *= invlen;
    n[1] *= invlen;
    n[2] *= invlen;
    (a, n)
}

pub fn cot_<T>(
    p0: &[T],
    p1: &[T],
    p2: &[T]) -> [T; 3]
    where T: num_traits::Float + 'static,
          f32: num_traits::AsPrimitive<T>
{
    use crate::vec3;
    assert!(p0.len() == 3 && p1.len() == 3 && p2.len() == 3);
    let v0 = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
    let v1 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let v2 = [p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2]];
    let na = [
        v1[1] * v2[2] - v2[1] * v1[2],
        v1[2] * v2[0] - v2[2] * v1[0],
        v1[0] * v2[1] - v2[0] * v1[1]];
    let area: T = vec3::squared_norm_(&na).sqrt() * 0.5_f32.as_();
    let tmp: T = 0.25_f32.as_() / area;
    let l0 = vec3::squared_norm_(&v0);
    let l1 = vec3::squared_norm_(&v1);
    let l2 = vec3::squared_norm_(&v2);
    [
        (l1 + l2 - l0) * tmp,
        (l2 + l0 - l1) * tmp,
        (l0 + l1 - l2) * tmp
    ]
}

// ----------

/// Möller–Trumbore ray-triangle intersection algorithm
/// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
/// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
pub fn ray_triangle_intersection_(
    ray_org: &[f32],
    ray_dir: &[f32],
    p0: &[f32],
    p1: &[f32],
    p2: &[f32]) -> Option<f32> {
    use crate::vec3;
    // const EPSILON: f32 = 1.0e-4;
    let edge1 = vec3::sub_(p1, p0);
    let edge2 = vec3::sub_(p2, p0);
    let pvec = vec3::cross_(&ray_dir, &edge2);
    let det = vec3::dot_(&edge1, &pvec);
    // if det > -EPSILON && det < EPSILON { return None; }
    let invdet = 1.0 / det;
    let tvec = vec3::sub_(ray_org, p0);
    let u = invdet * vec3::dot_(&tvec, &pvec);
    if u < 0.0 || u > 1.0 { return None; }
    let qvec = vec3::cross_(&tvec, &edge1);
    let v = invdet * vec3::dot_(ray_dir, &qvec);
    if v < 0.0 || u + v > 1.0 { return None; }
    // At this stage we can compute t to find out where the intersection point is on the line.
    let t = invdet * vec3::dot_(&edge2, &qvec);
    return Some(t);
}

pub fn nearest_to_point3_(
    ps: &[f32],
    q0: &[f32],
    q1: &[f32],
    q2: &[f32]) -> ([f32; 3], f32, f32) {
    use crate::{tet, edge3, vec3};
    let (_area, n012) = area_and_unorm_(q0, q1, q2);
    let pe = [ps[0] + n012[0], ps[1] + n012[1], ps[2] + n012[2]];
    let v012 = tet::volume_(ps, q0, q1, q2);
    if v012.abs() > 1.0e-10 {
        let sign: f32 = if v012 > 0_f32 { 1_f32 } else { -1_f32 };
        let v0: f32 = tet::volume_(ps, q1, q2, &pe) * sign;
        let v1: f32 = tet::volume_(ps, q2, q0, &pe) * sign;
        let v2: f32 = tet::volume_(ps, q0, q1, &pe) * sign;
        assert!((v0 + v1 + v2).abs() > 1.0e-10);
        let inv_v012 = 1.0 / (v0 + v1 + v2);
        let r0 = v0 * inv_v012;
        let r1 = v1 * inv_v012;
        let r2 = 1.0 - r0 - r1;
        let tol = 1.0e-4;
        if r0 > -tol && r1 > -tol && r2 > -tol {
            let nearp = [
                q0[0] * r0 + q1[0] * r1 + q2[0] * r2,
                q0[1] * r0 + q1[1] * r1 + q2[1] * r2,
                q0[2] * r0 + q1[2] * r1 + q2[2] * r2];
            return (nearp, r0, r1);
        }
    }
    let r12: [f32; 3] = edge3::nearest_point3_(ps, q1, q2);
    let r20: [f32; 3] = edge3::nearest_point3_(ps, q2, q0);
    let r01: [f32; 3] = edge3::nearest_point3_(ps, q0, q1);
    let d12 = vec3::distance_(&r12, ps);
    let d20 = vec3::distance_(&r20, ps);
    let d01 = vec3::distance_(&r01, ps);
    if d12 < d20 {
        if d12 < d01 { // 12 is the smallest
            let nearp = [r12[0], r12[1], r12[2]];
            let r0 = 0_f32;
            let r1 = vec3::distance_(&nearp, q2) / vec3::distance_(q1, q2);
            return (nearp, r0, r1);
        }
    } else {
        if d20 < d01 { // d20 is the smallest
            let nearp = [r20[0], r20[1], r20[2]];
            let r0 = vec3::distance_(&nearp, q2) / vec3::distance_(q0, q2);
            let r1 = 0_f32;
            return (nearp, r0, r1);
        }
    }
    let nearp = [r01[0], r01[1], r01[2]];
    let r0 = vec3::distance_(&nearp, q1) / vec3::distance_(q0, q1);
    let r1 = 1_f32 - r0;
    return (nearp, r0, r1);
}

/* ------------------------------------ */

pub fn height<T>(
    q: &nalgebra::Vector3::<T>,
    p0: &nalgebra::Vector3::<T>,
    p1: &nalgebra::Vector3::<T>) -> T
    where T: nalgebra::RealField + 'static + Copy + num_traits::Float,
        f64: AsPrimitive<T>
{
    let a = area_(q.as_slice(), p0.as_slice(), p1.as_slice());
    a * 2.0.as_() / (p0 - p1).norm()
}