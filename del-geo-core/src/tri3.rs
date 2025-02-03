//! methods for 3d triangle

use num_traits::AsPrimitive;

/// clamp barycentric coordinates inside a triangle
pub fn clamp<T>(r0: T, r1: T, r2: T) -> (T, T, T)
where
    T: num_traits::Float,
{
    // vertex projection
    if r0 <= T::zero() && r1 <= T::zero() {
        return (T::zero(), T::zero(), T::one());
    }
    if r1 <= T::zero() && r2 <= T::zero() {
        return (T::one(), T::zero(), T::zero());
    }
    if r2 <= T::zero() && r0 <= T::zero() {
        return (T::zero(), T::one(), T::zero());
    }
    // edge projection
    if r0 <= T::zero() {
        return (T::zero(), r1 / (r1 + r2), r2 / (r1 + r2));
    }
    if r1 <= T::zero() {
        return (r0 / (r0 + r2), T::zero(), r2 / (r0 + r2));
    }
    if r2 <= T::zero() {
        return (r0 / (r0 + r1), r1 / (r0 + r1), T::zero());
    }
    (r0, r1, r2)
}

// ----------------------

/// normal vector of a 3D triangle (coordinates given by stack-allocated arrays)
pub fn normal<T>(v1: &[T; 3], v2: &[T; 3], v3: &[T; 3]) -> [T; 3]
where
    T: std::ops::Sub<Output = T> + std::ops::Mul<Output = T> + std::ops::Sub + Copy,
{
    [
        (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v2[2] - v1[2]) * (v3[1] - v1[1]),
        (v2[2] - v1[2]) * (v3[0] - v1[0]) - (v2[0] - v1[0]) * (v3[2] - v1[2]),
        (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0]),
    ]
}

/// area of a 3D triangle (coordinates given by stack-allocated arrays)
pub fn area<T>(v1: &[T; 3], v2: &[T; 3], v3: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let na = normal(v1, v2, v3);
    let half = T::one() / (T::one() + T::one());
    na.squared_norm().sqrt() * half
}

/// height of triangle vertex `p2` against the edge connecting `p0` and `p1`
pub fn height<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3]) -> T
where
    T: num_traits::Float + 'static + Copy + num_traits::Float,
    f64: AsPrimitive<T>,
{
    use crate::vec3::Vec3;
    let a = area(p2, p0, p1);
    a * 2.0.as_() / p0.sub(p1).norm()
}

pub fn unit_normal_area<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3]) -> ([T; 3], T)
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let n = normal(p0, p1, p2);
    let half = T::one() / (T::one() + T::one());
    let a = n.norm() * half;
    let invlen: T = half / a;
    ([n[0] * invlen, n[1] * invlen, n[2] * invlen], a)
}

/// compute cotangents of the three angles of a triangle
pub fn cot<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    assert!(p0.len() == 3 && p1.len() == 3 && p2.len() == 3);
    let v0 = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
    let v1 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let v2 = [p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2]];
    let half = T::one() / (T::one() + T::one());
    let onefourth = T::one() / (T::one() + T::one() + T::one() + T::one());
    let area = {
        let na = [
            v1[1] * v2[2] - v2[1] * v1[2],
            v1[2] * v2[0] - v2[2] * v1[0],
            v1[0] * v2[1] - v2[0] * v1[1],
        ];
        na.squared_norm().sqrt() * half
    };
    let tmp: T = onefourth / area;
    let l0 = v0.squared_norm();
    let l1 = v1.squared_norm();
    let l2 = v2.squared_norm();
    [
        (l1 + l2 - l0) * tmp, // cot0 = cos0/sin0 = {(l1*l1+l2*l2-l0*l0)/(2*l1*l2)} / {2*area/(l1*l2)}
        (l2 + l0 - l1) * tmp, // cot1 = cos1/sin1 = {(l2*l2+l0*l0-l1*l1)/(2*l2*l0)} / {2*area/(l2*l0)}
        (l0 + l1 - l2) * tmp, // cot2 = cos2/sin2 = {(l0*l0+l1*l1-l2*l2)/(2*l0*l1)} / {2*area/(l0*l1)}
    ]
}

pub fn emat_cotangent_laplacian<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3]) -> [[[T; 1]; 3]; 3]
where
    T: num_traits::Float,
{
    let cots = cot(p0, p1, p2);
    let half = T::one() / (T::one() + T::one());
    [
        [
            [(cots[1] + cots[2]) * half],
            [-cots[2] * half],
            [-cots[1] * half],
        ],
        [
            [-cots[2] * half],
            [(cots[2] + cots[0]) * half],
            [-cots[0] * half],
        ],
        [
            [-cots[1] * half],
            [-cots[0] * half],
            [(cots[0] + cots[1]) * half],
        ],
    ]
}

pub fn emat_graph_laplacian<T>(l: T) -> [[[T; 1]; 3]; 3]
where
    T: num_traits::Float,
{
    let vo = -T::one() * l;
    let vd = (T::one() + T::one()) * l;
    [[[vd], [vo], [vo]], [[vo], [vd], [vo]], [[vo], [vo], [vd]]]
}

/// Möller–Trumbore ray-triangle intersection algorithm
///
/// <https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm>
/// <https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection>
pub fn intersection_against_line<T>(
    p0: &[T; 3],
    p1: &[T; 3],
    p2: &[T; 3],
    ray_org: &[T; 3],
    ray_dir: &[T; 3],
) -> Option<T>
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let eps: T = T::epsilon();
    let edge1 = p1.sub(p0);
    let edge2 = p2.sub(p0);
    let pvec = ray_dir.cross(&edge2);
    let det = edge1.dot(&pvec);
    if det > -eps && det < eps {
        return None;
    }
    let invdet = T::one() / det;
    let tvec = ray_org.sub(p0);
    let u = invdet * tvec.dot(&pvec);
    if u < T::zero() || u > T::one() {
        return None;
    }
    let qvec = tvec.cross(&edge1);
    let v = invdet * ray_dir.dot(&qvec);
    if v < T::zero() || u + v > T::one() {
        return None;
    }
    // At this stage we can compute t to find out where the intersection point is on the line.
    let t = invdet * edge2.dot(&qvec);
    Some(t)
}

/// ray triangle intersection.
/// * `dir` - any nonzero vector (not necessary to be a unit vector)
/// * `t` - ratio of `dir` vector from
/// * `u` - barycentric coordinate
/// * `v` - barycentric coordinate
///
/// `org + t * dir = (1 - u - v) * p0 + u * p1 + v * p2`
#[allow(clippy::too_many_arguments)]
#[allow(clippy::type_complexity)]
pub fn intersection_against_line_bwd_wrt_tri<Real>(
    p0: &[Real; 3],
    p1: &[Real; 3],
    p2: &[Real; 3],
    org: &[Real; 3],
    dir: &[Real; 3],
    d_t: Real,
    d_u: Real,
    d_v: Real,
) -> (Real, Real, Real, [Real; 3], [Real; 3], [Real; 3])
where
    Real: num_traits::Float + Copy,
{
    use crate::vec3::Vec3;
    let e1 = p0.sub(p1);
    let e2 = p2.sub(p0);
    let n = e1.cross(&e2);
    let n_dot_dir = n.dot(dir);
    /*
    if n_dot_dir.is_zero() {
        return None;
    }
     */
    let inv_det = Real::one() / n_dot_dir;
    let c = p0.sub(org);
    let n_dot_c = n.dot(&c);
    let t = inv_det * n_dot_c;
    /*
    if t < Real::zero() {
        return None;
    }
     */
    let r = dir.cross(&c);
    //
    let r_dot_e2 = r.dot(&e2);
    let u = inv_det * r_dot_e2;
    /*
    if u < Real::zero() {
        return None;
    }
     */
    //
    let r_dot_e1 = r.dot(&e1);
    let v = inv_det * r_dot_e1;
    /*
    if v < Real::zero() {
        return None;
    }
    if u + v >= Real::one() {
        return None;
    }
     */
    // --------------
    // below: bwd
    let d_n_dot_c = d_t * inv_det;
    let d_r_dot_e2 = d_u * inv_det;
    let d_r_dot_e1 = d_v * inv_det;
    let d_inv_det = d_t * n_dot_c + d_u * r_dot_e2 + d_v * r_dot_e1;
    //
    let mut d_n = c.scale(d_n_dot_c);
    let mut d_c = n.scale(d_n_dot_c);
    let mut d_e2 = r.scale(d_r_dot_e2);
    let mut d_e1 = r.scale(d_r_dot_e1);
    let d_r = e2.scale(d_r_dot_e2).add(&e1.scale(d_r_dot_e1));
    //
    let d_n_dot_dir = -d_inv_det / n_dot_dir / n_dot_dir;
    d_n.add_in_place(&dir.scale(d_n_dot_dir));
    d_c.add_in_place(&d_r.cross(dir));
    d_e2.add_in_place(&d_n.cross(&e1));
    d_e1.add_in_place(&e2.cross(&d_n));
    //
    let d_p0 = d_e1.sub(&d_e2).add(&d_c);
    let d_p1 = d_e1.scale(-Real::one());
    let d_p2 = d_e2;
    (t, u, v, d_p0, d_p1, d_p2)
}

#[test]
fn test_dw_ray_triangle_intersection() {
    use crate::vec3::Vec3;
    type Real = f64;
    let p0: [[Real; 3]; 3] = [[-13., -5., 8.], [14., -5., 8.], [1., 3., -3.]];

    // let origin = nalgebra::Vector3::<Real>::new(9., 11., 12.);
    let origin = [8., 11., 10.];
    let dir = [1., 0., 2.].sub(&origin);
    // d_t, d_u, d_u are back-propagated from the loss
    let d_t = 0.7;
    let d_u = 1.3;
    let d_v = 1.1;

    let (t0, u0, v0, d_p0, d_p1, d_p2) =
        intersection_against_line_bwd_wrt_tri(&p0[0], &p0[1], &p0[2], &origin, &dir, d_t, d_u, d_v);
    {
        let t1 = intersection_against_line(&p0[0], &p0[1], &p0[2], &origin, &dir).unwrap();
        assert!((t0 - t1).abs() < 1.0e-5);
    }

    let ha = p0[0]
        .scale(1. - u0 - v0)
        .add(&p0[1].scale(u0))
        .add(&p0[2].scale(v0));
    assert!(crate::tet::volume(&p0[0], &p0[1], &p0[2], &ha).abs() < 1.0e-8);
    let hb = crate::vec3::axpy(t0, &dir, &origin);
    assert!(ha.sub(&hb).norm() < 1.0e-6);

    let dp = [d_p0, d_p1, d_p2];

    let eps = 1.0e-5;
    for (i_node, i_dim) in itertools::iproduct!(0..3, 0..3) {
        let p1 = {
            let mut p1 = p0;
            p1[i_node][i_dim] += eps;
            p1
        };
        let (t1, u1, v1, _d_p0, _d_p1, _d_p2) = intersection_against_line_bwd_wrt_tri(
            &p1[0], &p1[1], &p1[2], &origin, &dir, d_t, d_u, d_v,
        );
        let dloss = ((t1 - t0) * d_t + (u1 - u0) * d_u + (v1 - v0) * d_v) / eps;
        let dloss_analytic = dp[i_node][i_dim];
        // dbg!(dloss, dloss_analytic);
        assert!(
            (dloss - dloss_analytic).abs() < 6.0e-7,
            "{} {} {}",
            dloss_analytic,
            dloss,
            dloss - dloss_analytic
        );
    }
}

// ----------------------------

pub fn to_barycentric_coords<T>(p0: &[T; 3], p1: &[T; 3], p2: &[T; 3], q: &[T; 3]) -> [T; 3]
where
    T: num_traits::Float,
{
    let a0 = area(q, p1, p2);
    let a1 = area(q, p2, p0);
    let a2 = area(q, p0, p1);
    let sum_inv = T::one() / (a0 + a1 + a2);
    [a0 * sum_inv, a1 * sum_inv, a2 * sum_inv]
}

pub fn position_from_barycentric_coords<T>(
    p0: &[T; 3],
    p1: &[T; 3],
    p2: &[T; 3],
    bc: &[T; 3],
) -> [T; 3]
where
    T: num_traits::Float,
{
    [
        bc[0] * p0[0] + bc[1] * p1[0] + bc[2] * p2[0],
        bc[0] * p0[1] + bc[1] * p1[1] + bc[2] * p2[1],
        bc[0] * p0[2] + bc[1] * p1[2] + bc[2] * p2[2],
    ]
}

fn wdw_inverse_distance_cubic_integrated_over_wedge(x: &[f64; 3], b: f64) -> (f64, [f64; 3]) {
    let l = crate::vec3::norm(x);
    let c = 1. / (b * 0.5).tan();
    let a = (l - x[0]) * c - x[1];
    let t = {
        let t = (x[2].abs() / a).atan();
        if t > 0. {
            t
        } else {
            t + core::f64::consts::PI
        }
    };
    let signz = if x[2] < 0. { -1. } else { 1. };
    let w = t * 2. / x[2].abs();
    let d = 1.0 / (x[2] * x[2] + a * a);
    let dwdx = 2. * (1. - x[0] / l) * c * d;
    let dwdy = 2. * (1. - x[2] * c / l) * d;
    let dwdz = -t / (x[2] * x[2]) + a * d / x[2].abs() - x[2].abs() * c * d / l;
    (w, [dwdx, dwdy, dwdz * signz * 2.])
}
pub fn wdw_integral_of_inverse_distance_cubic(
    p0: &[f64; 3],
    p1: &[f64; 3],
    p2: &[f64; 3],
    q: &[f64; 3],
) -> (f64, [f64; 3]) {
    use crate::vec3::Vec3;
    let (vz, _) = unit_normal_area(p0, p1, p2);
    let z = q.sub(p0).dot(&vz);
    let u10 = p0.sub(p1).normalize();
    let u21 = p1.sub(p2).normalize();
    let u02 = p2.sub(p0).normalize();
    //
    let vy0 = vz.cross(&u02);
    let beta0 = u02.dot(&u10).acos();
    let q0 = [q.sub(p0).dot(&u02), q.sub(p0).dot(&vy0), z];
    let (w0, dw0) = wdw_inverse_distance_cubic_integrated_over_wedge(&q0, beta0);
    let dw0dq =
        crate::vec3::add_three_vectors(&u02.scale(dw0[0]), &vy0.scale(dw0[1]), &vz.scale(dw0[2]));
    //
    let vy1 = vz.cross(&u10);
    let beta1 = u10.dot(&u21).acos();
    let q1 = [q.sub(p1).dot(&u10), q.sub(p1).dot(&vy1), z];
    let (w1, dw1) = wdw_inverse_distance_cubic_integrated_over_wedge(&q1, beta1);
    let dw1dq =
        crate::vec3::add_three_vectors(&u10.scale(dw1[0]), &vy1.scale(dw1[1]), &vz.scale(dw1[2]));
    //
    let vy2 = vz.cross(&u21);
    let beta2 = u21.dot(&u02).acos();
    let q2 = [q.sub(p2).dot(&u21), q.sub(p2).dot(&vy2), z];
    let (w2, dw2) = wdw_inverse_distance_cubic_integrated_over_wedge(&q2, beta2);
    let dw2dq =
        crate::vec3::add_three_vectors(&u21.scale(dw2[0]), &vy2.scale(dw2[1]), &vz.scale(dw2[2]));
    //
    let w = core::f64::consts::PI * 2_f64 / z.abs() - w0 - w1 - w2;
    let signz = if z < 0. { -1. } else { 1. };
    let dwdq = crate::vec3::add_three_vectors(&dw0dq, &dw1dq, &dw2dq);
    let dw = vz
        .scale(-signz * core::f64::consts::PI * 2_f64 / (z * z))
        .sub(&dwdq);
    (w, dw)
}

#[test]
fn test_w_inverse_distance_cubic_integrated_over_wedge() {
    use crate::vec3::Vec3;
    use rand::SeedableRng;
    let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    for _ in 0..1000 {
        let x = crate::vec3::sample_unit_cube::<_, f64>(&mut reng).sub(&[0.5, 0.5, 0.5]);
        if x[2].abs() < 0.1 {
            continue;
        }
        let b = core::f64::consts::PI * 0.5; // 90 degree
        let n_r = 200;
        let n_t = 100;
        let rad = 20.;
        //
        let mut val_nmr = 0.;
        for i_r in 0..n_r {
            for i_t in 0..n_t {
                let r = (i_r as f64 + 0.5) * rad / n_r as f64;
                let t = (i_t as f64 + 0.5) * b / n_t as f64;
                let pos = [r * t.cos(), r * t.sin(), 0.];
                let area = rad / n_r as f64 * r * b / n_t as f64;
                val_nmr += area * pos.sub(&x).norm().powi(-3);
            }
        }
        use crate::tri3::wdw_inverse_distance_cubic_integrated_over_wedge;
        let (v_anl, _) = wdw_inverse_distance_cubic_integrated_over_wedge(&x, b);
        assert!((v_anl - val_nmr).abs() < v_anl * 0.1);
        // dbg!(v_anl, val_nmr);
    }
}

pub fn nearest_to_point3<T>(q0: &[T; 3], q1: &[T; 3], q2: &[T; 3], ps: &[T; 3]) -> ([T; 3], T, T)
where
    T: num_traits::Float + Copy + 'static,
    f64: AsPrimitive<T>,
{
    use crate::vec3::Vec3;
    let (n012, _area) = unit_normal_area(q0, q1, q2);
    let pe = ps.add(&n012);
    let v012 = crate::tet::volume(ps, q0, q1, q2);
    if v012.abs() > T::zero() {
        let v0 = crate::tet::volume(ps, q1, q2, &pe);
        let v1 = crate::tet::volume(ps, q2, q0, &pe);
        let v2 = crate::tet::volume(ps, q0, q1, &pe);
        assert!((v0 + v1 + v2).abs() > T::zero());
        let inv_v012 = T::one() / (v0 + v1 + v2);
        let r0 = v0 * inv_v012;
        let r1 = v1 * inv_v012;
        let r2 = T::one() - r0 - r1;
        if r0 >= T::zero() && r1 >= T::zero() && r2 >= T::zero() {
            let nearp = crate::vec3::add_three_vectors(&q0.scale(r0), &q1.scale(r1), &q2.scale(r2));
            return (nearp, r0, r1);
        }
    }
    let r12 = crate::edge3::nearest_to_point3(q1, q2, ps);
    let r20 = crate::edge3::nearest_to_point3(q2, q0, ps);
    let r01 = crate::edge3::nearest_to_point3(q0, q1, ps);
    let d12 = r12.0;
    let d20 = r20.0;
    let d01 = r01.0;
    let r12 = q1.add(&q2.sub(q1).scale(r12.1));
    let r20 = q2.add(&q0.sub(q2).scale(r20.1));
    let r01 = q0.add(&q1.sub(q0).scale(r01.1));
    if d12 < d20 {
        if d12 < d01 {
            // 12 is the smallest
            let r0 = T::zero();
            let r1 = r12.sub(q2).norm() / q1.sub(q2).norm();
            return (r12, r0, r1);
        }
    } else if d20 < d01 {
        // d20 is the smallest
        let r0 = r20.sub(q2).norm() / q0.sub(q2).norm();
        let r1 = T::zero();
        return (r20, r0, r1);
    }
    let r0 = r01.sub(q1).norm() / q0.sub(q1).norm();
    let r1 = T::one() - r0;
    (r01, r0, r1)
}

/// if the triangle share a point, set the point as `p0` and `q0`
pub fn is_intersection_tri3<T>(
    p0: &[T; 3],
    p1: &[T; 3],
    p2: &[T; 3],
    q0: &[T; 3],
    q1: &[T; 3],
    q2: &[T; 3],
) -> Option<([T; 3], [T; 3])>
where
    T: num_traits::Float + Copy,
{
    use crate::vec3::Vec3;
    let sec = |p0: &[T; 3], p1: &[T; 3], dp0: T, dp1: T| {
        let r1 = dp0 / (dp0 - dp1);
        let r0 = T::one() - r1;
        p0.scale(r0).add(&p1.scale(r1))
    };
    let sgn = |v: T| {
        if v == T::zero() {
            1
        } else if v < T::zero() {
            0
        } else {
            2
        }
    };
    let nq = normal(q0, q1, q2);
    // heights of (p0,p1,p2) against the plane span by (q0,q1,q2)
    let np = normal(p0, p1, p2);
    // dbg!(p0,p1,p2,q0,q1,q2);
    let vz = np.cross(&nq);
    //
    let (ps, pe) = {
        let dp0 = p0.sub(q0).dot(&nq);
        let dp1 = p1.sub(q0).dot(&nq);
        let dp2 = p2.sub(q0).dot(&nq);
        let (sp0, sp1, sp2) = (sgn(dp0), sgn(dp1), sgn(dp2));
        if sp0 == 0 && sp1 == 0 && sp2 == 0 {
            return None;
        } // all negative side
        if sp0 == 2 && sp1 == 2 && sp2 == 2 {
            return None;
        } // all positive side
        if sp0 + sp1 + sp2 == 1 || sp0 + sp1 + sp2 == 5 {
            return None;
        } // sharing point but not intersecting
        if sp0 == 1 && sp1 == 1 && sp2 == 1 {
            return None;
        } // degenerate case inside same plane
          // intersection of the lines connecting (p0,p1),(p1,p2),(p2,p0) and the plane span by (q0,q1,q2)
        let mut ap = Vec::<[T; 3]>::with_capacity(2);
        if (sp0 == 0 && sp1 == 2) || (sp0 == 2 && sp1 == 0) {
            ap.push(sec(p0, p1, dp0, dp1));
        }
        if (sp1 == 0 && sp2 == 2) || (sp1 == 2 && sp2 == 0) {
            ap.push(sec(p1, p2, dp1, dp2));
        }
        if (sp2 == 0 && sp0 == 2) || (sp2 == 2 && sp0 == 0) {
            ap.push(sec(p2, p0, dp2, dp0));
        }
        if sp0 == 1 {
            ap.push(*p0);
        }
        if sp1 == 1 {
            ap.push(*p1);
        }
        if sp2 == 1 {
            ap.push(*p2);
        }
        assert_eq!(ap.len(), 2);
        (ap[0], ap[1])
    };
    let zps = ps.dot(&vz);
    let zpe = pe.dot(&vz);
    let (ps, pe, zps, zpe) = if zps > zpe {
        (pe, ps, zpe, zps)
    } else {
        (ps, pe, zps, zpe)
    };
    assert!(zps <= zpe);
    //
    let (qs, qe) = {
        // intersection line between triangle p and plane span by (q0,q1,q2)
        let dq0 = q0.sub(p0).dot(&np); // projection of q0
        let dq1 = q1.sub(p0).dot(&np); // projection of q1
        let dq2 = q2.sub(p0).dot(&np); // projection of q2
        let (sq0, sq1, sq2) = (sgn(dq0), sgn(dq1), sgn(dq2));
        if sq0 == 0 && sq1 == 0 && sq2 == 0 {
            return None;
        }
        if sq0 == 2 && sq1 == 2 && sq2 == 2 {
            return None;
        }
        if sq0 + sq1 + sq2 == 1 || sq0 + sq1 + sq2 == 5 {
            return None;
        } // sharing point
        if sq0 == 1 && sq1 == 1 && sq2 == 1 {
            return None;
        }
        // intersection of the lines connecting (q0,q1),(q1,q2),(q2,q0) and the plane span by (p0,p1,p2)
        let mut aq = Vec::<[T; 3]>::with_capacity(2);
        if (sq0 == 0 && sq1 == 2) || (sq0 == 2 && sq1 == 0) {
            aq.push(sec(q0, q1, dq0, dq1));
        }
        if (sq1 == 0 && sq2 == 2) || (sq1 == 2 && sq2 == 0) {
            aq.push(sec(q1, q2, dq1, dq2));
        }
        if (sq2 == 0 && sq0 == 2) || (sq2 == 2 && sq0 == 0) {
            aq.push(sec(q2, q0, dq2, dq0));
        }
        if sq0 == 1 {
            aq.push(*q0);
        }
        if sq1 == 1 {
            aq.push(*q1);
        }
        if sq2 == 1 {
            aq.push(*q2);
        }
        assert_eq!(aq.len(), 2);
        (aq[0], aq[1])
    };
    let zqs = qs.dot(&vz);
    let zqe = qe.dot(&vz);
    let (qs, qe, zqs, zqe) = if zqs > zqe {
        (qe, qs, zqe, zqs)
    } else {
        (qs, qe, zqs, zqe)
    };
    assert!(zqs <= zqe);
    //
    if zps >= zqe || zqs >= zpe {
        return None;
    } // no overlap or overlap at point
    let s = if zps < zqs { qs } else { ps };
    let e = if zpe < zqe { pe } else { qe };
    Some((s, e))
}

// -------------------------
#[derive(Debug, Copy, Clone)]
pub struct Tri3<'a, Real> {
    pub p0: &'a [Real; 3],
    pub p1: &'a [Real; 3],
    pub p2: &'a [Real; 3],
}

impl<Real> Tri3<'_, Real>
where
    Real: num_traits::Float,
{
    pub fn intersection_against_ray(
        &self,
        ray_org: &[Real; 3],
        ray_dir: &[Real; 3],
    ) -> Option<Real> {
        intersection_against_line(self.p0, self.p1, self.p2, ray_org, ray_dir)
            .filter(|&t| t >= Real::zero())
    }

    pub fn intersection_against_line(
        &self,
        line_org: &[Real; 3],
        line_dir: &[Real; 3],
    ) -> Option<Real> {
        intersection_against_line(self.p0, self.p1, self.p2, line_org, line_dir)
    }

    pub fn area(&self) -> Real {
        area(self.p0, self.p1, self.p2)
    }

    pub fn cog(&self) -> [Real; 3] {
        use crate::vec3::Vec3;
        let one = Real::one();
        let one3rd = one / (one * one * one);
        [
            (self.p0[0] + self.p1[0] + self.p2[0]),
            (self.p0[1] + self.p1[1] + self.p2[1]),
            (self.p0[2] + self.p1[2] + self.p2[2]),
        ]
        .scale(one3rd)
    }

    pub fn normal(&self) -> [Real; 3] {
        normal(self.p0, self.p1, self.p2)
    }

    pub fn unit_normal(&self) -> [Real; 3] {
        let n = normal(self.p0, self.p1, self.p2);
        use crate::vec3::Vec3;
        n.normalize()
    }

    pub fn position_from_barycentric_coordinates(&self, r0: Real, r1: Real) -> [Real; 3] {
        let r2 = Real::one() - r0 - r1;
        [
            r0 * self.p0[0] + r1 * self.p1[0] + r2 * self.p2[0],
            r0 * self.p0[1] + r1 * self.p1[1] + r2 * self.p2[1],
            r0 * self.p0[2] + r1 * self.p1[2] + r2 * self.p2[2],
        ]
    }
}
