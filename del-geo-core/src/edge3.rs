//! methods for 3D edge (line segment)

use num_traits::AsPrimitive;
/// trait for 3D edge (line segment)
pub trait Edge3Trait<T> {
    fn length(&self, other: &Self) -> T;
    fn squared_length(&self, other: &Self) -> T;
    fn nearest_to_point3(&self, p1: &Self, point_pos: &Self) -> (T, T);
}

impl<Real> Edge3Trait<Real> for [Real; 3]
where
    Real: num_traits::Float + 'static,
    f64: AsPrimitive<Real>,
{
    fn length(&self, other: &Self) -> Real {
        length(self, other)
    }
    fn squared_length(&self, other: &Self) -> Real {
        squared_length(self, other)
    }
    fn nearest_to_point3(&self, p1: &Self, point_pos: &Self) -> (Real, Real) {
        nearest_to_point3(self, p1, point_pos)
    }
}

// ------------------------------

pub fn length<T>(p0: &[T; 3], p1: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    let x = p0[0] - p1[0];
    let y = p0[1] - p1[1];
    let z = p0[2] - p1[2];
    (x * x + y * y + z * z).sqrt()
}

pub fn squared_length<T>(p0: &[T; 3], p1: &[T; 3]) -> T
where
    T: num_traits::Float,
{
    let x = p0[0] - p1[0];
    let y = p0[1] - p1[1];
    let z = p0[2] - p1[2];
    x * x + y * y + z * z
}

/// `ratio==0` should output `p0`
pub fn position_from_ratio<T>(p0: &[T; 3], p1: &[T; 3], ratio: T) -> [T; 3]
where
    T: num_traits::Float,
{
    let one = T::one();
    [
        (one - ratio) * p0[0] + ratio * p1[0],
        (one - ratio) * p0[1] + ratio * p1[1],
        (one - ratio) * p0[2] + ratio * p1[2],
    ]
}

pub fn lerp<T>(p0: &[T; 3], p1: &[T; 3], s: T) -> [T; 3]
where
    T: num_traits::Float,
{
    [
        p0[0] + s * (p1[0] - p0[0]),
        p0[1] + s * (p1[1] - p0[1]),
        p0[2] + s * (p1[2] - p0[2]),
    ]
}

pub fn wdw_integral_of_inverse_distance_cubic<T>(
    q: &[T; 3],
    p0: &[T; 3],
    p1: &[T; 3],
) -> (T, [T; 3])
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let one = T::one();
    let two = one + one;
    let three = two + one;
    let len = p1.sub(p0).norm();
    let lsinv = one / (len * len);
    // dist^2 = er^2+2br+c
    let d = p0.sub(p1).dot(&q.sub(p0)) * lsinv;
    let a = q.sub(p0).squared_norm() * lsinv - d * d;
    // dist^2 = e{ x^2 + a2}, x = r + d
    // \int 1/sqrt(x^2+a)^3 dx = x/(a\sqrt{a+x^2})
    let f = |x| x / (a * (a + x * x).sqrt());
    let v = (f(d + one) - f(d)) * lsinv;
    //
    let dd = p0.sub(p1).scale(lsinv);
    let da = q.sub(p0).scale(two).scale(lsinv).sub(&dd.scale(two * d));
    // these formula was calculated by WolframAlpha
    let dfdx = |x| one / (a + x * x).powf(three / two);
    let dfda = |x| -(x * (three * a + two * x * x)) / (two * a * a * (a + x * x).powf(three / two));
    let t0 = dd.scale(dfdx(d + one) - dfdx(d));
    let t1 = da.scale(dfda(d + one) - dfda(d));
    let dv = t0.add(&t1);
    (v, dv.scale(lsinv))
}

#[cfg(test)]
mod tests {
    use crate::edge3::position_from_ratio;

    fn numerical(q: &[f64; 3], p0: &[f64; 3], p1: &[f64; 3], n: usize, p: usize) -> f64 {
        use crate::vec3::Vec3;
        use num_traits::Pow;
        let len = p1.sub(p0).norm();
        let mut ret = 0.;
        for i_seg in 0..n {
            let r0 = i_seg as f64 / n as f64;
            let r1 = (i_seg + 1) as f64 / n as f64;
            let pr0q = position_from_ratio(p0, p1, r0).sub(q);
            let pr1q = position_from_ratio(p0, p1, r1).sub(q);
            let dist0 = pr0q.norm();
            let dist1 = pr1q.norm();
            let v0 = 1. / dist0.pow(p as i32);
            let v1 = 1. / dist1.pow(p as i32);
            let v = (v0 + v1) * 0.5;
            ret += v;
        }
        ret *= len / (n as f64);
        ret
    }

    #[test]
    fn test_wdw_integral_of_inverse_distance_cubic() {
        use crate::vec3::Vec3;
        use rand::SeedableRng;
        let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
        for _i in 0..10000 {
            let p0 = crate::vec3::sample_unit_cube(&mut reng);
            let p1 = crate::vec3::sample_unit_cube(&mut reng);
            let q = crate::vec3::sample_unit_cube(&mut reng);
            let len = p0.sub(&p1).norm();
            let height = crate::tri3::height(&p0, &p1, &q);
            if height < 0.1 {
                continue;
            }
            if len < 0.1 {
                continue;
            }
            if p0.sub(&q).norm() < 0.1 {
                continue;
            }
            if p1.sub(&q).norm() < 0.1 {
                continue;
            }
            // dbg!(numerical(&q, &p0, &p1, 10, 3));
            // dbg!(numerical(&q, &p0, &p1, 100, 3));
            // dbg!(numerical(&q, &p0, &p1, 1000, 3));
            let (v0, dv0) = crate::edge3::wdw_integral_of_inverse_distance_cubic(&q, &p0, &p1);
            assert!((v0 - numerical(&q, &p0, &p1, 1000, 3)).abs() < 1.0e-4 * v0.abs());
            let eps = 1.0e-4_f64;
            let qex = [q[0] + eps, q[1], q[2]];
            let qey = [q[0], q[1] + eps, q[2]];
            let qez = [q[0], q[1], q[2] + eps];
            let vx = (numerical(&qex, &p0, &p1, 1000, 3) - v0) / eps;
            let vy = (numerical(&qey, &p0, &p1, 1000, 3) - v0) / eps;
            let vz = (numerical(&qez, &p0, &p1, 1000, 3) - v0) / eps;
            let dv1 = [vx, vy, vz];
            // dbg!(p0, p1, q);
            assert!(dv0.sub(&dv1).norm() < 0.03 * (dv0.norm() + 1.0));
        }
    }
}

// ----------------------------------
// below proximity

pub fn nearest_to_origin3<T>(p0: &[T; 3], p1: &[T; 3]) -> ([T; 3], T, T)
where
    T: num_traits::Float,
{
    let origin = &[T::zero(); 3];
    let (_dist, t) = nearest_to_point3(p0, p1, origin);
    let p = position_from_ratio(p0, p1, t);
    let s0 = T::one() - t;
    (p, s0, t)
}

/// * Returns `(dist, ratio)`
///   - `dist` : distance
///   - `ratio`: ratio
pub fn nearest_to_point3<T>(p0: &[T; 3], p1: &[T; 3], point_pos: &[T; 3]) -> (T, T)
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let zero = T::zero();
    let one = T::one();
    let half = one / (one + one);
    let d = p1.sub(p0);
    let t = {
        if d.dot(&d) > T::epsilon() {
            let ps = std::array::from_fn(|i| p0[i] - point_pos[i]);
            let a = d.dot(&d);
            let b = d.dot(&ps);
            (-b / a).clamp(zero, one)
        } else {
            half
        }
    };
    let p = crate::vec3::axpy(t, &d, p0);
    let dist = length(&p, point_pos);
    (dist, t)
}

pub fn nearest_to_edge3<T>(p0: &[T; 3], p1: &[T; 3], q0: &[T; 3], q1: &[T; 3]) -> (T, T, T)
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let zero = T::zero();
    let one = T::one();
    let half = one / (one + one);
    let vp = p1.sub(p0);
    let vq = q1.sub(q0);
    assert!(vp.norm() > T::zero());
    assert!(vq.norm() > T::zero());
    if vp.cross(&vq).norm() < T::epsilon() {
        // handling parallel edge
        let pq0 = p0.sub(q0);
        let uvp = vp.normalize();
        // a vector vertical to vp and vq and in the plane of vp and vq
        let vert = pq0.sub(&uvp.scale(pq0.dot(&uvp)));
        let dist = vert.norm(); // distance betwen two edges
        let lp0 = p0.dot(&uvp);
        let lp1 = p1.dot(&uvp);
        let lq0 = q0.dot(&uvp);
        let lq1 = q1.dot(&uvp);
        let (lp_min, lp_max, p_min, p_max, rp_min, rp_max) =
            (lp0, lp1, p0, p1, T::zero(), T::one());
        assert!(lp_min < lp_max);
        let (lq_min, lq_max, q_min, q_max, rq_min, rq_max) = if lq0 < lq1 {
            (lq0, lq1, q0, q1, T::zero(), T::one())
        } else {
            (lq1, lq0, q1, q0, T::one(), T::zero())
        };
        if lp_max < lq_min {
            return (p_max.sub(q_min).norm(), rp_max, rq_min);
        }
        if lq_max < lp_min {
            return (q_max.sub(p_min).norm(), rp_min, rq_max);
        }
        let lm_min = lp_min.max(lq_min);
        let lm_max = lp_max.min(lq_max);
        let lm = (lm_min + lm_max) * half;
        let ratio_p = (lm - lp0) / (lp1 - lp0);
        let ratio_q = (lm - lq0) / (lq1 - lq0);
        return (dist, ratio_p, ratio_q);
    }
    let (rp1, rq1) = {
        // line-line intersection
        let t0 = vp.dot(&vp);
        let t1 = vq.dot(&vq);
        let t2 = vp.dot(&vq);
        let t3 = vp.dot(&q0.sub(p0));
        let t4 = vq.dot(&q0.sub(p0));
        let det = t0 * t1 - t2 * t2;
        let invdet = one / det;
        let rp1 = (t1 * t3 - t2 * t4) * invdet;
        let rq1 = (t2 * t3 - t0 * t4) * invdet;
        (rp1, rq1)
    };
    if zero <= rp1 && rp1 <= one && zero <= rq1 && rq1 <= one {
        // both in range
        let pc = p0.add(&vp.scale(rp1));
        let qc = q0.add(&vq.scale(rq1));
        return (pc.sub(&qc).norm(), rp1, rq1);
    }
    if (zero <= rp1 && rp1 <= one) && (rq1 <= zero || one <= rq1) {
        // p in range
        let rq1 = num_traits::clamp(rq1, zero, one);
        let qc = crate::vec3::axpy(rq1, &vq, q0);
        let (dist, rp1) = nearest_to_point3(p0, p1, &qc);
        return (dist, rp1, rq1);
    }
    if (zero <= rq1 && rq1 <= one) && (rp1 <= zero || one <= rp1) {
        // q in range
        let rp1 = num_traits::clamp(rp1, zero, one);
        let pc = crate::vec3::axpy(rp1, &vp, p0);
        let (dist, rq1) = nearest_to_point3(q0, q1, &pc);
        return (dist, rp1, rq1);
    }
    // convex projection technique
    let rp1 = num_traits::clamp(rp1, zero, one);
    let pc = p0.add(&vp.scale(rp1));
    let (_dist, rq1) = nearest_to_point3(q0, q1, &pc);
    let qc = q0.add(&q1.sub(q0).scale(rq1));
    let (_dist, rp1) = nearest_to_point3(p0, p1, &qc);
    let pc = p0.add(&p1.sub(p0).scale(rp1));
    let (dist, rq1) = nearest_to_point3(q0, q1, &pc);
    (dist, rp1, rq1)
}

#[test]
fn test_nearest_to_edge3() {
    use crate::vec3::Vec3;
    use crate::vec3::axpy;
    use rand::SeedableRng;
    let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    let eps = 1.0e-4;
    for _i in 0..10000 {
        let p0 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let p1 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let q0 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let q1 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let (dist, rp, rq) = nearest_to_edge3(&p0, &p1, &q0, &q1);
        {
            let [sqdist, rp0, rq0] = nearest_to_edge3_accurate(&p0, &p1, &q0, &q1);
            println!("{} {} {}", dist * dist - sqdist, rp - rp0, rq - rq0);
        }
        //
        let vp = p1.sub(&p0);
        //let pc0 = p0 + f64::clamp(rp - eps, 0.0, 1.0) * vp;
        let pc0 = axpy(f64::clamp(rp - eps, 0.0, 1.0), &vp, &p0);
        let pc1 = axpy(rp, &vp, &p0);
        let pc2 = axpy(f64::clamp(rp + eps, 0.0, 1.0), &vp, &p0);
        //
        let vq = q1.sub(&q0);
        let qc0 = axpy(f64::clamp(rq - eps, 0.0, 1.0), &vq, &q0);
        let qc1 = axpy(rq, &vq, &q0);
        let qc2 = axpy(f64::clamp(rq + eps, 0.0, 1.0), &vq, &q0);
        assert!((dist - (pc1.sub(&qc1)).norm()).abs() < 1.0e-5);
        assert!(dist <= pc0.sub(&qc0).norm());
        assert!(dist <= pc0.sub(&qc1).norm());
        assert!(dist <= pc0.sub(&qc2).norm());
        assert!(dist <= pc1.sub(&qc0).norm());
        assert!(dist <= pc1.sub(&qc2).norm());
        assert!(dist <= pc2.sub(&qc0).norm());
        assert!(dist <= pc2.sub(&qc1).norm());
        assert!(dist <= pc2.sub(&qc2).norm());
    }
}

/// Closest-point barycentric coefficients of two segments:
///
/// A(s) = ea0 + s * (ea1 - ea0)
/// B(t) = eb0 + t * (eb1 - eb0)
///
/// Returns:
///
/// [sqdist, param_a, param_b]
pub fn nearest_to_edge3_accurate<T: num_traits::Float>(
    ea0: &[T; 3],
    ea1: &[T; 3],
    eb0: &[T; 3],
    eb1: &[T; 3],
) -> [T; 3] {
    let clamp_unit = |v: T| {
        if v > T::zero() {
            if v < T::one() { v } else { T::one() }
        } else {
            T::zero()
        }
    };

    use crate::vec3::{cross, dot, sub};
    let r0 = sub(ea1, ea0);
    let r1 = sub(eb1, eb0);
    let d = sub(ea0, eb0);

    let a = dot(&r0, &r0);
    let e = dot(&r1, &r1);
    let b = dot(&r0, &r1);
    let c = dot(&r0, &d);
    let f = dot(&r1, &d);

    let zero = T::zero();
    let one = T::one();

    // If an edge has collapsed to a point, every parameter on that
    // edge represents the same point.
    let inv_a = if a > zero { one / a } else { zero };
    let inv_e = if e > zero { one / e } else { zero };

    // Interior stationary-point candidate.
    //
    // Use |r0 x r1|^2 instead of a*e - b*b for better numerical
    // conditioning for nearly parallel edges.
    let n = cross(&r0, &r1);
    let nn = dot(&n, &n);

    let inv_nn = if nn > zero { one / nn } else { zero };

    let r1_cross_d = cross(&r1, &d);

    let s_stat = clamp_unit(dot(&n, &r1_cross_d) * inv_nn);

    let t_stat = clamp_unit((b * s_stat + f) * inv_e);

    // Four boundary candidates.
    //
    // s = 0: project ea0 onto edge B
    // s = 1: project ea1 onto edge B
    // t = 0: project eb0 onto edge A
    // t = 1: project eb1 onto edge A
    let t_a0 = clamp_unit(f * inv_e);
    let t_a1 = clamp_unit((f + b) * inv_e);

    let s_b0 = clamp_unit(-c * inv_a);
    let s_b1 = clamp_unit((b - c) * inv_a);

    let cand_s = [s_stat, zero, one, s_b0, s_b1];
    let cand_t = [t_stat, t_a0, t_a1, zero, one];

    let mut best_s = zero;
    let mut best_t = zero;
    let mut best_sqdist = T::max_value();

    for i in 0..5 {
        let s = cand_s[i];
        let t = cand_t[i];

        // Difference between the two candidate closest points:
        //
        // A(s) - B(t)
        // = d + s*r0 - t*r1
        let v = [
            d[0] + s * r0[0] - t * r1[0],
            d[1] + s * r0[1] - t * r1[1],
            d[2] + s * r0[2] - t * r1[2],
        ];

        let sqdist = dot(&v, &v);
        if sqdist < best_sqdist {
            best_sqdist = sqdist;
            best_s = s;
            best_t = t;
        }
    }

    [best_sqdist, best_s, best_t]
}

/// the two edges need to be co-planar
pub fn intersection_edge3_when_coplanar<T>(
    p0: &[T; 3],
    p1: &[T; 3],
    q0: &[T; 3],
    q1: &[T; 3],
) -> Option<(T, T, T, T)>
where
    T: num_traits::Float + Copy + 'static,
    f64: AsPrimitive<T>,
{
    use crate::vec3::Vec3;
    let n = {
        let n0 = p1.sub(p0).cross(&q0.sub(p0));
        let n1 = p1.sub(p0).cross(&q1.sub(p0));
        if n0.squared_norm() < n1.squared_norm() {
            n1
        } else {
            n0
        }
    };
    let p2 = p0.add(&n);
    let rq1 = crate::tet::volume(p0, p1, &p2, q0);
    let rq0 = crate::tet::volume(p0, p1, &p2, q1);
    let rp1 = crate::tet::volume(q0, q1, &p2, p0);
    let rp0 = crate::tet::volume(q0, q1, &p2, p1);
    if (rp0 - rp1).abs() <= T::zero() {
        return None;
    }
    if (rq0 - rq1).abs() <= T::zero() {
        return None;
    }
    let t = T::one() / (rp0 - rp1);
    let (rp0, rp1) = (rp0 * t, -rp1 * t);
    let t = T::one() / (rq0 - rq1);
    let (rq0, rq1) = (rq0 * t, -rq1 * t);
    Some((rp0, rp1, rq0, rq1))
}

pub fn squared_distance_against_edge(
    a0: &[f64; 3],
    a1: &[f64; 3],
    b0: &[f64; 3],
    b1: &[f64; 3],
) -> f64 {
    use crate::vec3::{dot, sub};
    const EPS_LEN_SQ: f64 = 1.0e-12;
    const EPS_DENOM: f64 = 1.0e-12;
    let d1 = sub(a1, a0);
    let d2 = sub(b1, b0);
    let r = sub(a0, b0);

    let a = dot(&d1, &d1);
    let e = dot(&d2, &d2);
    let f = dot(&d2, &r);

    let (s, t);
    if a < EPS_LEN_SQ && e < EPS_LEN_SQ {
        let d = sub(a0, b0);
        return dot(&d, &d);
    } else if a < EPS_LEN_SQ {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else if e < EPS_LEN_SQ {
        t = 0.0;
        s = (-dot(&d1, &r) / a).clamp(0.0, 1.0);
    } else {
        let b_val = dot(&d1, &d2);
        let c = dot(&d1, &r);
        // Gram determinant |d1 x d2|^2, a quartic (L^4) quantity, not a length^2.
        let denom = a * e - b_val * b_val;
        let s_init = if denom.abs() > EPS_DENOM {
            ((b_val * f - c * e) / denom).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let t_init = (b_val * s_init + f) / e;

        if t_init < 0.0 {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else if t_init > 1.0 {
            t = 1.0;
            s = ((b_val - c) / a).clamp(0.0, 1.0);
        } else {
            s = s_init;
            t = t_init;
        }
    }

    let closest_a = [a0[0] + s * d1[0], a0[1] + s * d1[1], a0[2] + s * d1[2]];
    let closest_b = [b0[0] + t * d2[0], b0[1] + t * d2[1], b0[2] + t * d2[2]];
    let diff = sub(&closest_a, &closest_b);
    dot(&diff, &diff)
}

pub fn squared_distance_against_point(e0: &[f64; 3], e1: &[f64; 3], p: &[f64; 3]) -> f64 {
    use crate::vec3::{dot, sub};
    const EPS_LEN_SQ: f64 = 1.0e-12;
    // const EPS_DENOM: f64 = 1.0e-12;
    let edge = sub(e1, e0);
    let edge_len_sq = dot(&edge, &edge);
    if edge_len_sq < EPS_LEN_SQ {
        let d = sub(p, e0);
        return dot(&d, &d);
    }
    let t = (dot(&sub(p, e0), &edge) / edge_len_sq).clamp(0.0, 1.0);
    let closest = [
        e0[0] + t * edge[0],
        e0[1] + t * edge[1],
        e0[2] + t * edge[2],
    ];
    let d = sub(p, &closest);
    dot(&d, &d)
}

// end: proximity
// ----------------------------------------------

#[allow(clippy::type_complexity)]
pub fn wdwddw_squared_length_difference<T>(
    node2xyz_def: &[[T; 3]; 2],
    stiffness: T,
    edge_length_ini: T,
) -> (T, [[T; 3]; 2], [[[T; 9]; 2]; 2])
where
    T: num_traits::Float,
{
    use crate::mat3_col_major::Mat3ColMajor;
    use crate::vec3::Vec3;
    //
    let one = T::one();
    let half = one / (one + one);
    let v = node2xyz_def[0].sub(&node2xyz_def[1]);
    let l = v.norm();
    let c = edge_length_ini - l;
    let dw = [v.scale(-c * stiffness / l), v.scale(c * stiffness / l)];
    let m = {
        let mvv = crate::mat3_col_major::from_scaled_outer_product(one, &v, &v);
        let t0 = stiffness * edge_length_ini / (l * l * l);
        let t1 = stiffness * (l - edge_length_ini) / l;
        let t2 = crate::mat3_col_major::from_identity().scale(t1);
        mvv.scale(t0).add(&t2)
    };
    let ddw = [[m, m.scale(-one)], [m.scale(-one), m]];
    let w = half * stiffness * c * c;
    (w, dw, ddw)
}

pub fn w_squared_length_difference<T>(
    node2xyz_def: &[[T; 3]; 2],
    stiffness: T,
    edge_length_ini: T,
) -> T
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    //
    let one = T::one();
    let half = one / (one + one);
    let v = node2xyz_def[0].sub(&node2xyz_def[1]);
    let l = v.norm();
    let c = edge_length_ini - l;
    half * stiffness * c * c
}
