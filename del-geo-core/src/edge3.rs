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

/// `ratio==1` should output `p0`
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

pub fn wdw_integral_of_inverse_distance_cubic(
    q: &[f64; 3],
    p0: &[f64; 3],
    p1: &[f64; 3],
) -> (f64, [f64; 3]) {
    use crate::vec3::Vec3;
    let len = p1.sub(p0).norm();
    let lsinv = 1.0 / (len * len);
    // dist^2 = er^2+2br+c
    let d = p0.sub(p1).dot(&q.sub(p0)) * lsinv;
    let a = q.sub(p0).squared_norm() * lsinv - d * d;
    // dist^2 = e{ x^2 + a2}, x = r + d
    // \int 1/sqrt(x^2+a)^3 dx = x/(a\sqrt{a+x^2})
    let f = |x: f64| x / (a * (a + x * x).sqrt());
    let v = (f(d + 1.) - f(d)) * lsinv;
    //
    let dd = p0.sub(p1).scale(lsinv);
    let da = q.sub(p0).scale(2_f64).scale(lsinv).sub(&dd.scale(2.0 * d));
    // these formula was calculated by WolframAlpha
    let dfdx = |x: f64| 1_f64 / (a + x * x).powf(1.5);
    let dfda = |x: f64| -(x * (3. * a + 2. * x * x)) / (2. * a * a * (a + x * x).powf(1.5));
    let t0 = dd.scale(dfdx(d + 1.) - dfdx(d));
    let t1 = da.scale(dfda(d + 1.) - dfda(d));
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
    use crate::vec3::axpy;
    use crate::vec3::Vec3;
    use rand::SeedableRng;
    let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    let eps = 1.0e-4;
    for _i in 0..10000 {
        let p0 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let p1 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let q0 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let q1 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let (dist, rp, rq) = crate::edge3::nearest_to_edge3(&p0, &p1, &q0, &q1);
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
