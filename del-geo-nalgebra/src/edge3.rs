//! functions for 3D edge (line segment)

use num_traits::AsPrimitive;

/// return  distance and parameter
pub fn nearest_to_point3<T>(
    edge_pos0: &nalgebra::Vector3<T>,
    edge_pos1: &nalgebra::Vector3<T>,
    point_pos: &nalgebra::Vector3<T>,
) -> (T, T)
where
    T: nalgebra::RealField + 'static + Copy + PartialOrd,
    f64: AsPrimitive<T>,
{
    let d = edge_pos1 - edge_pos0;
    let dsq = d.norm_squared();
    let t = if dsq != T::zero() {
        let r: T = -d.dot(&(edge_pos0 - point_pos)) / dsq;
        nalgebra::clamp(r, T::zero(), T::one())
    } else {
        0.5_f64.as_()
    };
    let distance = (edge_pos0 + d.scale(t) - point_pos).norm();
    (distance, t)
}

/// the two edges need to be co-planar
pub fn intersection_edge3_when_coplanar<T>(
    p0: &nalgebra::Vector3<T>,
    p1: &nalgebra::Vector3<T>,
    q0: &nalgebra::Vector3<T>,
    q1: &nalgebra::Vector3<T>,
) -> Option<(T, T, T, T)>
where
    T: nalgebra::RealField + Copy + 'static,
    f64: AsPrimitive<T>,
{
    let n = {
        let n0 = (p1 - p0).cross(&(q0 - p0));
        let n1 = (p1 - p0).cross(&(q1 - p0));
        if n0.norm_squared() < n1.norm_squared() {
            n1
        } else {
            n0
        }
    };
    let p2 = p0 + n;
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
    let (rp0, rp1) = (rp0.scale(t), -rp1.scale(t));
    let t = T::one() / (rq0 - rq1);
    let (rq0, rq1) = (rq0.scale(t), -rq1.scale(t));
    Some((rp0, rp1, rq0, rq1))
}

pub fn nearest_to_line3(
    edge_start: &nalgebra::Vector3<f64>,
    edge_end: &nalgebra::Vector3<f64>,
    line_origin: &nalgebra::Vector3<f64>,
    line_direction: &nalgebra::Vector3<f64>,
) -> (nalgebra::Vector3<f64>, nalgebra::Vector3<f64>) {
    let (scale, scaled_ratio_edge, _, scaled_nearest_edge, scaled_nearest_line) =
        crate::line3::nearest_to_line3(
            edge_start,
            &(edge_end - edge_start),
            line_origin,
            line_direction,
        );
    if scale.abs() < 1.0e-10 {
        // pararell
        let nearest_edge = (edge_start + edge_end) * 0.5;
        let (nearest_line, _) =
            crate::line::nearest_to_point(&nearest_edge, line_origin, line_direction);
        return (nearest_edge, nearest_line);
    }
    let ratio_edge = scaled_ratio_edge / scale;
    if ratio_edge > 0_f64 && ratio_edge < 1_f64 {
        // nearst point is inside the segment
        let nearest_edge = scaled_nearest_edge / scale;
        let nearest_line = scaled_nearest_line / scale;
        return (nearest_edge, nearest_line);
    }
    //
    let (p1, _) = crate::line::nearest_to_point(edge_start, line_origin, line_direction);
    let (p2, _) = crate::line::nearest_to_point(edge_end, line_origin, line_direction);
    let dist1 = (p1 - edge_start).norm();
    let dist2 = (p2 - edge_end).norm();
    if dist1 < dist2 {
        let nearest_edge = edge_start;
        let nearest_line = p1;
        return (*nearest_edge, nearest_line);
    }
    let nearest_edge = edge_end;
    let nearest_line = p2;
    (*nearest_edge, nearest_line)
}

pub fn wdw_integral_of_inverse_distance_cubic(
    q: &nalgebra::Vector3<f64>,
    p0: &nalgebra::Vector3<f64>,
    p1: &nalgebra::Vector3<f64>,
) -> (f64, nalgebra::Vector3<f64>) {
    let len = (p1 - p0).norm();
    let lsinv = 1.0 / (len * len);
    // dist^2 = er^2+2br+c
    let d = (p0 - p1).dot(&(q - p0)) * lsinv;
    let a = (q - p0).norm_squared() * lsinv - d * d;
    // dist^2 = e{ x^2 + a2}, x = r + d
    // \int 1/sqrt(x^2+a)^3 dx = x/(a\sqrt{a+x^2})
    let f = |x: f64| x / (a * (a + x * x).sqrt());
    let v = (f(d + 1.) - f(d)) * lsinv;
    //
    let dd = (p0 - p1) * lsinv;
    let da = (q - p0).scale(2_f64) * lsinv - 2.0 * d * dd;
    // these formula was calculated by WolframAlpha
    let dfdx = |x: f64| 1_f64 / (a + x * x).powf(1.5);
    let dfda = |x: f64| -(x * (3. * a + 2. * x * x)) / (2. * a * a * (a + x * x).powf(1.5));
    let dv = (dfdx(d + 1.) - dfdx(d)) * dd + (dfda(d + 1.) - dfda(d)) * da;
    (v, dv * lsinv)
}

pub fn nearest_to_edge3<T>(
    p0: &nalgebra::Vector3<T>,
    p1: &nalgebra::Vector3<T>,
    q0: &nalgebra::Vector3<T>,
    q1: &nalgebra::Vector3<T>,
) -> (T, T, T)
where
    T: nalgebra::RealField + Copy + 'static,
    f64: AsPrimitive<T>,
{
    use nalgebra::clamp;
    let vp = p1 - p0;
    let vq = q1 - q0;
    assert!(vp.norm() > T::zero());
    assert!(vq.norm() > T::zero());
    if vp.cross(&vq).norm() < T::default_epsilon() {
        // handling parallel edge
        let pq0 = p0 - q0;
        let uvp = vp.normalize();
        // a vector vertical to vp and vq and in the plane of vp and vq
        let vert = pq0 - uvp.scale(pq0.dot(&uvp));
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
            return ((p_max - q_min).norm(), rp_max, rq_min);
        }
        if lq_max < lp_min {
            return ((q_max - p_min).norm(), rp_min, rq_max);
        }
        let lm_min = lp_min.max(lq_min);
        let lm_max = lp_max.min(lq_max);
        let lm = (lm_min + lm_max) * 0.5f64.as_();
        let ratio_p = (lm - lp0) / (lp1 - lp0);
        let ratio_q = (lm - lq0) / (lq1 - lq0);
        return (dist, ratio_p, ratio_q);
    }
    let (rp1, rq1) = {
        // line-line intersection
        let t0 = vp.dot(&vp);
        let t1 = vq.dot(&vq);
        let t2 = vp.dot(&vq);
        let t3 = vp.dot(&(q0 - p0));
        let t4 = vq.dot(&(q0 - p0));
        let det = t0 * t1 - t2 * t2;
        let invdet = T::one() / det;
        let rp1 = (t1 * t3 - t2 * t4) * invdet;
        let rq1 = (t2 * t3 - t0 * t4) * invdet;
        (rp1, rq1)
    };
    if T::zero() <= rp1 && rp1 <= T::one() && T::zero() <= rq1 && rq1 <= T::one() {
        // both in range
        let pc = p0 + vp.scale(rp1);
        let qc = q0 + vq.scale(rq1);
        return ((pc - qc).norm(), rp1, rq1);
    }
    if (T::zero() <= rp1 && rp1 <= T::one()) && (rq1 <= T::zero() || T::one() <= rq1) {
        // p in range
        let rq1 = clamp(rq1, T::zero(), T::one());
        let qc = q0 + vq.scale(rq1);
        let (dist, rp1) = nearest_to_point3(p0, p1, &qc);
        return (dist, rp1, rq1);
    }
    if (T::zero() <= rq1 && rq1 <= T::one()) && (rp1 <= T::zero() || T::one() <= rp1) {
        // q in range
        let rp1 = clamp(rp1, T::zero(), T::one());
        let pc = p0 + vp.scale(rp1);
        let (dist, rq1) = nearest_to_point3(q0, q1, &pc);
        return (dist, rp1, rq1);
    }
    // convex projection technique
    let rp1 = clamp(rp1, T::zero(), T::one());
    let pc = p0 + vp.scale(rp1);
    let (_dist, rq1) = nearest_to_point3(q0, q1, &pc);
    let qc = q0 + (q1 - q0).scale(rq1);
    let (_dist, rp1) = nearest_to_point3(p0, p1, &qc);
    let pc = p0 + (p1 - p0).scale(rp1);
    let (dist, rq1) = nearest_to_point3(q0, q1, &pc);
    (dist, rp1, rq1)
}

#[cfg(test)]
mod tests {
    fn numerical(
        q: &nalgebra::Vector3<f64>,
        p0: &nalgebra::Vector3<f64>,
        p1: &nalgebra::Vector3<f64>,
        n: usize,
        p: usize,
    ) -> f64 {
        use num_traits::Pow;
        let len = (p1 - p0).norm();
        let mut ret = 0.;
        for i_seg in 0..n {
            let r0 = i_seg as f64 / n as f64;
            let r1 = (i_seg + 1) as f64 / n as f64;
            let pr0q = p0.scale(1. - r0) + p1.scale(r0) - q;
            let pr1q = p0.scale(1. - r1) + p1.scale(r1) - q;
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
        use rand::SeedableRng;
        let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
        for _i in 0..10000 {
            let p0 = crate::vec3::sample_unit_cube(&mut reng);
            let p1 = crate::vec3::sample_unit_cube(&mut reng);
            let q = crate::vec3::sample_unit_cube(&mut reng);
            let len = (p0 - p1).norm();
            let height = crate::tri3::height(&p0, &p1, &q);
            if height < 0.1 {
                continue;
            }
            if len < 0.1 {
                continue;
            }
            if (p0 - q).norm() < 0.1 {
                continue;
            }
            if (p1 - q).norm() < 0.1 {
                continue;
            }
            // dbg!(numerical(&q, &p0, &p1, 10, 3));
            // dbg!(numerical(&q, &p0, &p1, 100, 3));
            // dbg!(numerical(&q, &p0, &p1, 1000, 3));
            let (v0, dv0) = crate::edge3::wdw_integral_of_inverse_distance_cubic(&q, &p0, &p1);
            assert!((v0 - numerical(&q, &p0, &p1, 1000, 3)).abs() < 1.0e-4 * v0.abs());
            let eps = 1.0e-4_f64;
            let qex = nalgebra::Vector3::new(q.x + eps, q.y, q.z);
            let qey = nalgebra::Vector3::new(q.x, q.y + eps, q.z);
            let qez = nalgebra::Vector3::new(q.x, q.y, q.z + eps);
            let vx = (numerical(&qex, &p0, &p1, 1000, 3) - v0) / eps;
            let vy = (numerical(&qey, &p0, &p1, 1000, 3) - v0) / eps;
            let vz = (numerical(&qez, &p0, &p1, 1000, 3) - v0) / eps;
            let dv1 = nalgebra::Vector3::<f64>::new(vx, vy, vz);
            // dbg!(p0, p1, q);
            assert!((dv0 - dv1).norm() < 0.03 * (dv0.norm() + 1.0));
        }
    }

    #[test]
    fn test_distance() {
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
            let vp = p1 - p0;
            //let pc0 = p0 + f64::clamp(rp - eps, 0.0, 1.0) * vp;
            let pc0 = p0 + f64::clamp(rp - eps, 0.0, 1.0) * vp;
            let pc1 = p0 + rp * vp;
            let pc2 = p0 + f64::clamp(rp + eps, 0.0, 1.0) * vp;
            //
            let vq = q1 - q0;
            let qc0 = q0 + f64::clamp(rq - eps, 0.0, 1.0) * vq;
            let qc1 = q0 + rq * vq;
            let qc2 = q0 + f64::clamp(rq + eps, 0.0, 1.0) * vq;
            assert!((dist - (pc1 - qc1).norm()).abs() < 1.0e-5);
            assert!(dist <= (pc0 - qc0).norm());
            assert!(dist <= (pc0 - qc1).norm());
            assert!(dist <= (pc0 - qc2).norm());
            assert!(dist <= (pc1 - qc0).norm());
            assert!(dist <= (pc1 - qc2).norm());
            assert!(dist <= (pc2 - qc0).norm());
            assert!(dist <= (pc2 - qc1).norm());
            assert!(dist <= (pc2 - qc2).norm());
        }
    }
}
