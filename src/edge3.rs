use num_traits::AsPrimitive;

pub fn nearest_point3_<T>(
    point_pos: &[T],
    edge_pos0: &[T],
    edge_pos1: &[T]) -> [T; 3]
    where T: num_traits::Float + 'static + Copy + PartialOrd,
          f64: num_traits::AsPrimitive<T>
{
    use crate::vec3;
    let d = [
        edge_pos1[0] - edge_pos0[0],
        edge_pos1[1] - edge_pos0[1],
        edge_pos1[2] - edge_pos0[2]];
    let t = {
        if vec3::dot_(&d, &d) > 1.0e-20_f64.as_() {
            let ps = [
                edge_pos0[0] - point_pos[0],
                edge_pos0[1] - point_pos[1],
                edge_pos0[2] - point_pos[2]];
            let a = vec3::dot_(&d, &d);
            let b = vec3::dot_(&d, &ps);
            let mut r: T = -b / a;
            r = if r < 0_f64.as_() { 0_f64.as_() } else { r };
            r = if r > 1_f64.as_() { 1_f64.as_() } else { r };
            r
        } else {
            0.5_f64.as_()
        }
    };
    [edge_pos0[0] + t * d[0],
        edge_pos0[1] + t * d[1],
        edge_pos0[2] + t * d[2]]
}

/* -------------------------- */

pub fn nearest_to_line3(
    edge_start: &nalgebra::Vector3<f64>,
    edge_end: &nalgebra::Vector3<f64>,
    line_origin: &nalgebra::Vector3<f64>,
    line_direction: &nalgebra::Vector3<f64>) -> (nalgebra::Vector3<f64>, nalgebra::Vector3<f64>)
{
    let (scale, scaled_ratio_edge, _, scaled_nearest_edge, scaled_nearest_line)
        = crate::line3::nearest_to_line3(
        edge_start, &(edge_end - edge_start),
        line_origin, line_direction);
    if scale.abs() < 1.0e-10 { // pararell
        let nearest_edge = (edge_start + edge_end) * 0.5;
        let (nearest_line, _) = crate::line::nearest_to_point(
            &nearest_edge, line_origin, line_direction);
        return (nearest_edge, nearest_line);
    }
    let ratio_edge = scaled_ratio_edge / scale;
    if ratio_edge > 0_f64 && ratio_edge < 1_f64 { // nearst point is inside the segment
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
    q: &nalgebra::Vector3::<f64>,
    p0: &nalgebra::Vector3::<f64>,
    p1: &nalgebra::Vector3::<f64>) -> (f64, nalgebra::Vector3::<f64>) {
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


pub fn nearest_to_edge3(
    p0: &nalgebra::Vector3::<f64>,
    p1: &nalgebra::Vector3::<f64>,
    q0: &nalgebra::Vector3::<f64>,
    q1: &nalgebra::Vector3::<f64>) -> (f64, f64, f64)
{
    let vp = p1 - p0;
    let vq = q1 - q0;
    if vp.cross(&vq).norm() < 1.0e-10 { // handling parallel edge
        let pq0 = p0 - q0;
        let nvp = vp;
        let nvp = nvp.normalize();
        let vert = pq0 - pq0.dot(&nvp) * nvp;
        let dist = vert.norm();
        let lp0 = p0.dot(&nvp);
        let lp1 = p1.dot(&nvp);
        let lq0 = q0.dot(&nvp);
        let lq1 = q1.dot(&nvp);
        let p_min = if lp0 < lp1 { lp0 } else { lp1 };
        let p_max = if lp0 > lp1 { lp0 } else { lp1 };
        let q_min = if lq0 < lq1 { lq0 } else { lq1 };
        let q_max = if lq0 > lq1 { lq0 } else { lq1 };
        let lm;
        if p_max < q_min {
            lm = (p_max + q_min) * 0.5;
        } else if q_max < p_min {
            lm = (q_max + p_min) * 0.5;
        } else if p_max < q_max {
            lm = (p_max + q_min) * 0.5;
        } else {
            lm = (q_max + p_min) * 0.5;
        }
        let ratio_p = (lm - lp0) / (lp1 - lp0);
        let ratio_q = (lm - lq0) / (lq1 - lq0);
        return (dist, ratio_p, ratio_q);
    }
    let t0 = vp.dot(&vp);
    let t1 = vq.dot(&vq);
    let t2 = vp.dot(&vq);
    let t3 = vp.dot(&(q0 - p0));
    let t4 = vq.dot(&(q0 - p0));
    let det = t0 * t1 - t2 * t2;
    let invdet = 1.0 / det;
    let ratio_p = (t1 * t3 - t2 * t4) * invdet;
    let ratio_q = (t2 * t3 - t0 * t4) * invdet;
    let pc = p0 + ratio_p * vp;
    let qc = q0 + ratio_q * vq;
    ((pc - qc).norm(), ratio_p, ratio_q)
}

#[cfg(test)]
mod tests {
    fn numerical(
        q: &nalgebra::Vector3::<f64>,
        p0: &nalgebra::Vector3::<f64>,
        p1: &nalgebra::Vector3::<f64>,
        n: usize,
        p: usize) -> f64 {
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
        for _i in 0..10000 {
            let p0 = crate::vec3::sample_unit_cube();
            let p1 = crate::vec3::sample_unit_cube();
            let q = crate::vec3::sample_unit_cube();
            let len = (p0 - p1).norm();
            let height = crate::tri3::height(&q, &p0, &p1);
            if height < 0.1 { continue; }
            if len < 0.1 { continue; }
            if (p0 - q).norm() < 0.1 { continue; }
            if (p1 - q).norm() < 0.1 { continue; }
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
        let eps = 10e-3;
        for _i in 0..10000 {
            let p0 = crate::vec3::sample_unit_cube();
            let p1 = crate::vec3::sample_unit_cube();
            let q0 = crate::vec3::sample_unit_cube();
            let q1 = crate::vec3::sample_unit_cube();
            let (dist, rp, rq) = crate::edge3::nearest_to_edge3(&p0, &p1, &q0, &q1);
            let vp = p1 - p0;
            let vq = q1 - q0;
            let pc1 = p0 + rp * vp;
            let pc0 = p0 + f64::clamp(rp - eps, 0.0, 1.0) * vp;
            let pc2 = p0 + f64::clamp(rp + eps, 0.0, 1.0) * vp;
            let qc1 = q0 + rq * vq;
            let qc0 = q0 + f64::clamp(rq - eps, 0.0, 1.0) * vq;
            let qc2 = q0 + f64::clamp(rq + eps, 0.0, 1.0) * vq;
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