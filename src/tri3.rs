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
    vec3::squared_norm_(&na).sqrt() * 0.5_f64.as_()
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
    let pvec = vec3::cross_(ray_dir, &edge2);
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
    Some(t)
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
    } else if d20 < d01 { // d20 is the smallest
        let nearp = [r20[0], r20[1], r20[2]];
        let r0 = vec3::distance_(&nearp, q2) / vec3::distance_(q0, q2);
        let r1 = 0_f32;
        return (nearp, r0, r1);
    }
    let nearp = [r01[0], r01[1], r01[2]];
    let r0 = vec3::distance_(&nearp, q1) / vec3::distance_(q0, q1);
    let r1 = 1_f32 - r0;
    (nearp, r0, r1)
}

// above: w/o nalgebra
/* ------------------------------------ */
// below: w/ nalgebra

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

pub fn unit_normal<T>(
    p0: &nalgebra::Vector3::<T>,
    p1: &nalgebra::Vector3::<T>,
    p2: &nalgebra::Vector3::<T>) -> nalgebra::Vector3::<T>
    where T: nalgebra::RealField
{
    let n = (p1 - p0).cross(&(p2 - p0));
    n.normalize()
}

pub fn area<T>(
    p0: &nalgebra::Vector3::<T>,
    p1: &nalgebra::Vector3::<T>,
    p2: &nalgebra::Vector3::<T>) -> T
    where T: nalgebra::RealField + 'static + Copy,
          f64: AsPrimitive<T>
{
    (p1 - p0).cross(&(p2 - p0)).norm() * 0.5_f64.as_()
}

fn wdw_inverse_distance_cubic_integrated_over_wedge(
    x: nalgebra::Vector3::<f64>,
    b: f64) -> (f64, nalgebra::Vector3::<f64>)
{
    let l = x.norm();
    let c = 1. / (b * 0.5).tan();
    let a = (l - x.x) * c - x.y;
    let t = {
        let t = (x.z.abs() / a).atan();
        if t > 0. { t } else { t + core::f64::consts::PI }
    };
    let signz = if x.z < 0. { -1. } else { 1. };
    let w = t * 2. / x.z.abs();
    let d = 1.0 / (x.z * x.z + a * a);
    let dwdx = 2. * (1. - x.x / l) * c * d;
    let dwdy = 2. * (1. - x.y * c / l) * d;
    let dwdz = -t / (x.z * x.z) + a * d / x.z.abs() - x.z.abs() * c * d / l;
    (w, nalgebra::Vector3::<f64>::new(dwdx, dwdy, dwdz * signz * 2.))
}

pub fn wdw_integral_of_inverse_distance_cubic(
    p0: &nalgebra::Vector3::<f64>,
    p1: &nalgebra::Vector3::<f64>,
    p2: &nalgebra::Vector3::<f64>,
    q: &nalgebra::Vector3::<f64>) -> (f64, nalgebra::Vector3::<f64>)
{
    let vz = unit_normal(p0, p1, p2);
    let z = (q - p0).dot(&vz);
    let u10 = (p0 - p1).normalize();
    let u21 = (p1 - p2).normalize();
    let u02 = (p2 - p0).normalize();
    //
    let vy0 = vz.cross(&u02);
    let beta0 = u02.dot(&u10).acos();
    let q0 = nalgebra::Vector3::<f64>::new((q - p0).dot(&u02), (q - p0).dot(&vy0), z);
    let (w0, dw0) = wdw_inverse_distance_cubic_integrated_over_wedge(q0, beta0);
    let dw0dq = u02.scale(dw0.x) + vy0.scale(dw0.y) + vz.scale(dw0.z);
    //
    let vy1 = vz.cross(&u10);
    let beta1 = u10.dot(&u21).acos();
    let q1 = nalgebra::Vector3::<f64>::new((q - p1).dot(&u10), (q - p1).dot(&vy1), z);
    let (w1, dw1) = wdw_inverse_distance_cubic_integrated_over_wedge(q1, beta1);
    let dw1dq = u10.scale(dw1.x) + vy1.scale(dw1.y) + vz.scale(dw1.z);
    //
    let vy2 = vz.cross(&u21);
    let beta2 = u21.dot(&u02).acos();
    let q2 = nalgebra::Vector3::<f64>::new((q - p2).dot(&u21), (q - p2).dot(&vy2), z);
    let (w2, dw2) = wdw_inverse_distance_cubic_integrated_over_wedge(q2, beta2);
    let dw2dq = u21.scale(dw2.x) + vy2.scale(dw2.y) + vz.scale(dw2.z);
    //
    let w = core::f64::consts::PI * 2_f64 / z.abs() - w0 - w1 - w2;
    let signz = if z < 0. { -1. } else { 1. };
    let dw = -vz.scale(signz * core::f64::consts::PI * 2_f64 / (z * z))
        - dw0dq - dw1dq - dw2dq;
    (w, dw)
}

pub fn integrate_numerically<F>(
    p0: &nalgebra::Vector3::<f64>,
    p1: &nalgebra::Vector3::<f64>,
    p2: &nalgebra::Vector3::<f64>,
    integrand: F,
    n: usize) -> f64
    where F: Fn(f64, f64) -> f64
{
    let area = crate::tri3::area(p0, p1, p2);
    let jacobian = area / (n * n) as f64;
    let mut val_num = 0_f64;
    for i in 0..n {
        for j in 0..i * 2 + 1 {
            let j0 = j / 2;
            let (u, v) = match j % 2 {
                0 => {
                    let v = (n - i) + (n - i - 1) * 2;
                    let u = j0 * 2 + j0 + 1;
                    (u, v)
                }
                1 => {
                    let v = (n - i) * 2 + (n - i - 1);
                    let u = (j0 + 1) * 2 + j0;
                    (u, v)
                }
                _ => { panic!(); }
            };
            let (u, v) = (u as f64 / (n * 3) as f64, v as f64 / (n * 3) as f64);
            let dist = integrand(u, v);
            val_num += jacobian * dist;
        }
    }
    val_num
}

#[cfg(test)]
mod tests {
    use rand::Rng;
    use crate::tri3::{area, height, integrate_numerically, wdw_integral_of_inverse_distance_cubic};

    #[test]
    fn test_w_inverse_distance_cubic_integrated_over_wedge() {
        for _ in 0..1000 {
            let x = crate::vec3::sample_unit_cube() - nalgebra::Vector3::new(0.5, 0.5, 0.5);
            if x.z.abs() < 0.1 { continue; }
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
                    let pos = nalgebra::Vector3::<f64>::new(r * t.cos(), r * t.sin(), 0.);
                    let area = rad / n_r as f64 * r * b / n_t as f64;
                    val_nmr += area * (pos - x).norm().powi(-3);
                }
            }
            use crate::tri3::wdw_inverse_distance_cubic_integrated_over_wedge;
            let (v_anl, _) = wdw_inverse_distance_cubic_integrated_over_wedge(x, b);
            // dbg!(v_anl, val_nmr);
            assert!((v_anl - val_nmr).abs() < v_anl * 0.1);
        }
    }

    #[test]
    fn test_dw_inverse_distance_cubic_integrated_over_wedge() {
        use crate::tri3::wdw_inverse_distance_cubic_integrated_over_wedge;
        // let x0 = nalgebra::Vector3::<f64>::new(10., 1., -2.);
        for _ in 0..100000 {
            let x0 = crate::vec3::sample_unit_cube() - nalgebra::Vector3::new(0.5, 0.5, 0.5);
            if x0.z.abs() < 0.2 { continue; }
            // dbg!(x0);
            let b = core::f64::consts::PI * (rand::thread_rng().gen::<f64>() * 0.8 + 0.1); // 90 degree
            let eps = 1.0e-4_f64;
            let (w0, dw) = wdw_inverse_distance_cubic_integrated_over_wedge(x0, b);
            let x1x = x0 + nalgebra::Vector3::<f64>::new(eps, 0., 0.);
            let (w1x, _) = wdw_inverse_distance_cubic_integrated_over_wedge(x1x, b);
            // dbg!((w1x-w0)/eps, dw.x);
            let x1y = x0 + nalgebra::Vector3::<f64>::new(0., eps, 0.);
            let (w1y, _) = wdw_inverse_distance_cubic_integrated_over_wedge(x1y, b);
            // dbg!((w1y-w0)/eps, dw.y);
            let x1z = x0 + nalgebra::Vector3::<f64>::new(0., 0., eps);
            let (w1z, _) = wdw_inverse_distance_cubic_integrated_over_wedge(x1z, b);
            //dbg!((w1z-w0)/eps, dw.z);
            // dbg!(b,w0,w1x,w1y,w1z);
            // dbg!(((w1y - w0) / eps - dw.y).abs(), 1.0e-2 * dw.y.abs());
            assert!(((w1x - w0) / eps - dw.x).abs() < 2.0e-2 * (dw.x.abs() + 1.0e-1));
            assert!(((w1y - w0) / eps - dw.y).abs() < 3.0e-2 * (dw.y.abs() + 5.0e-1));
            assert!(((w1z - w0) / eps - dw.z).abs() < 2.0e-2 * (dw.z.abs() + 1.0e-1));
        }
    }

    #[test]
    fn test_w_inverse_distance_cubic_integrated() {
        for _ in 0..1000 {
            let p0 = crate::vec3::sample_unit_cube();
            let p1 = crate::vec3::sample_unit_cube();
            let p2  = crate::vec3::sample_unit_cube();
            let q  = crate::vec3::sample_unit_cube();
            {
                let area = area(&p0, &p1, &p2);
                let h0 = height(&p0, &p1, &p2);
                let h1 = height(&p1, &p2, &p0);
                let h2 = height(&p2, &p0, &p1);
                // dbg!(area, h0, h1, h2);
                if area < 0.1 || h0 < 0.1 || h1 < 0.1 || h2 < 0.1 { continue; }
                let height = crate::tet::height(&p0, &p1, &p2, &q);
                if height.abs() < 0.1 { continue; }
            }
            if q.z.abs() < 0.2 { continue; }
            // let q = nalgebra::Vector3::<f64>::new(1.5, 0.5, 1.);
            let (val_anly, _) = wdw_integral_of_inverse_distance_cubic(&p0, &p1, &p2, &q);
            let integrand = |u: f64, v: f64| {
                let p = (1. - u - v) * p0 + u * p1 + v * p2;
                let dist = (p - q).norm();
                1. / (dist * dist * dist)
            };
            let val_num = integrate_numerically(&p0, &p1, &p2, integrand, 100);
            // dbg!(val/area, (p0+p1+p2)/3.0);
            // dbg!(val_num, val_anly);
            assert!((val_num - val_anly).abs() < 3.0e-3 * (val_anly + 0.01));
        }
    }

    #[test]
    fn test_wdw_integral_of_inverse_distance_cubic() {
        use crate::tri3::wdw_integral_of_inverse_distance_cubic;
        for _ in 0..1000 {
            let p0 = crate::vec3::sample_unit_cube();
            let p1 = crate::vec3::sample_unit_cube();
            let p2 = crate::vec3::sample_unit_cube();
            let q0 =
                crate::vec3::sample_unit_cube() - nalgebra::Vector3::<f64>::new(0.5, 0.5, 0.5);
            {
                let area = area(&p0, &p1, &p2);
                let h0 = height(&p0, &p1, &p2);
                let h1 = height(&p1, &p2, &p0);
                let h2 = height(&p2, &p0, &p1);
                if area < 0.1 || h0 < 0.1 || h1 < 0.1 || h2 < 0.1 { continue; }
                let h = crate::tet::height(&p0, &p1, &p2, &q0);
                if h.abs() < 0.1 { continue; }
            }
            let (w0, dw) = wdw_integral_of_inverse_distance_cubic(&p0, &p1, &p2, &q0);
            let eps = 1.0e-4_f64;
            //
            let q1x = q0 + nalgebra::Vector3::<f64>::new(eps, 0., 0.);
            let (w1x, _) = wdw_integral_of_inverse_distance_cubic(&p0, &p1, &p2, &q1x);
            // dbg!((w1x-w0)/eps, dw.x);
            assert!(((w1x - w0) / eps - dw.x).abs() < 5.0e-2 * (dw.x.abs() + 0.1));
            //
            let q1y = q0 + nalgebra::Vector3::<f64>::new(0., eps, 0.);
            let (w1y, _) = wdw_integral_of_inverse_distance_cubic(&p0, &p1, &p2, &q1y);
            // dbg!((w1y-w0)/eps, dw.y);
            assert!(((w1y - w0) / eps - dw.y).abs() < 5.0e-2 * (dw.y.abs() + 0.1));
            //
            let q1z = q0 + nalgebra::Vector3::<f64>::new(0., 0., eps);
            let (w1z, _) = wdw_integral_of_inverse_distance_cubic(&p0, &p1, &p2, &q1z);
            // dbg!((w1z-w0)/eps, dw.z);
            assert!(((w1z - w0) / eps - dw.z).abs() < 5.0e-2 * (dw.z.abs() + 0.1));
        }
    }
}
