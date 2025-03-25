//! methods for 3D triangle (parameterized by three corner points)

use num_traits::AsPrimitive;

/// height of triangle vertex `p2` against the edge connecting `p0` and `p1`
pub fn height<T>(
    p0: &nalgebra::Vector3<T>,
    p1: &nalgebra::Vector3<T>,
    p2: &nalgebra::Vector3<T>,
) -> T
where
    T: nalgebra::RealField + 'static + Copy + num_traits::Float,
    f64: AsPrimitive<T>,
{
    let a = area(p2, p0, p1);
    a * 2.0.as_() / (p0 - p1).norm()
}

/// normal of a triangle
pub fn normal<T>(
    p0: &nalgebra::Vector3<T>,
    p1: &nalgebra::Vector3<T>,
    p2: &nalgebra::Vector3<T>,
) -> nalgebra::Vector3<T>
where
    T: nalgebra::RealField,
{
    (p1 - p0).cross(&(p2 - p0))
}

pub fn dw_normal<T>(
    p0: &nalgebra::Vector3<T>,
    p1: &nalgebra::Vector3<T>,
    p2: &nalgebra::Vector3<T>,
) -> [nalgebra::Matrix3<T>; 3]
where
    T: nalgebra::RealField + Copy,
{
    [
        crate::mat3::skew(&(p2 - p1)),
        crate::mat3::skew(&(p0 - p2)),
        crate::mat3::skew(&(p1 - p0)),
    ]
}

#[test]
fn test_normal() {
    type V3 = nalgebra::Vector3<f64>;
    {
        let p0 = [
            V3::new(0.0, 0.0, 0.0),
            V3::new(1.0, 0.0, 0.0),
            V3::new(0.0, 1.0, 0.0),
        ];
        let pn = normal(&p0[0], &p0[1], &p0[2]);
        assert!((pn - V3::new(0., 0., 1.)).norm() < 1.0e-10);
    }
    {
        let p0 = V3::new(0.1, 0.5, 0.4);
        let p1 = V3::new(1.2, 0.3, 0.2);
        let p2 = V3::new(0.3, 1.4, 0.1);
        let pn = normal(&p0, &p1, &p2);
        assert!((p0 - p1).dot(&pn) < 1.0e-10);
        assert!((p1 - p2).dot(&pn) < 1.0e-10);
        assert!((p2 - p0).dot(&pn) < 1.0e-10);
    }
    {
        let p0 = [
            V3::new(0.1, 0.4, 0.2),
            V3::new(1.2, 0.3, 0.7),
            V3::new(0.3, 1.5, 0.3),
        ];
        let cc0 = normal(&p0[0], &p0[1], &p0[2]);
        let dcc = dw_normal(&p0[0], &p0[1], &p0[2]);
        let eps = 1.0e-5;
        for i_node in 0..3 {
            for i_dim in 0..3 {
                let p1 = {
                    let mut p1 = p0;
                    p1[i_node][i_dim] += eps;
                    p1
                };
                let cc1 = normal(&p1[0], &p1[1], &p1[2]);
                let dcc_num = (cc1 - cc0) / eps;
                let mut b = V3::zeros();
                b[i_dim] = 1.0;
                let dcc_ana = dcc[i_node] * b;
                let diff = (dcc_num - dcc_ana).norm();
                assert!(diff < 1.0e-4);
                println!(
                    "normal {}, {} --> {:?}, {:?}, {:?}",
                    i_node, i_dim, dcc_num, dcc_ana, diff
                );
            }
        }
    }
}

/// unit normal of a triangle
pub fn unit_normal<T>(
    p0: &nalgebra::Vector3<T>,
    p1: &nalgebra::Vector3<T>,
    p2: &nalgebra::Vector3<T>,
) -> nalgebra::Vector3<T>
where
    T: nalgebra::RealField,
{
    let n = (p1 - p0).cross(&(p2 - p0));
    n.normalize()
}

pub fn area<T>(p0: &nalgebra::Vector3<T>, p1: &nalgebra::Vector3<T>, p2: &nalgebra::Vector3<T>) -> T
where
    T: nalgebra::RealField + 'static + Copy,
    f64: AsPrimitive<T>,
{
    (p1 - p0).cross(&(p2 - p0)).norm() * 0.5_f64.as_()
}

pub fn nearest_to_point3<T>(
    q0: &nalgebra::Vector3<T>,
    q1: &nalgebra::Vector3<T>,
    q2: &nalgebra::Vector3<T>,
    ps: &nalgebra::Vector3<T>,
) -> (nalgebra::Vector3<T>, T, T)
where
    T: nalgebra::RealField + Copy + 'static,
    f64: AsPrimitive<T>,
{
    let n012 = unit_normal(q0, q1, q2);
    let pe = ps + n012;
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
            let nearp = q0.scale(r0) + q1.scale(r1) + q2.scale(r2);
            return (nearp, r0, r1);
        }
    }
    let r12 = crate::edge3::nearest_to_point3(q1, q2, ps);
    let r20 = crate::edge3::nearest_to_point3(q2, q0, ps);
    let r01 = crate::edge3::nearest_to_point3(q0, q1, ps);
    let d12 = r12.0;
    let d20 = r20.0;
    let d01 = r01.0;
    let r12 = q1 + (q2 - q1).scale(r12.1);
    let r20 = q2 + (q0 - q2).scale(r20.1);
    let r01 = q0 + (q1 - q0).scale(r01.1);
    if d12 < d20 {
        if d12 < d01 {
            // 12 is the smallest
            let r0 = T::zero();
            let r1 = (r12 - q2).norm() / (q1 - q2).norm();
            return (r12, r0, r1);
        }
    } else if d20 < d01 {
        // d20 is the smallest
        let r0 = (r20 - q2).norm() / (q0 - q2).norm();
        let r1 = T::zero();
        return (r20, r0, r1);
    }
    let r0 = (r01 - q1).norm() / (q0 - q1).norm();
    let r1 = T::one() - r0;
    (r01, r0, r1)
}

#[test]
fn test_triangle_point_nearest() {
    use rand::SeedableRng;
    let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    let eps = 1.0 - 5f64;
    for _ in 0..10000 {
        let q0 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let q1 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let q2 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let ps = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
        let (qs, r0, r1) = nearest_to_point3(&q0, &q1, &q2, &ps);
        let dist0 = (qs - ps).norm();
        use del_geo_core::tri3::clamp_barycentric_coords;
        let (r0a, r1a, r2a) =
            clamp_barycentric_coords(r0 + eps, r1 + eps, 1. - r0 - eps - r1 - eps);
        let qa = q0.scale(r0a) + q1.scale(r1a) + q2.scale(r2a);
        assert!((qa - ps).norm() >= dist0);
        let (r0b, r1b, r2b) =
            clamp_barycentric_coords(r0 + eps, r1 - eps, 1. - r0 - eps - r1 + eps);
        let qb = q0.scale(r0b) + q1.scale(r1b) + q2.scale(r2b);
        assert!((qb - ps).norm() >= dist0);
        let (r0c, r1c, r2c) =
            clamp_barycentric_coords(r0 - eps, r1 + eps, 1. - r0 + eps - r1 - eps);
        let qc = q0.scale(r0c) + q1.scale(r1c) + q2.scale(r2c);
        assert!((qc - ps).norm() >= dist0);
        let (r0d, r1d, r2d) =
            clamp_barycentric_coords(r0 - eps, r1 - eps, 1. - r0 + eps - r1 + eps);
        let qd = q0.scale(r0d) + q1.scale(r1d) + q2.scale(r2d);
        assert!((qd - ps).norm() >= dist0);
    }
}

fn wdw_inverse_distance_cubic_integrated_over_wedge(
    x: nalgebra::Vector3<f64>,
    b: f64,
) -> (f64, nalgebra::Vector3<f64>) {
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
    (
        w,
        nalgebra::Vector3::<f64>::new(dwdx, dwdy, dwdz * signz * 2.),
    )
}

pub fn wdw_integral_of_inverse_distance_cubic(
    p0: &nalgebra::Vector3<f64>,
    p1: &nalgebra::Vector3<f64>,
    p2: &nalgebra::Vector3<f64>,
    q: &nalgebra::Vector3<f64>,
) -> (f64, nalgebra::Vector3<f64>) {
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
    let dw = -vz.scale(signz * core::f64::consts::PI * 2_f64 / (z * z)) - dw0dq - dw1dq - dw2dq;
    (w, dw)
}

pub fn numerical_integration<F>(
    p0: &nalgebra::Vector3<f64>,
    p1: &nalgebra::Vector3<f64>,
    p2: &nalgebra::Vector3<f64>,
    integrand: F,
    n: usize,
) -> f64
where
    F: Fn(f64, f64) -> f64,
{
    let area = area(p0, p1, p2);
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
                _ => {
                    panic!();
                }
            };
            let (u, v) = (u as f64 / (n * 3) as f64, v as f64 / (n * 3) as f64);
            let dist = integrand(u, v);
            val_num += jacobian * dist;
        }
    }
    val_num
}

pub fn is_intersection_tri3_sat(
    p0: &nalgebra::Vector3<f32>,
    p1: &nalgebra::Vector3<f32>,
    p2: &nalgebra::Vector3<f32>,
    q0: &nalgebra::Vector3<f32>,
    q1: &nalgebra::Vector3<f32>,
    q2: &nalgebra::Vector3<f32>,
) -> bool {
    let ps = nalgebra::Matrix3::<f32>::from_columns(&[*p0, *p1, *p2]);
    let qs = nalgebra::Matrix3::<f32>::from_columns(&[*q0, *q1, *q2]);
    let sep = |dir: nalgebra::Vector3<f32>| {
        let prj0 = dir.transpose() * ps;
        let prj1 = dir.transpose() * qs;
        prj0.min() > prj1.max() || prj0.max() < prj1.min()
    };
    if sep((p1 - p0).cross(&(p2 - p0))) {
        return false;
    }
    if sep((q1 - q0).cross(&(q2 - q0))) {
        return false;
    }
    if sep((p0 - p1).cross(&(q0 - q1))) {
        return false;
    }
    if sep((p0 - p1).cross(&(q1 - q2))) {
        return false;
    }
    if sep((p0 - p1).cross(&(q2 - q0))) {
        return false;
    }
    if sep((p1 - p2).cross(&(q0 - q1))) {
        return false;
    }
    if sep((p1 - p2).cross(&(q1 - q2))) {
        return false;
    }
    if sep((p1 - p2).cross(&(q2 - q0))) {
        return false;
    }
    if sep((p2 - p0).cross(&(q0 - q1))) {
        return false;
    }
    if sep((p2 - p0).cross(&(q1 - q2))) {
        return false;
    }
    if sep((p2 - p0).cross(&(q2 - q0))) {
        return false;
    }
    true
}

/// if the triangle share a point, set the point as `p0` and `q0`
pub fn is_intersection_tri3<T>(
    p0: &nalgebra::Vector3<T>,
    p1: &nalgebra::Vector3<T>,
    p2: &nalgebra::Vector3<T>,
    q0: &nalgebra::Vector3<T>,
    q1: &nalgebra::Vector3<T>,
    q2: &nalgebra::Vector3<T>,
) -> Option<(nalgebra::Vector3<T>, nalgebra::Vector3<T>)>
where
    T: nalgebra::RealField + Copy,
{
    let sec = |p0: &nalgebra::Vector3<T>, p1: &nalgebra::Vector3<T>, dp0: T, dp1: T| {
        let r1 = dp0 / (dp0 - dp1);
        let r0 = T::one() - r1;
        p0.scale(r0) + p1.scale(r1)
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
        let dp0 = (p0 - q0).dot(&nq);
        let dp1 = (p1 - q0).dot(&nq);
        let dp2 = (p2 - q0).dot(&nq);
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
        let mut ap = Vec::<nalgebra::Vector3<T>>::with_capacity(2);
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
        let dq0 = (q0 - p0).dot(&np); // projection of q0
        let dq1 = (q1 - p0).dot(&np); // projection of q1
        let dq2 = (q2 - p0).dot(&np); // projection of q2
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
        let mut aq = Vec::<nalgebra::Vector3<T>>::with_capacity(2);
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

#[test]
fn test_triangle_intersection() {
    use rand::SeedableRng;
    let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    for _ in 0..10000 {
        let p0 = crate::vec3::sample_unit_cube(&mut reng);
        let p1 = crate::vec3::sample_unit_cube(&mut reng);
        let p2 = crate::vec3::sample_unit_cube(&mut reng);
        if area(&p0, &p1, &p2) < 1.0e-2 {
            continue;
        }
        let q0 = crate::vec3::sample_unit_cube(&mut reng);
        let q1 = crate::vec3::sample_unit_cube(&mut reng);
        let q2 = crate::vec3::sample_unit_cube(&mut reng);
        if area(&q0, &q1, &q2) < 1.0e-2 {
            continue;
        }
        let res = is_intersection_tri3(&p0, &p1, &p2, &q0, &q1, &q2);
        {
            let res_sat = is_intersection_tri3_sat(&p0, &p1, &p2, &q0, &q1, &q2);
            assert_eq!(res_sat, res.is_some());
        }
        if let Some((r0, r1)) = res {
            assert!(crate::tet::height(&p0, &p1, &p2, &r0).abs() < 2.0e-6);
            assert!(crate::tet::height(&p0, &p1, &p2, &r1).abs() < 2.0e-6);
            assert!(crate::tet::height(&q0, &q1, &q2, &r0).abs() < 2.0e-6);
            assert!(crate::tet::height(&q0, &q1, &q2, &r1).abs() < 2.0e-6);
            let bp0 = barycentric(&p0, &p1, &p2, &r0);
            let bp1 = barycentric(&p0, &p1, &p2, &r1);
            let bq0 = barycentric(&q0, &q1, &q2, &r0);
            let bq1 = barycentric(&q0, &q1, &q2, &r1);
            assert!(bp0.min() > -2.0e-5);
            assert!(bp1.min() > -2.0e-5);
            assert!(bq0.min() > -2.0e-5);
            assert!(bq1.min() > -2.0e-5);
            assert!(bp0.min().min(bq0.min()) < 3.0e-5);
            assert!(bp1.min().min(bq1.min()) < 3.0e-5);
        }
    }
}

/// q must be coplanar to the triangle (p0,p1,p2)
pub fn barycentric<T>(
    p0: &nalgebra::Vector3<T>,
    p1: &nalgebra::Vector3<T>,
    p2: &nalgebra::Vector3<T>,
    q: &nalgebra::Vector3<T>,
) -> nalgebra::Vector3<T>
where
    T: nalgebra::RealField + Copy,
    f64: AsPrimitive<T>,
{
    let n = (p1 - p0).cross(&(p2 - p0));
    let s = q + n;
    let r0 = crate::tet::volume(p1, p2, q, &s);
    let r1 = crate::tet::volume(p2, p0, q, &s);
    let r2 = crate::tet::volume(p0, p1, q, &s);
    let tmp = T::one() / (r0 + r1 + r2);
    nalgebra::Vector3::<T>::new(r0 * tmp, r1 * tmp, r2 * tmp)
}

// -----------------------------------

#[cfg(test)]
mod tests {
    use rand::Rng;

    use crate::tri3::{
        area, barycentric, numerical_integration, wdw_integral_of_inverse_distance_cubic,
    };

    #[test]
    fn test_w_inverse_distance_cubic_integrated_over_wedge() {
        use rand::SeedableRng;
        let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
        for _ in 0..1000 {
            let x = crate::vec3::sample_unit_cube::<_, f64>(&mut reng)
                - nalgebra::Vector3::new(0.5, 0.5, 0.5);
            if x.z.abs() < 0.1 {
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
                    let pos = nalgebra::Vector3::<f64>::new(r * t.cos(), r * t.sin(), 0.);
                    let area = rad / n_r as f64 * r * b / n_t as f64;
                    val_nmr += area * (pos - x).norm().powi(-3);
                }
            }
            use crate::tri3::wdw_inverse_distance_cubic_integrated_over_wedge;
            let (v_anl, _) = wdw_inverse_distance_cubic_integrated_over_wedge(x, b);
            assert!((v_anl - val_nmr).abs() < v_anl * 0.1);
        }
    }

    #[test]
    fn test_dw_inverse_distance_cubic_integrated_over_wedge() {
        use rand::SeedableRng;
        let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
        use crate::tri3::wdw_inverse_distance_cubic_integrated_over_wedge;
        for _ in 0..100000 {
            let x0 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng)
                - nalgebra::Vector3::new(0.5, 0.5, 0.5);
            if x0.z.abs() < 0.2 {
                continue;
            }
            let b = core::f64::consts::PI * (rand::rng().random::<f64>() * 0.8 + 0.1); // 90 degree
            let eps = 1.0e-4_f64;
            let (w0, dw) = wdw_inverse_distance_cubic_integrated_over_wedge(x0, b);
            let x1x = x0 + nalgebra::Vector3::<f64>::new(eps, 0., 0.);
            let (w1x, _) = wdw_inverse_distance_cubic_integrated_over_wedge(x1x, b);
            let x1y = x0 + nalgebra::Vector3::<f64>::new(0., eps, 0.);
            let (w1y, _) = wdw_inverse_distance_cubic_integrated_over_wedge(x1y, b);
            let x1z = x0 + nalgebra::Vector3::<f64>::new(0., 0., eps);
            let (w1z, _) = wdw_inverse_distance_cubic_integrated_over_wedge(x1z, b);
            assert!(((w1x - w0) / eps - dw.x).abs() < 2.0e-2 * (dw.x.abs() + 1.0e-1));
            assert!(((w1y - w0) / eps - dw.y).abs() < 3.0e-2 * (dw.y.abs() + 5.0e-1));
            assert!(((w1z - w0) / eps - dw.z).abs() < 2.0e-2 * (dw.z.abs() + 1.0e-1));
        }
    }

    #[test]
    fn test_w_inverse_distance_cubic_integrated() {
        use rand::SeedableRng;
        let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
        for _ in 0..1000 {
            let p0 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
            let p1 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
            let p2 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
            let q = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
            {
                let area = area(&p0, &p1, &p2);
                let h0 = crate::tri3::height(&p1, &p2, &p0);
                let h1 = crate::tri3::height(&p2, &p0, &p1);
                let h2 = crate::tri3::height(&p0, &p1, &p2);
                // dbg!(area, h0, h1, h2);
                if area < 0.1 || h0 < 0.1 || h1 < 0.1 || h2 < 0.1 {
                    continue;
                }
                let height = crate::tet::height(&p0, &p1, &p2, &q);
                if height.abs() < 0.1 {
                    continue;
                }
            }
            if q.z.abs() < 0.2 {
                continue;
            }
            let (val_anly, _) = wdw_integral_of_inverse_distance_cubic(&p0, &p1, &p2, &q);
            let integrand = |u: f64, v: f64| {
                let p = (1. - u - v) * p0 + u * p1 + v * p2;
                let dist = (p - q).norm();
                1. / (dist * dist * dist)
            };
            let val_num = numerical_integration(&p0, &p1, &p2, integrand, 100);
            assert!((val_num - val_anly).abs() < 3.0e-3 * (val_anly + 0.01));
        }
    }

    #[test]
    fn test_wdw_integral_of_inverse_distance_cubic() {
        use rand::SeedableRng;
        let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
        use crate::tri3::wdw_integral_of_inverse_distance_cubic;
        for _ in 0..10000 {
            let p0 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
            let p1 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
            let p2 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
            let q0 = crate::vec3::sample_unit_cube::<_, f64>(&mut reng);
            {
                let area = area(&p0, &p1, &p2);
                let h0 = crate::tri3::height(&p1, &p2, &p0);
                let h1 = crate::tri3::height(&p2, &p0, &p1);
                let h2 = crate::tri3::height(&p0, &p1, &p2);
                if area < 0.1 || h0 < 0.1 || h1 < 0.1 || h2 < 0.1 {
                    continue;
                }
                let h = crate::tet::height(&p0, &p1, &p2, &q0);
                if h.abs() < 0.1 {
                    continue;
                }
                let b = barycentric(&p0, &p1, &p2, &q0);
                if b.abs().min() < 1.0e-2 {
                    continue;
                }
            }
            // dbg!(p0-q0,p1-q0,p2-q0);
            let (w0, dw) = wdw_integral_of_inverse_distance_cubic(&p0, &p1, &p2, &q0);
            let eps = 1.0e-4_f64;
            //
            let q1x = q0 + nalgebra::Vector3::<f64>::new(eps, 0., 0.);
            let (w1x, _) = wdw_integral_of_inverse_distance_cubic(&p0, &p1, &p2, &q1x);
            // dbg!((w1x-w0)/eps, dw.x);
            assert!(((w1x - w0) / eps - dw.x).abs() < 7.0e-2 * (dw.x.abs() + 1.));
            //
            let q1y = q0 + nalgebra::Vector3::<f64>::new(0., eps, 0.);
            let (w1y, _) = wdw_integral_of_inverse_distance_cubic(&p0, &p1, &p2, &q1y);
            // dbg!((w1y-w0)/eps, dw.y);
            assert!(((w1y - w0) / eps - dw.y).abs() < 7.0e-2 * (dw.y.abs() + 1.));
            //
            let q1z = q0 + nalgebra::Vector3::<f64>::new(0., 0., eps);
            let (w1z, _) = wdw_integral_of_inverse_distance_cubic(&p0, &p1, &p2, &q1z);
            // dbg!((w1z-w0)/eps, dw.z);
            assert!(((w1z - w0) / eps - dw.z).abs() < 7.0e-2 * (dw.z.abs() + 1.));
        }
    }
}
