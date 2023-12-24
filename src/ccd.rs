//! methods related to Continuous Collision Detection (CCD)
//!

use num_traits::AsPrimitive;

/// find root of quadratic function
/// f(x) = c0 + c1*x + c2*x^2
fn quadratic_root<T>(
    c0: T,
    c1: T,
    c2: T) -> Option<(T, T)>
    where T: num_traits::Float + 'static + Copy + std::fmt::Debug,
          i64: num_traits::AsPrimitive<T>
{
    let det = c1 * c1 - 4.as_() * c2 * c0;
    if det < T::zero() {
        return None;
    }
    let two = 2.as_();
    let sgnb: T = if c1 < T::zero() { -T::one() } else { T::one() };
    let tmp = -(c1 + sgnb * det.sqrt());
    let x2 = tmp / (two * c2);
    let x1 = (two * c0) / tmp;
    let (x1, x2) = if x1 < x2 { (x1, x2) } else { (x2, x1) };
    assert!(x1 <= x2);
    Some((x1, x2))
}

#[test]
fn test_quadratic_root() {
    use rand::Rng;
    let mut rng = rand::thread_rng();
    for _ in 0..1000 {
        let x0: f64 = 4. * rng.gen::<f64>() - 2.;
        let x1: f64 = 4. * rng.gen::<f64>() - 2.;
        let (x0, x1) = if x0 > x1 { (x1, x0) } else { (x0, x1) };
        let c2: f64 = 4. * rng.gen::<f64>() - 2.;
        let c1: f64 = -(x0 + x1) * c2;
        let c0: f64 = x0 * x1 * c2;
        let res = quadratic_root(c0, c1, c2);
        assert!(res.is_some());
        let (y0, y1) = res.unwrap();
        assert!((x0 - y0).abs() < 1.0e-8);
        assert!((x1 - y1).abs() < 1.0e-8);
        //
        let c0 = c2 * (x0 * x0 + rng.gen::<f64>() + std::f64::EPSILON);
        let c1 = -2. * x0 * c2;
        let res = quadratic_root(c0, c1, c2);
        assert!(res.is_none());
    }
}

/// f(x) = c0 + c1*x + c2*x^2 + c3*x^3
fn cubic_root_in_range_zero_to_one_closest_zero<T>(
    c0: T,
    c1: T,
    c2: T,
    c3: T,
    epsilon: T) -> Option<T>
    where T: num_traits::Float + 'static + Copy + std::fmt::Debug,
          i64: AsPrimitive<T>
{
    let f0 = c0;
    if f0.abs() < T::epsilon() { return Some(T::zero()); }
    let f1 = c0 + c1 + c2 + c3;
    if c3.abs() < T::epsilon() {
        if c2.abs() < T::epsilon() {
            if c1.abs() < T::epsilon() { // constant function
                if c0.abs() < T::epsilon() {
                    return Some(T::zero());
                }
                else {
                    return None;
                }
            }
            else { // linear function
                return if (f0 < T::zero() && f1 > T::zero()) || (f0 > T::zero() && f1 < T::zero()) {
                    Some(f0 / (f0 - f1))
                } else {
                    None
                }
            }
        }
        else { // quadratic function
            panic!();
        }
    }
    //
    let two = 2.as_();
    let three = 3.as_();
    //
    let eval_f = |r| ((c3 * r + c2) * r + c1) * r + c0;
    let newton = |x0: T, x1: T, f0: T, f1: T| {
        if (f0 < T::zero() && f1 < T::zero()) || (f0 > T::zero() && f1 > T::zero()) {
            return None;
        }
        let mut r = (f0 * x1 - f1 * x0) / (f0 - f1);
        assert!(r >= T::zero() && r <= T::one());
        for _i in 0..20 {
            let fr = eval_f(r);
            if fr.abs() < epsilon {
                break;
            }
            let dfr = c1 + two * c2 * r + three * c3 * r * r;
            r = r - fr / dfr;
        }
        if r < T::zero() || r > T::one() { return None; }
        Some(r)
    };
    // f'(x) = c1 + 2*c2*x + 3*c3*x^2
    if let Some((e0, e1)) = quadratic_root(c1, two * c2, three * c3) {
        let fe0 = eval_f(e0);
        let fe1 = eval_f(e1);
        if T::zero() < e0 && e0 < T::one() { // overlap [0.1] and [-\infty,e0]
            if let Some(r) = newton(T::zero(), e0, f0, fe0) { return Some(r); }
        }
        if e0 < T::one() && T::zero() < e1 { // overlap [0,1] and [e0,e1]
            let a0 = T::zero().max(e0);
            let a1 = T::one().min(e1);
            let fa0 = eval_f(a0);
            let fa1 = eval_f(a1);
            if let Some(r) = newton(a0, a1, fa0, fa1) { return Some(r); }
        }
        if T::zero() < e1 && e1 < T::one() { // overlap [0.1] and [-\infty,e0]
            if let Some(r) = newton(e1, T::one(), fe1, f1) { return Some(r); }
        }
    } else { // monotonic
        if let Some(r) = newton(T::zero(), T::one(), f0, f1) { return Some(r); }
    }
    None
}

#[test]
fn test_cubic_root() {
    use rand::Rng;
    let mut rng = rand::thread_rng();
    let eps = 1.0e-8;
    for _ in 0..10000 {
        let c0: f64 = 4. * rng.gen::<f64>() - 2.;
        let c1: f64 = 4. * rng.gen::<f64>() - 2.;
        let c2: f64 = 4. * rng.gen::<f64>() - 2.;
        let c3: f64 = 4. * rng.gen::<f64>() - 2.;
        if let Some(r) = crate::ccd::cubic_root_in_range_zero_to_one_closest_zero(c0, c1, c2, c3, eps) {
            let fr: f64 = c0 + c1 * r + c2 * r * r + c3 * r * r * r;
            assert!(fr.abs() < eps);
        } else {}
    }
}

/// compute time where four points gets co-planar
fn coplanar_time<T>(
    s0: &nalgebra::Vector3::<T>,
    s1: &nalgebra::Vector3::<T>,
    s2: &nalgebra::Vector3::<T>,
    s3: &nalgebra::Vector3::<T>,
    e0: &nalgebra::Vector3::<T>,
    e1: &nalgebra::Vector3::<T>,
    e2: &nalgebra::Vector3::<T>,
    e3: &nalgebra::Vector3::<T>,
    epsilon: T) -> Option<T>
    where T: nalgebra::RealField + std::marker::Copy + num_traits::Float,
          i64: AsPrimitive<T>
{
    let x1 = s1 - s0;
    let x2 = s2 - s0;
    let x3 = s3 - s0;
    let v1 = e1 - e0 - x1;
    let v2 = e2 - e0 - x2;
    let v3 = e3 - e0 - x3;
    // compute coefficient for cubic function
    use crate::vec3::scalar_triple_product;
    let k0 = scalar_triple_product(&x3, &x1, &x2);
    let k1 = scalar_triple_product(&v3, &x1, &x2)
        + scalar_triple_product(&x3, &v1, &x2)
        + scalar_triple_product(&x3, &x1, &v2);
    let k2 = scalar_triple_product(&v3, &v1, &x2)
        + scalar_triple_product(&v3, &x1, &v2)
        + scalar_triple_product(&x3, &v1, &v2);
    let k3 = scalar_triple_product(&v3, &v1, &v2);
    // cubic function is f(x) = k0 + k1*x + k2*x^2 + k3*x^3
    cubic_root_in_range_zero_to_one_closest_zero(
        k0, k1, k2, k3,
        epsilon)
}

pub fn intersecting_time_fv<T>(
    f0s: &nalgebra::Vector3::<T>,
    f1s: &nalgebra::Vector3::<T>,
    f2s: &nalgebra::Vector3::<T>,
    vs: &nalgebra::Vector3::<T>,
    f0e: &nalgebra::Vector3::<T>,
    f1e: &nalgebra::Vector3::<T>,
    f2e: &nalgebra::Vector3::<T>,
    ve: &nalgebra::Vector3::<T>,
    epsilon: T) -> Option<T>
    where T: nalgebra::RealField + Copy + num_traits::Float,
          i64: AsPrimitive<T>,
          f64: AsPrimitive<T>
{
    let te = coplanar_time(f0s, f1s, f2s, vs, f0e, f1e, f2e, ve, epsilon);
    let Some(te) = te else { return None; };
    let ts = T::one() - te;
    let f0 = f0s.scale(ts) + f0e.scale(te);
    let f1 = f1s.scale(ts) + f1e.scale(te);
    let f2 = f2s.scale(ts) + f2e.scale(te);
    let v = vs.scale(ts) + ve.scale(te);
    let coord = crate::tri3::barycentric(&f0, &f1, &f2, &v);
    if coord.x < T::zero() || coord.y < T::zero() || coord.z < T::zero() { return None; }
    Some(te)
}

pub fn intersecting_time_ee<T>(
    a0s: &nalgebra::Vector3::<T>,
    a1s: &nalgebra::Vector3::<T>,
    b0s: &nalgebra::Vector3::<T>,
    b1s: &nalgebra::Vector3::<T>,
    a0e: &nalgebra::Vector3::<T>,
    a1e: &nalgebra::Vector3::<T>,
    b0e: &nalgebra::Vector3::<T>,
    b1e: &nalgebra::Vector3::<T>,
    epsilon: T) -> Option<T>
    where T: nalgebra::RealField + Copy + num_traits::Float,
          i64: AsPrimitive<T>,
          f64: AsPrimitive<T>
{
    let te = coplanar_time(
        a0s, a1s, b0s, b1s,
        a0e, a1e, b0e, b1e,
        epsilon);
    let Some(te) = te else { return None; };
    let ts = T::one() - te;
    let a0 = a0s.scale(ts) + a0e.scale(te);
    let a1 = a1s.scale(ts) + a1e.scale(te);
    let b0 = b0s.scale(ts) + b0e.scale(te);
    let b1 = b1s.scale(ts) + b1e.scale(te);
    let coord = crate::edge3::intersection_edge3(&a0, &a1, &b0, &b1);
    let Some(coord) = coord else { return None; }; // coplanar case
    if coord.0 < T::zero() || coord.1 < T::zero() || coord.2 < T::zero() || coord.3 < T::zero() { return None; }
    /*
    {
        let n = (a1-a0).cross(&(b0-a0));
        dbg!(a0.dot(&n), a1.dot(&n), b0.dot(&n), b1.dot(&n));
        let a2 = a0 + n;
        let rq1 = crate::tet::volume(&a0,&a1,&a2, &b0);
        let rq0 = crate::tet::volume(&a0,&a1,&a2, &b1);
        let rp1 = crate::tet::volume(&a0,&a1,&a2, &b0);
        let rp0 = crate::tet::volume(&a0,&a1,&a2, &b1);
        dbg!(rq1,rq0,rp1,rp0);
    }
    dbg!((a0.scale(coord.0)+a1.scale(coord.1), b0.scale(coord.2)+b1.scale(coord.3)));
    dbg!(&coord,te);
     */
    Some(te)
}
