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
fn cubic_roots_in_range_zero_to_t<T>(
    c0: T,
    c1: T,
    c2: T,
    c3: T,
    t: T,
    epsilon: T) -> Vec<T>
    where T: num_traits::Float + 'static + Copy + std::fmt::Debug,
          i64: AsPrimitive<T>
{
    assert!(t>T::zero());
    let mut result = vec!(T::zero();0);
    let eval_f = |r| ((c3 * r + c2) * r + c1) * r + c0;
    let f0 = c0;
    let ft = eval_f(t);
    if c3.abs() < T::epsilon() { // lesser than cubic function
        if c2.abs() < T::epsilon() { // lesser than quadratic function
            if c1.abs() < T::epsilon() { // lesser than linear function
                // constant function
                if c0.abs() < T::epsilon() { result.push(T::zero()); }
                return result;
            }
            else { // linear function
                if (f0 <= T::zero() && ft >= T::zero()) || (f0 >= T::zero() && ft <= T::zero()) {
                    assert_ne!(f0,ft);
                    result.push(f0 / (f0 - ft));
                }
                return result;
            }
        }
        else { // quadratic function
            if let Some((e0,e1)) = quadratic_root(c0,c1,c2) {
                let (e0, e1) = if e0 < e1 { (e0, e1) } else { (e1, e0) };
                if e0 >= T::zero() && e0 <= t { result.push(e0); }
                if e1 >= T::zero() && e1 <= t { result.push(e1); }
            }
            return result;
        }
    }
    //
    let two = 2.as_();
    let three = 3.as_();
    //
    let newton = |xs: T, xe: T, fs: T, fe: T| {
        if (fs < T::zero() && fe < T::zero()) || (fs > T::zero() && fe > T::zero()) {
            return None;
        }
        assert_ne!(fs,fe);
        let mut r = (fs * xe - fe * xs) / (fs - fe);
        assert!(r >= T::zero() && r <= t);
        for _i in 0..20 {
            let fr = eval_f(r);
            if fr.abs() < epsilon {
                break;
            }
            let dfr = c1 + two * c2 * r + three * c3 * r * r;
            r = r - fr / dfr;
            r = num_traits::clamp(r, T::zero(), t);
        }
        if r < T::zero() || r > t { return None; }
        Some(r)
    };
    // f'(x) = c1 + 2*c2*x + 3*c3*x^2
    if let Some((e0, e1)) = quadratic_root(c1, two * c2, three * c3) {
        assert!(e0<=e1);
        if T::zero() <= e0 { // overlap [0,t] and [-\infty,e0]
            let a1 = t.min(e0);
            let fa1 = eval_f(a1);
            if let Some(r) = newton(T::zero(), a1, f0, fa1) { result.push(r); }
        }
        if e0 <= t && T::zero() <= e1 { // overlap [0,t] and [e0,e1]
            let a0 = T::zero().max(e0);
            let a1 = t.min(e1);
            let fa0 = eval_f(a0);
            let fa1 = eval_f(a1);
            assert!(a0<=a1);
            if let Some(r) = newton(a0, a1, fa0, fa1) { result.push(r); }
        }
        if e1 <= t { // overlap [0,t] and [e1,\infty]
            let a0 = T::zero().max(e1);
            let fa0 = eval_f(a0);
            if let Some(r) = newton(a0, t, fa0, ft) { result.push(r); }
        }
    } else { // monotonic
        if let Some(r) = newton(T::zero(), t, f0, ft) { result.push(r); }
    }
    return result;
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
        let list_time = crate::ccd::cubic_roots_in_range_zero_to_t(
            c0, c1, c2, c3,1.0,eps);
        for t in list_time {
            let fr: f64 = c0 + c1 * t + c2 * t * t + c3 * t * t * t;
            assert!(fr.abs() < eps);
        }
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
    epsilon: T) -> Vec<T>
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
    cubic_roots_in_range_zero_to_t(
        k0, k1, k2, k3, T::one(),
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
    let list_te = coplanar_time(
        f0s, f1s, f2s, vs,
        f0e, f1e, f2e, ve,
        epsilon);
    for te in list_te {
        let ts = T::one() - te;
        let f0 = f0s.scale(ts) + f0e.scale(te);
        let f1 = f1s.scale(ts) + f1e.scale(te);
        let f2 = f2s.scale(ts) + f2e.scale(te);
        let v = vs.scale(ts) + ve.scale(te);
        let coord = crate::tri3::barycentric(&f0, &f1, &f2, &v);
        if coord.x >= T::zero() && coord.y >= T::zero() && coord.z >= T::zero() {
            return Some(te);
        }
    }
    return None;
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
    let list_te = coplanar_time(
        a0s, a1s, b0s, b1s,
        a0e, a1e, b0e, b1e,
        epsilon);
    for te in list_te {
        let ts = T::one() - te;
        let a0 = a0s.scale(ts) + a0e.scale(te);
        let a1 = a1s.scale(ts) + a1e.scale(te);
        let b0 = b0s.scale(ts) + b0e.scale(te);
        let b1 = b1s.scale(ts) + b1e.scale(te);
        let coord
            = crate::edge3::intersection_edge3_when_coplanar(
            &a0, &a1, &b0, &b1);
        let Some(coord) = coord else { continue; }; // coplanar case
        if coord.0 >= T::zero() && coord.1 >= T::zero() &&
            coord.2 >= T::zero() && coord.3 >= T::zero() { return Some(te); }
    }
    return None;
}
