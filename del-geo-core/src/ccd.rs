//! methods related to Continuous Collision Detection (CCD)

use num_traits::AsPrimitive;

/// find root of quadratic function
/// f(x) = c0 + c1*x + c2*x^2
fn quadratic_root<T>(c0: T, c1: T, c2: T) -> Option<(T, T)>
where
    T: num_traits::Float + 'static + Copy + std::fmt::Debug,
    i64: AsPrimitive<T>,
{
    assert_ne!(c2, T::zero());
    let det = c1 * c1 - 4.as_() * c2 * c0;
    if det < T::zero() {
        return None;
    }
    let two = 2.as_();
    let sgnb: T = if c1 < T::zero() { -T::one() } else { T::one() };
    let tmp = -(c1 + sgnb * det.sqrt());
    let x2 = tmp / (two * c2);
    let x1 = if tmp != T::zero() {
        (two * c0) / tmp
    } else {
        // c1 ==0, det == 0
        x2
    };
    let (x1, x2) = if x1 < x2 { (x1, x2) } else { (x2, x1) };
    if x1 <= x2 {
    } else {
        dbg!(c0, c1, c2, det, tmp, x1, x2);
    }
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
        let c0 = c2 * (x0 * x0 + rng.gen::<f64>() + f64::EPSILON);
        let c1 = -2. * x0 * c2;
        let res = quadratic_root(c0, c1, c2);
        assert!(res.is_none());
    }
}

/// f(x) = c0 + c1*x + c2*x^2 + c3*x^3
pub fn cubic_roots_in_range_zero_to_t<T>(c0: T, c1: T, c2: T, c3: T, t: T, epsilon: T) -> Vec<T>
where
    T: num_traits::Float + 'static + Copy + std::fmt::Debug + std::fmt::Display,
    i64: AsPrimitive<T>,
{
    assert!(t > T::zero());
    let mut result = vec![];
    let eval_f = |r| ((c3 * r + c2) * r + c1) * r + c0;
    let f0 = c0;
    let ft = eval_f(t);
    if c3.abs() < T::epsilon() {
        // lesser than cubic function
        return if c2.abs() < T::epsilon() {
            // lesser than quadratic function
            if c1.abs() < T::epsilon() {
                // lesser than linear function
                // constant function
                if c0.abs() < T::epsilon() {
                    result.push(T::zero());
                }
                result
            } else {
                // linear function
                if (f0 <= T::zero() && ft >= T::zero()) || (f0 >= T::zero() && ft <= T::zero()) {
                    assert_ne!(f0, ft);
                    result.push(f0 / (f0 - ft));
                }
                result
            }
        } else {
            // quadratic function
            if let Some((e0, e1)) = quadratic_root(c0, c1, c2) {
                let (e0, e1) = if e0 < e1 { (e0, e1) } else { (e1, e0) };
                if e0 >= T::zero() && e0 <= t {
                    result.push(e0);
                }
                if e1 >= T::zero() && e1 <= t {
                    result.push(e1);
                }
            }
            result
        };
    }
    //
    let two = 2.as_();
    let three = 3.as_();
    //
    let newton = |xs: T, xe: T, fs: T, fe: T| {
        if (fs < T::zero() && fe < T::zero()) || (fs > T::zero() && fe > T::zero()) {
            return None;
        }
        if xs == xe {
            return if fs == T::zero() { Some(xs) } else { None };
        }
        assert_ne!(fs, fe, "hoge {} {} {} {} {} {}", xs, xe, c0, c1, c2, c3);
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
        if r < T::zero() || r > t {
            return None;
        }
        Some(r)
    };
    // f'(x) = c1 + 2*c2*x + 3*c3*x^2
    if let Some((e0, e1)) = quadratic_root(c1, two * c2, three * c3) {
        assert!(e0 <= e1);
        if T::zero() <= e0 {
            // overlap [0,t] and [-\infty,e0]
            let a1 = t.min(e0);
            let fa1 = eval_f(a1);
            if let Some(r) = newton(T::zero(), a1, f0, fa1) {
                result.push(r);
            }
        }
        if e0 <= t && T::zero() <= e1 {
            // overlap [0,t] and [e0,e1]
            let a0 = T::zero().max(e0);
            let a1 = t.min(e1);
            let fa0 = eval_f(a0);
            let fa1 = eval_f(a1);
            assert!(a0 <= a1);
            if let Some(r) = newton(a0, a1, fa0, fa1) {
                result.push(r);
            }
        }
        if e1 <= t {
            // overlap [0,t] and [e1,\infty]
            let a0 = T::zero().max(e1);
            let fa0 = eval_f(a0);
            if let Some(r) = newton(a0, t, fa0, ft) {
                result.push(r);
            }
        }
    } else {
        // monotonic
        if let Some(r) = newton(T::zero(), t, f0, ft) {
            result.push(r);
        }
    }
    result
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
        let list_time = cubic_roots_in_range_zero_to_t(c0, c1, c2, c3, 1.0, eps);
        for t in list_time {
            let fr: f64 = c0 + c1 * t + c2 * t * t + c3 * t * t * t;
            assert!(fr.abs() < eps);
        }
    }
}
