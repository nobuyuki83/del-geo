//! methods for 3D Continuous Collision Detection (CCD)

use num_traits::{AsPrimitive, Float, FromPrimitive};

pub struct FourPoints<'a, T> {
    pub p0: &'a [T; 3],
    pub p1: &'a [T; 3],
    pub p2: &'a [T; 3],
    pub p3: &'a [T; 3],
}

/// compute time when four points gets co-planar
fn coplanar_time<T>(s: FourPoints<T>, e: FourPoints<T>, epsilon: T) -> Vec<T>
where
    T: Copy + num_traits::Float + 'static + std::fmt::Debug + std::fmt::Display,
    i64: AsPrimitive<T>,
{
    use crate::vec3::Vec3;
    let x1 = s.p1.sub(s.p0);
    let x2 = s.p2.sub(s.p0);
    let x3 = s.p3.sub(s.p0);
    let v1 = e.p1.sub(e.p0).sub(&x1);
    let v2 = e.p2.sub(e.p0).sub(&x2);
    let v3 = e.p3.sub(e.p0).sub(&x3);
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
    crate::polynomial_root::cubic_roots_in_range_zero_to_t(k0, k1, k2, k3, T::one(), epsilon)
}

pub struct FaceVertex<'a, T> {
    pub f0: &'a [T; 3],
    pub f1: &'a [T; 3],
    pub f2: &'a [T; 3],
    pub v: &'a [T; 3],
}

pub fn intersecting_time_fv<T>(s: FaceVertex<T>, e: FaceVertex<T>, epsilon: T) -> Option<T>
where
    T: Copy + num_traits::Float + 'static + std::fmt::Debug + std::fmt::Display,
    i64: AsPrimitive<T>,
    f64: AsPrimitive<T>,
{
    use crate::vec3::Vec3;
    let list_te = coplanar_time(
        FourPoints {
            p0: s.f0,
            p1: s.f1,
            p2: s.f2,
            p3: s.v,
        },
        FourPoints {
            p0: e.f0,
            p1: e.f1,
            p2: e.f2,
            p3: e.v,
        },
        epsilon,
    );
    for te in list_te {
        let ts = T::one() - te;
        let f0 = s.f0.scale(ts).add(&e.f0.scale(te));
        let f1 = s.f1.scale(ts).add(&e.f1.scale(te));
        let f2 = s.f2.scale(ts).add(&e.f2.scale(te));
        let v = s.v.scale(ts).add(&e.v.scale(te));
        // println!("{:?}, {:?}, {:?}, {:?}", f0, f1, f2, v);
        // println!("stt, volume {}", crate::tet::volume(&s.f0, &s.f1, &s.f2, &s.v));
        // println!("end, volume {}", crate::tet::volume(&e.f0, &e.f1, &e.f2, &e.v));
        // println!("time {}, volume {}", te, crate::tet::volume(&f0, &f1, &f2, &v));
        let coord = crate::tri3::to_barycentric_coords(&f0, &f1, &f2, &v);
        // println!("coord {:?}", coord);
        if coord[0] >= T::zero() && coord[1] >= T::zero() && coord[2] >= T::zero() {
            return Some(te);
        }
    }
    None
}

pub struct EdgeEdge<'a, T> {
    pub a0: &'a [T; 3],
    pub a1: &'a [T; 3],
    pub b0: &'a [T; 3],
    pub b1: &'a [T; 3],
}

pub fn intersecting_time_ee<T>(s: EdgeEdge<T>, e: EdgeEdge<T>, epsilon: T) -> Option<T>
where
    T: Copy + num_traits::Float + 'static + std::fmt::Debug + std::fmt::Display,
    i64: AsPrimitive<T>,
    f64: AsPrimitive<T>,
{
    use crate::vec3::Vec3;
    let list_te = coplanar_time(
        FourPoints {
            p0: s.a0,
            p1: s.a1,
            p2: s.b0,
            p3: s.b1,
        },
        FourPoints {
            p0: e.a0,
            p1: e.a1,
            p2: e.b0,
            p3: e.b1,
        },
        epsilon,
    );
    for te in list_te {
        let ts = T::one() - te;
        let a0 = s.a0.scale(ts).add(&e.a0.scale(te));
        let a1 = s.a1.scale(ts).add(&e.a1.scale(te));
        let b0 = s.b0.scale(ts).add(&e.b0.scale(te));
        let b1 = s.b1.scale(ts).add(&e.b1.scale(te));
        let coord = crate::edge3::intersection_edge3_when_coplanar(&a0, &a1, &b0, &b1);
        let Some(coord) = coord else {
            continue;
        }; // coplanar case
        if coord.0 >= T::zero()
            && coord.1 >= T::zero()
            && coord.2 >= T::zero()
            && coord.3 >= T::zero()
        {
            return Some(te);
        }
    }
    None
}

// --------------------
// below CCD using mat3x4_array_of_cols

/// Maximum distance between any column in `0..SPLIT` and any column in `SPLIT..4`.
/// for CCD
pub fn max_relative_norm<T, const SPLIT: usize>(x: &[[T; 3]; 4]) -> T
where
    T: num_traits::Float,
{
    let mut max_sq = T::zero();
    for i in 0..SPLIT {
        for j in SPLIT..4 {
            let du = [x[i][0] - x[j][0], x[i][1] - x[j][1], x[i][2] - x[j][2]];
            let sq = du[0] * du[0] + du[1] * du[1] + du[2] * du[2];
            if sq > max_sq {
                max_sq = sq;
            }
        }
    }
    max_sq.sqrt()
}

pub fn normalize_centerized_configuration<T>(
    x: &mut [[T; 3]; 4],
    dx: &mut [[T; 3]; 4],
    max_t: T,
) -> Option<T>
where
    T: Float + FromPrimitive,
{
    // max over x(t=0) and x(t=max_t)
    let mut max_entry = T::zero();

    for i in 0..4 {
        for j in 0..3 {
            max_entry = max_entry.max(x[i][j].abs());

            let x1 = x[i][j] + max_t * dx[i][j];
            max_entry = max_entry.max(x1.abs());
        }
    }

    let c099 = T::from_f64(0.99).unwrap();
    let mut scale = c099 / max_entry;

    if scale <= T::zero() || !scale.is_finite() {
        return None;
    }

    // gradually scale with 8 to avoid floating point error
    let mut scaled = T::one();
    let s = T::from_f64(8.0).unwrap();

    loop {
        let scale_now = if scale > s {
            scale = scale / s;
            s
        } else {
            scale
        };

        crate::mat3x4_array_of_cols::scale_in_place(x, scale_now);
        crate::mat3x4_array_of_cols::scale_in_place(dx, scale_now);
        scaled = scaled * scale_now;

        if scale_now != s {
            break;
        }
    }

    Some(scaled)
}

pub fn coplanar_time_mat3x4<T>(x0: &[[T; 3]; 4], dx: &[[T; 3]; 4]) -> Vec<T>
where
    T: Copy + num_traits::Float + 'static + std::fmt::Debug + std::fmt::Display,
    i64: AsPrimitive<T>,
{
    use crate::vec3::Vec3;
    let x1 = x0[1].sub(&x0[0]);
    let x2 = x0[2].sub(&x0[0]);
    let x3 = x0[3].sub(&x0[0]);
    let v1 = dx[1].sub(&dx[0]);
    let v2 = dx[2].sub(&dx[0]);
    let v3 = dx[3].sub(&dx[0]);
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
    crate::polynomial_root::cubic_roots_in_range_zero_to_t(k0, k1, k2, k3, T::one(), T::zero())
}

pub trait SquareDistance<T> {
    fn eval_with_dir(&self, x: &[[T; 3]; 4]) -> (T, [T; 3]);

    fn eval(&self, x: &[[T; 3]; 4]) -> T;
}

pub struct EdgeEdgeSquaredDist;

impl<T> SquareDistance<T> for EdgeEdgeSquaredDist
where
    T: Float,
{
    fn eval_with_dir(&self, x: &[[T; 3]; 4]) -> (T, [T; 3]) {
        let [_sqdist, s, t] = crate::edge3::nearest_to_edge3_accurate(&x[0], &x[1], &x[2], &x[3]);
        let pa = crate::edge3::lerp(&x[0], &x[1], s);
        let pb = crate::edge3::lerp(&x[2], &x[3], t);
        let dir = crate::vec3::sub(&pb, &pa);
        let sqdist = crate::vec3::dot(&dir, &dir);
        (sqdist, dir)
    }

    fn eval(&self, x: &[[T; 3]; 4]) -> T {
        self.eval_with_dir(x).0
    }
}

fn directional_advance<const SPLIT: usize>(
    x: &[[f32; 3]; 4],
    dx: &[[f32; 3]; 4],
    w: &[f32; 3],
    w_norm: f32,
    dip: f32,
    park: f32,
    max_t: f32,
) -> f32 {
    const ROUND_SLACK: f32 = 2.0 * f32::EPSILON;
    const COORD_SLACK: f32 = 6.2e-7;

    let w_abs = [w[0].abs(), w[1].abs(), w[2].abs()];
    let dip_w = dip * w_norm;
    let park_w = park * w_norm;
    let mut step = f32::MAX;

    for i in 0..SPLIT {
        for j in SPLIT..4 {
            let gap = [x[j][0] - x[i][0], x[j][1] - x[i][1], x[j][2] - x[i][2]];
            let rel = [
                dx[i][0] - dx[j][0],
                dx[i][1] - dx[j][1],
                dx[i][2] - dx[j][2],
            ];

            let gap_abs = [gap[0].abs(), gap[1].abs(), gap[2].abs()];
            let rel_abs = [rel[0].abs(), rel[1].abs(), rel[2].abs()];

            let w_dot_gap_abs =
                w_abs[0] * gap_abs[0] + w_abs[1] * gap_abs[1] + w_abs[2] * gap_abs[2];
            let w_dot_rel_abs =
                w_abs[0] * rel_abs[0] + w_abs[1] * rel_abs[1] + w_abs[2] * rel_abs[2];

            let err = ROUND_SLACK * (w_dot_gap_abs + max_t * w_dot_rel_abs + park_w)
                + COORD_SLACK * w_norm;

            let projected = w[0] * gap[0] + w[1] * gap[1] + w[2] * gap[2] - err;

            if projected <= dip_w {
                return 0.0;
            }

            let rate = w[0] * rel[0] + w[1] * rel[1] + w[2] * rel[2];
            if rate > 0.0 {
                step = step.min((projected - park_w).max(0.0) / rate);
            }
        }
    }
    step
}

pub fn maximum_safe_time_accurate<F, const SPLIT: usize>(
    x0: &[[f32; 3]; 4],
    dx: &[[f32; 3]; 4],
    u_max: f32,
    square_dist_func: &F,
    offset: f32,
    park_clearance: f32,
    max_t: f32,
) -> f32
where
    F: SquareDistance<f32>,
{
    let mut x = *x0;

    let mut lower_t = 0f32;
    let mut toi = 0f32;

    let (mut sqdist, mut w) = square_dist_func.eval_with_dir(x0);

    let offset_squared = offset * offset;

    // The pair must start separated.
    if sqdist <= offset_squared {
        return 0f32;
    }

    let park = offset + park_clearance;
    let dip = offset + park_clearance * 0.5;

    // Equivalent to fminf(park * park, nextafterf(d2, 0.0f))
    let park_squared = park * park; //.min(sqdist.next_down());

    let inv_u_max = 1.0 / u_max;

    const MAX_PROBES: usize = 4096;
    let mut probes = 0usize;

    while sqdist > park_squared {
        lower_t = toi;

        let dist = sqdist.sqrt();

        // Direction-agnostic conservative advance.
        let step = {
            let mut step = (dist - dip) * inv_u_max;
            if toi + step <= max_t {
                let directional_step = directional_advance::<SPLIT>(
                    &x,
                    dx,
                    &w,
                    dist * (1.0 + 2.0 * f32::EPSILON),
                    dip,
                    park,
                    max_t,
                );
                step = step.max(directional_step);
            }
            step
        };

        toi += step;

        if toi > max_t {
            return max_t;
        }

        probes += 1;
        if probes > MAX_PROBES {
            return lower_t;
        }

        x = crate::mat3x4_array_of_cols::add_scaled(x0, dx, toi);
        (sqdist, w) = square_dist_func.eval_with_dir(&x);
    }

    // Bisection to pin down lower_t (last safe time) vs upper_t (first contact).
    let mut upper_t = toi;
    let mut window = upper_t - lower_t;
    loop {
        toi = 0.5 * (upper_t + lower_t);
        x = crate::mat3x4_array_of_cols::add_scaled(x0, dx, toi);
        sqdist = square_dist_func.eval(&x);
        if sqdist > park_squared {
            lower_t = toi;
        } else {
            upper_t = toi;
        }
        let new_window = upper_t - lower_t;
        if new_window == window {
            break;
        }
        window = new_window;
    }

    {
        //dbg!(lower_t,upper_t,upper_t - lower_t,park_squared);
        let x_lower = crate::mat3x4_array_of_cols::add_scaled(x0, dx, lower_t);
        let x_upper = crate::mat3x4_array_of_cols::add_scaled(x0, dx, upper_t);
        let d2_lower = square_dist_func.eval(&x_lower);
        let d2_upper = square_dist_func.eval(&x_upper);
        assert!(d2_lower >= park_squared);
        assert!(park_squared >= d2_upper);
    }

    lower_t.max(0.0)
}

#[test]
fn hoge() {
    let x0 = [
        [-1., -1., -1.],
        [1., 1., -1.],
        [-1., 1., -1.],
        [-1., -1., 1.],
    ];
    let dx = [
        [2.1, 0., 2.3],
        [-2.1, 0., 2.1],
        [2.4, 0., 2.0],
        [2.2, 0., -2.1],
    ];
    //
    {
        let t = crate::ccd3::coplanar_time_mat3x4(&x0, &dx);
        dbg!(t);
    }
    //
    let u_max = max_relative_norm::<f32, 2>(&dx);
    let dist = EdgeEdgeSquaredDist;
    let dhat = 0.001;
    let t1 = maximum_safe_time_accurate::<EdgeEdgeSquaredDist, 2>(
        &x0,
        &dx,
        u_max,
        &dist,
        0.005,
        dhat * 0.01,
        1.0,
    );
    dbg!(&t1);
}
