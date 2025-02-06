//! methods for 3D Continuous Collision Detection (CCD)

use num_traits::AsPrimitive;

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
