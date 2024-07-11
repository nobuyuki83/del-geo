//! methods for 3D Axis-aligned Bounding Box (AABB)

use num_traits::AsPrimitive;
use rand::distributions::{Distribution, Standard};

pub fn center<T>(aabb: &[T; 6]) -> [T; 3]
where
    T: num_traits::Float + 'static + Copy,
    f64: AsPrimitive<T>,
{
    [
        (aabb[0] + aabb[3]) * 0.5f64.as_(),
        (aabb[1] + aabb[4]) * 0.5f64.as_(),
        (aabb[2] + aabb[5]) * 0.5f64.as_(),
    ]
}

pub fn max_edge_size<T>(aabb: &[T; 6]) -> T
where
    T: num_traits::Float,
{
    let lx = aabb[3] - aabb[0];
    let ly = aabb[4] - aabb[1];
    let lz = aabb[5] - aabb[2];
    if lx > ly {
        return if lx > lz { lx } else { lz };
    }
    if ly > lz {
        return ly;
    }
    lz
}

pub fn from_two_aabbs<T>(i0: &[T; 6], i1: &[T; 6]) -> [T; 6]
where
    T: num_traits::Float,
{
    assert_eq!(i0.len(), 6);
    assert_eq!(i1.len(), 6);
    let mut o = [T::zero(); 6];
    for i in 0..3 {
        o[i] = if i0[i] < i1[i] { i0[i] } else { i1[i] };
        o[i + 3] = if i0[i + 3] > i1[i + 3] {
            i0[i + 3]
        } else {
            i1[i + 3]
        };
    }
    o
}

pub fn is_active<T>(i0: &[T; 6]) -> bool
where
    T: PartialOrd,
{
    i0[0] <= i0[3]
}

pub fn is_intersect<T>(i0: &[T; 6], i1: &[T; 6]) -> bool
where
    T: PartialOrd,
{
    assert_eq!(i0.len(), 6);
    assert_eq!(i1.len(), 6);
    if !is_active(i0) {
        return false;
    }
    if !is_active(i1) {
        return false;
    }
    if i0[0] > i1[3] {
        return false;
    }
    if i0[1] > i1[4] {
        return false;
    }
    if i0[2] > i1[5] {
        return false;
    }
    if i0[3] < i1[0] {
        return false;
    }
    if i0[4] < i1[1] {
        return false;
    }
    if i0[5] < i1[2] {
        return false;
    }
    true
}

/// transform aabb to unit square (0,1)^3 while preserving aspect ratio
/// return 4x4 homogeneous transformation matrix in **column major** order
pub fn to_transformation_world2unit_ortho_preserve_asp(aabb_world: &[f32; 6]) -> [f32; 16] {
    let cntr = [
        (aabb_world[0] + aabb_world[3]) * 0.5,
        (aabb_world[1] + aabb_world[4]) * 0.5,
        (aabb_world[2] + aabb_world[5]) * 0.5,
    ];
    let size = max_edge_size(aabb_world);
    let a = 1f32 / size;
    let b = -a * cntr[0] + 0.5;
    let c = -a * cntr[1] + 0.5;
    let d = -a * cntr[2] + 0.5;
    [a, 0., 0., 0., 0., a, 0., 0., 0., 0., a, 0., b, c, d, 1.]
}

pub fn sample<Reng, T>(aabb: &[T; 6], reng: &mut Reng) -> [T; 3]
where
    Reng: rand::Rng,
    T: num_traits::Float,
    Standard: Distribution<T>,
{
    let r0 = reng.gen::<T>();
    let r1 = reng.gen::<T>();
    let r2 = reng.gen::<T>();
    [
        aabb[0] + r0 * (aabb[3] - aabb[0]),
        aabb[1] + r1 * (aabb[4] - aabb[1]),
        aabb[2] + r2 * (aabb[5] - aabb[2]),
    ]
}

// --------------------------

pub type AABB3<'a, Real> = crate::aabb::AABB<'a, Real, 3, 6>;

pub fn from_slice<Real>(s: &[Real]) -> AABB3<Real>
where
    Real: num_traits::Float,
{
    AABB3 {
        aabb: arrayref::array_ref!(s, 0, 6),
    }
}

pub fn from_aabbs<Real>(aabbs: &[Real], i_aabb: usize) -> AABB3<Real>
where
    Real: num_traits::Float,
{
    AABB3 {
        aabb: arrayref::array_ref!(aabbs, i_aabb * 6, 6),
    }
}
