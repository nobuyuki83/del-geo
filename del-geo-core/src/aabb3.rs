//! methods for 3D Axis-aligned Bounding Box (AABB)

use rand::distributions::{Distribution, Standard};

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

// Above: from method
// ----------------------------------
// Below: to method

/// transform aabb to unit square (0,1)^3 while preserving aspect ratio
/// return 4x4 homogeneous transformation matrix in **column major** order
pub fn to_mat4_col_major_transf_into_unit_preserve_asp(aabb_world: &[f32; 6]) -> [f32; 16] {
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


/// transform aabb to unit square (0,1)^3
/// return 4x4 homogeneous transformation matrix in **column major** order
pub fn to_mat4_col_major_transf_into_unit(aabb_world: &[f32; 6]) -> [f32; 16] {
    let cntr = crate::aabb3::center(aabb_world);
    let size = crate::aabb3::size(aabb_world);
    let ax = 1f32 / size[0];
    let ay = 1f32 / size[1];
    let az = 1f32 / size[2];
    let b = -ax * cntr[0] + 0.5;
    let c = -ay * cntr[1] + 0.5;
    let d = -az * cntr[2] + 0.5;
    [ax, 0., 0., 0., 0., ay, 0., 0., 0., 0., az, 0., b, c, d, 1.]
}

// ----------------------------------

pub fn scale<T>(aabb: &[T; 6], s: T) -> [T; 6]
where
    T: num_traits::Float,
{
    let c = center(aabb);
    let size = size(aabb);
    let half = T::one() / (T::one() + T::one());
    [
        c[0] - size[0] * s * half,
        c[1] - size[1] * s * half,
        c[2] - size[2] * s * half,
        c[0] + size[0] * s * half,
        c[1] + size[1] * s * half,
        c[2] + size[2] * s * half,
    ]
}

pub fn center<T>(aabb: &[T; 6]) -> [T; 3]
where
    T: num_traits::Float,
{
    let half = T::one() / (T::one() + T::one());
    [
        (aabb[0] + aabb[3]) * half,
        (aabb[1] + aabb[4]) * half,
        (aabb[2] + aabb[5]) * half,
    ]
}

pub fn size<T>(aabb: &[T; 6]) -> [T; 3]
where
    T: num_traits::Float,
{
    [aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2]]
}

pub fn volume<T>(aabb: &[T; 6]) -> T
where T: num_traits::Float
{
    (aabb[3] - aabb[0])*(aabb[4] - aabb[1])*(aabb[5] - aabb[2])
}

pub fn xyz_from_hex_index<Real>(aabb: &[Real; 6], i_vtx: usize) -> [Real; 3]
where
    Real: num_traits::Float,
{
    match i_vtx {
        0 => [aabb[0], aabb[1], aabb[2]],
        1 => [aabb[3], aabb[1], aabb[2]],
        2 => [aabb[3], aabb[4], aabb[2]],
        3 => [aabb[0], aabb[4], aabb[2]],
        4 => [aabb[0], aabb[1], aabb[5]],
        5 => [aabb[3], aabb[1], aabb[5]],
        6 => [aabb[3], aabb[4], aabb[5]],
        7 => [aabb[0], aabb[4], aabb[5]],
        _ => panic!(),
    }
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


/// return a vec3 sampled inside a aabb
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

pub fn set_as_cube<Real>(aabb: &mut [Real; 6], xyz: &[Real; 3], eps: Real)
where
    Real: num_traits::Float,
{
    {
        let cgx = xyz[0];
        aabb[0] = cgx - eps;
        aabb[3] = cgx + eps;
    }
    {
        let cgy = xyz[1];
        aabb[1] = cgy - eps;
        aabb[4] = cgy + eps;
    }
    {
        let cgz = xyz[2];
        aabb[2] = cgz - eps;
        aabb[5] = cgz + eps;
    }
}

pub fn update<Real>(aabb: &mut [Real; 6], xyz: &[Real; 3], eps: Real)
where
    Real: num_traits::Float,
{
    aabb[0] = aabb[0].min(xyz[0] - eps);
    aabb[3] = aabb[3].max(xyz[0] + eps);
    aabb[1] = aabb[1].min(xyz[1] - eps);
    aabb[4] = aabb[4].max(xyz[1] + eps);
    aabb[2] = aabb[2].min(xyz[2] - eps);
    aabb[5] = aabb[5].max(xyz[2] + eps);
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
