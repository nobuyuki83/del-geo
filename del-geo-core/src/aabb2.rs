//! methods for 2D Axis-aligned Bounding Box (AABB)

use num_traits::AsPrimitive;
use rand::distributions::Standard;
use rand::prelude::Distribution;

pub fn from_two_points<T>(p0: &[T; 2], p1: &[T; 2], rad: T) -> [T; 4]
where
    T: num_traits::Float,
{
    [
        (p0[0] - rad).min(p1[0] - rad),
        (p0[1] - rad).min(p1[1] - rad),
        (p0[0] + rad).max(p1[0] + rad),
        (p0[1] + rad).max(p1[1] + rad),
    ]
}

pub fn from_two_aabbs<T>(i0: &[T; 4], i1: &[T; 4]) -> [T; 4]
where
    T: num_traits::Float,
{
    let mut o = [T::zero(); 4];
    for i in 0..2 {
        o[i] = if i0[i] < i1[i] { i0[i] } else { i1[i] };
        o[i + 2] = if i0[i + 2] > i1[i + 2] {
            i0[i + 2]
        } else {
            i1[i + 2]
        };
    }
    o
}

pub fn rasterize<T>(aabb: &[T; 4], img_size: &(usize, usize)) -> [usize; 4]
where
    T: num_traits::Float + AsPrimitive<usize> + 'static + Copy,
    usize: AsPrimitive<T>,
{
    let half = T::one() / (T::one() + T::one());
    [
        (aabb[0] - half).ceil().max(T::zero()).as_(),
        (aabb[1] - half).ceil().max(T::zero()).as_(),
        (aabb[2] - half).ceil().min(img_size.0.as_()).as_(),
        (aabb[3] - half).ceil().min(img_size.1.as_()).as_(),
    ]
}

pub fn center<T>(aabb: &[T; 4]) -> [T; 2]
where
    T: num_traits::Float,
{
    let half = T::one() / (T::one() + T::one());
    [(aabb[0] + aabb[2]) * half, (aabb[1] + aabb[3]) * half]
}

pub fn max_edge_size<T>(aabb: &[T; 4]) -> T
where
    T: num_traits::Float,
{
    let lx = aabb[2] - aabb[0];
    let ly = aabb[3] - aabb[1];
    lx.max(ly)
}

pub fn transform_homogeneous<T>(aabb: &[T; 4], transform: &[T; 9]) -> [T; 4]
where
    T: num_traits::Float,
{
    let p0 = crate::mat3::transform_homogeneous(transform, &[aabb[0], aabb[1]]).unwrap();
    let p1 = crate::mat3::transform_homogeneous(transform, &[aabb[0], aabb[3]]).unwrap();
    let p2 = crate::mat3::transform_homogeneous(transform, &[aabb[2], aabb[1]]).unwrap();
    let p3 = crate::mat3::transform_homogeneous(transform, &[aabb[2], aabb[3]]).unwrap();
    let ax = [p0[0], p1[0], p2[0], p3[0]];
    let ay = [p0[1], p1[1], p2[1], p3[1]];
    [
        ax.iter().fold(T::infinity(), |a, &b| a.min(b)),
        ay.iter().fold(T::infinity(), |a, &b| a.min(b)),
        ax.iter().fold(-T::infinity(), |a, &b| a.max(b)),
        ay.iter().fold(-T::infinity(), |a, &b| a.max(b)),
    ]
}

pub fn sample<Reng, T>(aabb: &[T; 4], reng: &mut Reng) -> [T; 2]
where
    Reng: rand::Rng,
    T: num_traits::Float,
    Standard: Distribution<T>,
{
    let r0 = reng.gen::<T>();
    let r1 = reng.gen::<T>();
    [
        aabb[0] + r0 * (aabb[2] - aabb[0]),
        aabb[1] + r1 * (aabb[3] - aabb[1]),
    ]
}

/// transform aabb to unit square (0,1)^2 while preserving aspect ratio
/// return 3x3 homogeneous transformation matrix in **column major** order
pub fn to_transformation_world2unit_ortho_preserve_asp(aabb_world: &[f32; 4]) -> [f32; 9] {
    let cntr = center(aabb_world);
    let size = max_edge_size(aabb_world);
    let a = 1f32 / size;
    let b = -a * cntr[0] + 0.5;
    let c = -a * cntr[1] + 0.5;
    [a, 0., 0., 0., a, 0., b, c, 1.]
}

#[test]
fn test_to_transformation_world2unit_ortho_preserve_asp() {
    let aabb = [-4.1, 0.3, 4.5, 3.3];
    let transf = to_transformation_world2unit_ortho_preserve_asp(&aabb);
    let mut reng = rand::thread_rng();
    for _ in 0..1000 {
        let p0 = sample(&aabb, &mut reng);
        let q0 = crate::mat3::transform_homogeneous(&transf, &p0).unwrap();
        assert!(q0[0] >= 0.0 && q0[0] <= 1.0);
        assert!(q0[1] >= 0.0 && q0[1] <= 1.0);
    }
}

// -------------------------------------------------------------------------------

pub type AABB2<'a, Real> = crate::aabb::AABB<'a, Real, 2, 4>;

pub fn from_slice<Real>(s: &[Real]) -> AABB2<Real>
where
    Real: num_traits::Float,
{
    AABB2 {
        aabb: arrayref::array_ref!(s, 0, 4),
    }
}

pub fn from_aabbs<Real>(aabbs: &[Real], i_aabb: usize) -> AABB2<Real>
where
    Real: num_traits::Float,
{
    AABB2 {
        aabb: arrayref::array_ref!(aabbs, i_aabb * 4, 4),
    }
}
