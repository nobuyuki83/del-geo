//! methods for 2D Axis-aligned Bounding Box (AABB)

use num_traits::AsPrimitive;
use rand::distributions::Standard;
use rand::prelude::Distribution;

pub fn from_vtx2xy<Real>(vtx2xy: &[Real]) -> [Real; 4]
where
    Real: num_traits::Float,
{
    let mut aabb = [vtx2xy[0], vtx2xy[1], vtx2xy[0], vtx2xy[1]];
    vtx2xy.chunks(2).skip(1).for_each(|v| {
        aabb[0] = if v[0] < aabb[0] { v[0] } else { aabb[0] };
        aabb[1] = if v[1] < aabb[1] { v[1] } else { aabb[1] };
        aabb[2] = if v[0] > aabb[2] { v[0] } else { aabb[2] };
        aabb[3] = if v[1] > aabb[3] { v[1] } else { aabb[3] };
    });
    aabb
}

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
    from_vtx2xy(&[p0[0], p0[1], p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]])
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

pub fn from_list_of_vertices<Index, T>(idx2vtx: &[Index], vtx2xy: &[T], eps: T) -> [T; 4]
where
    T: num_traits::Float,
    Index: AsPrimitive<usize>,
{
    assert!(!idx2vtx.is_empty());
    let mut aabb = [T::zero(); 4];
    {
        let i_vtx: usize = idx2vtx[0].as_();
        {
            let cgx = vtx2xy[i_vtx * 2];
            aabb[0] = cgx - eps;
            aabb[2] = cgx + eps;
        }
        {
            let cgy = vtx2xy[i_vtx * 2 + 1];
            aabb[1] = cgy - eps;
            aabb[3] = cgy + eps;
        }
    }
    for &i_vtx in idx2vtx.iter().skip(1) {
        let i_vtx = i_vtx.as_();
        {
            let cgx = vtx2xy[i_vtx * 2];
            aabb[0] = if cgx - eps < aabb[0] {
                cgx - eps
            } else {
                aabb[0]
            };
            aabb[2] = if cgx + eps > aabb[2] {
                cgx + eps
            } else {
                aabb[2]
            };
        }
        {
            let cgy = vtx2xy[i_vtx * 2 + 1];
            aabb[1] = if cgy - eps < aabb[1] {
                cgy - eps
            } else {
                aabb[1]
            };
            aabb[3] = if cgy + eps > aabb[3] {
                cgy + eps
            } else {
                aabb[3]
            };
        }
    }
    assert!(aabb[0] <= aabb[2]);
    assert!(aabb[1] <= aabb[3]);
    aabb
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
