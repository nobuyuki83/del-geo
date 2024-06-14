//! methods for 2D Axis-aligned Bounding Box (AABB)

use num_traits::AsPrimitive;

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
    if lx > ly {
        lx
    } else {
        ly
    }
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


/// transform aabb to unit square (0,1)^2 while preserving aspect ratio
/// return 3x3 homogeneous transformation matrix in **column major** order
pub fn transform_world2unit_ortho_preserve_asp(
    aabb_world: &[f32; 4],
) -> [f32; 9] {
    let width_world = aabb_world[2] - aabb_world[0];
    let height_world = aabb_world[3] - aabb_world[1];
    let cntr_world = [
        (aabb_world[0] + aabb_world[2]) * 0.5,
        (aabb_world[1] + aabb_world[3]) * 0.5,
    ];
    let aabb_world1 = if (width_world / height_world) > 1f32 {
        [
            aabb_world[0],
            cntr_world[1] - width_world * 0.5,
            aabb_world[2],
            cntr_world[1] + width_world * 0.5,
        ]
    } else {
        [
            cntr_world[0] - height_world * 0.5,
            aabb_world[1],
            cntr_world[0] + height_world * 0.5,
            aabb_world[3],
        ]
    };
    // dbg!(&aabb_world1);
    let p_00 = [aabb_world1[0], aabb_world1[1]];
    let p_11 = [aabb_world1[2], aabb_world1[3]];
    let a = 1f32 / (p_11[0] - p_00[0]);
    let c = -a * p_00[0];
    let b = 1f32 / (p_11[1] - p_00[1]);
    let d = -b * p_00[1];
    [a, 0., 0., 0., b, 0., c, d, 1.]
}


#[allow(clippy::identity_op)]
pub fn from_list_of_vertices<Index, T>(idx2vtx: &[Index], vtx2xy: &[T], eps: T) -> [T; 4]
    where
        T: num_traits::Float,
        Index: AsPrimitive<usize>
{
    assert!(!idx2vtx.is_empty());
    let mut aabb = [T::zero(); 4];
    {
        let i_vtx: usize = idx2vtx[0].as_();
        {
            let cgx = vtx2xy[i_vtx * 2 + 0];
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
            let cgx = vtx2xy[i_vtx * 2 + 0];
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

#[allow(clippy::identity_op)]
pub fn from_two_aabbs<T>(i0: &[T; 4], i1: &[T; 4]) -> [T; 4]
    where
        T: num_traits::Float,
{
    let mut o = [T::zero(); 4];
    for i in 0..2 {
        o[i + 0] = if i0[i + 0] < i1[i + 0] {
            i0[i + 0]
        } else {
            i1[i + 0]
        };
        o[i + 2] = if i0[i + 2] > i1[i + 2] {
            i0[i + 2]
        } else {
            i1[i + 2]
        };
    }
    o
}

// ----------------------
// below: interface includes nalgebra

/// signed distance from axis-aligned bounding box
/// * `pos_in` - where the signed distance is evaluated
/// * `x_min` - bounding box's x-coordinate minimum
/// * `x_max` - bounding box's x-coordinate maximum
/// * `y_min` - bounding box's y-coordinate minimum
/// * `y_max` - bounding box's y-coordinate maximum
/// * signed distance (inside is negative)
pub fn signed_distance_aabb<Real>(
    pos_in: nalgebra::Vector2<Real>,
    min0: nalgebra::Vector2<Real>,
    max0: nalgebra::Vector2<Real>,
) -> Real
where
    Real: nalgebra::RealField + Copy,
    f64: AsPrimitive<Real>,
{
    let half = 0.5_f64.as_();
    let x_center = (max0.x + min0.x) * half;
    let y_center = (max0.y + min0.y) * half;
    let x_dist = (pos_in.x - x_center).abs() - (max0.x - min0.x) * half;
    let y_dist = (pos_in.y - y_center).abs() - (max0.y - min0.y) * half;
    x_dist.max(y_dist)
}

pub fn from_vtx2vec<T>(vtx2vec: &[nalgebra::Vector2<T>]) -> [T; 4]
where
    T: nalgebra::RealField + Copy,
{
    let mut aabb = [vtx2vec[0][0], vtx2vec[0][1], vtx2vec[0][0], vtx2vec[0][1]];
    for xy in vtx2vec.iter().skip(1) {
        if xy[0] < aabb[0] {
            aabb[0] = xy[0];
        }
        if xy[1] < aabb[1] {
            aabb[1] = xy[1];
        }
        if xy[0] > aabb[2] {
            aabb[2] = xy[0];
        }
        if xy[1] > aabb[3] {
            aabb[3] = xy[1];
        }
    }
    aabb
}
