//! methods for 3D Axis-aligned Bounding Box (AABB)

use num_traits::AsPrimitive;

pub fn from_vtx2xyz<T>(vtx2xyz: &[T], eps: T) -> [T; 6]
where
    T: num_traits::Float,
{
    assert!(!vtx2xyz.is_empty());
    let mut aabb = [T::zero(); 6];
    {
        {
            let cgx = vtx2xyz[0];
            aabb[0] = cgx - eps;
            aabb[3] = cgx + eps;
        }
        {
            let cgy = vtx2xyz[1];
            aabb[1] = cgy - eps;
            aabb[4] = cgy + eps;
        }
        {
            let cgz = vtx2xyz[2];
            aabb[2] = cgz - eps;
            aabb[5] = cgz + eps;
        }
    }
    for i_vtx in 1..vtx2xyz.len() / 3 {
        {
            let cgx = vtx2xyz[i_vtx * 3];
            aabb[0] = if cgx - eps < aabb[0] {
                cgx - eps
            } else {
                aabb[0]
            };
            aabb[3] = if cgx + eps > aabb[3] {
                cgx + eps
            } else {
                aabb[3]
            };
        }
        {
            let cgy = vtx2xyz[i_vtx * 3 + 1];
            aabb[1] = if cgy - eps < aabb[1] {
                cgy - eps
            } else {
                aabb[1]
            };
            aabb[4] = if cgy + eps > aabb[4] {
                cgy + eps
            } else {
                aabb[4]
            };
        }
        {
            let cgz = vtx2xyz[i_vtx * 3 + 2];
            aabb[2] = if cgz - eps < aabb[2] {
                cgz - eps
            } else {
                aabb[2]
            };
            aabb[5] = if cgz + eps > aabb[5] {
                cgz + eps
            } else {
                aabb[5]
            };
        }
    }
    assert!(aabb[0] <= aabb[3]);
    assert!(aabb[1] <= aabb[4]);
    assert!(aabb[2] <= aabb[5]);
    aabb
}

pub fn from_list_of_vertices<Index, Real>(
    idx2vtx: &[Index],
    vtx2xyz: &[Real],
    eps: Real,
) -> [Real; 6]
where
    Real: num_traits::Float,
    Index: AsPrimitive<usize>,
{
    assert!(!idx2vtx.is_empty());
    let mut aabb = [Real::zero(); 6];
    {
        let i_vtx: usize = idx2vtx[0].as_();
        {
            let cgx = vtx2xyz[i_vtx * 3];
            aabb[0] = cgx - eps;
            aabb[3] = cgx + eps;
        }
        {
            let cgy = vtx2xyz[i_vtx * 3 + 1];
            aabb[1] = cgy - eps;
            aabb[4] = cgy + eps;
        }
        {
            let cgz = vtx2xyz[i_vtx * 3 + 2];
            aabb[2] = cgz - eps;
            aabb[5] = cgz + eps;
        }
    }
    for &i_vtx in idx2vtx.iter().skip(1) {
        let i_vtx: usize = i_vtx.as_();
        {
            let cgx = vtx2xyz[i_vtx * 3];
            aabb[0] = if cgx - eps < aabb[0] {
                cgx - eps
            } else {
                aabb[0]
            };
            aabb[3] = if cgx + eps > aabb[3] {
                cgx + eps
            } else {
                aabb[3]
            };
        }
        {
            let cgy = vtx2xyz[i_vtx * 3 + 1];
            aabb[1] = if cgy - eps < aabb[1] {
                cgy - eps
            } else {
                aabb[1]
            };
            aabb[4] = if cgy + eps > aabb[4] {
                cgy + eps
            } else {
                aabb[4]
            };
        }
        {
            let cgz = vtx2xyz[i_vtx * 3 + 2];
            aabb[2] = if cgz - eps < aabb[2] {
                cgz - eps
            } else {
                aabb[2]
            };
            aabb[5] = if cgz + eps > aabb[5] {
                cgz + eps
            } else {
                aabb[5]
            };
        }
    }
    assert!(aabb[0] <= aabb[3]);
    assert!(aabb[1] <= aabb[4]);
    assert!(aabb[2] <= aabb[5]);
    aabb
}

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
        if lx > lz {
            return lx;
        } else {
            return lz;
        }
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
    T: std::cmp::PartialOrd,
{
    i0[0] <= i0[3]
}

pub fn is_intersect<T>(i0: &[T; 6], i1: &[T; 6]) -> bool
where
    T: std::cmp::PartialOrd,
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
