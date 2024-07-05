//! method common in 2D or 3D Axis-Aligned Bounding Box (AABB)

use num_traits::Zero;
pub fn is_intersect_ray<const NDIM: usize, const SIZE_AABB: usize>(
    aabb: &[f32; SIZE_AABB],
    ray_org: &[f32; NDIM],
    ray_dir: &[f32; NDIM],
) -> bool {
    assert_eq!(NDIM * 2, SIZE_AABB);
    let mut tmin = -f32::INFINITY;
    let mut tmax = f32::INFINITY;
    for i_dim in 0..NDIM {
        if !ray_dir[i_dim].abs().is_zero() {
            let t1 = (aabb[i_dim] - ray_org[i_dim]) / ray_dir[i_dim];
            let t2 = (aabb[i_dim + NDIM] - ray_org[i_dim]) / ray_dir[i_dim];
            tmin = tmin.max(t1.min(t2));
            tmax = tmax.min(t1.max(t2));
        } else if ray_org[i_dim] < aabb[i_dim] || ray_org[i_dim] > aabb[i_dim + NDIM] {
            return false;
        }
    }
    tmax >= tmin && tmax >= 0.0
}

pub fn is_include_point<Real, const NDIM: usize, const SIZE_AABB: usize>(
    aabb: &[Real; SIZE_AABB],
    point: &[Real; NDIM],
) -> bool
where
    Real: num_traits::Float,
{
    assert_eq!(NDIM * 2, SIZE_AABB);
    for i_dim in 0..NDIM {
        if point[i_dim] < aabb[i_dim] || point[i_dim] > aabb[i_dim + NDIM] {
            return false;
        }
    }
    true
}
