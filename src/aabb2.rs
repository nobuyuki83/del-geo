//! methods for 2D Axis-aligned Bounding Box (AABB)

use num_traits::AsPrimitive;

/// signed distance from axis-aligned bounding box
/// * `pos_in` - where the signed distance is evaluated
/// * `x_min` - bounding box's x-coordinate minimum
/// * `x_max` - bounding box's x-coordinate maximum
/// * `y_min` - bounding box's y-coordinate minimum
/// * `y_max` - bounding box's y-coordinate maximum
/// * signed distance (inside is negative)
pub fn signed_distance_aabb<Real>(
    pos_in: nalgebra::Vector2::<Real>,
    min0: nalgebra::Vector2::<Real>,
    max0: nalgebra::Vector2::<Real>) -> Real
    where Real: nalgebra::RealField + Copy,
          f64: AsPrimitive<Real>
{
    let half = 0.5_f64.as_();
    let x_center = (max0.x + min0.x) * half;
    let y_center = (max0.y + min0.y) * half;
    let x_dist = (pos_in.x - x_center).abs() - (max0.x - min0.x) * half;
    let y_dist = (pos_in.y - y_center).abs() - (max0.y - min0.y) * half;
    x_dist.max(y_dist)
}