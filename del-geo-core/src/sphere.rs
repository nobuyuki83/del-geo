pub fn intersection_ray<T>(rad: T, center: &[T; 3], ray_src: &[T; 3], ray_dir: &[T; 3]) -> Option<T>
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    let op = ray_src.sub(center);
    let a = ray_dir.dot(ray_dir);
    let b = op.dot(ray_dir);
    let c = op.dot(&op) - rad * rad;
    let det: T = b * b - c * a;
    if det < T::zero() {
        None
    } else {
        let det = det.sqrt();
        if -b - det >= T::zero() {
            return Some((-b - det) / a);
        } else if -b + det >= T::zero() {
            return Some((-b + det) / a);
        }
        None
    }
}

pub fn area(r: f32) -> f32 {
    r * r * 4f32 * std::f32::consts::PI
}

/// <https://corysimon.github.io/articles/uniformdistn-on-sphere/>
pub fn sample(rnd: &[f32; 2]) -> [f32; 3] {
    let phi = (1. - 2. * rnd[0]).acos();
    let theta = 2. * std::f32::consts::PI * rnd[1];
    [theta.cos() * phi.sin(), theta.sin() * phi.sin(), phi.cos()]
}

/// uniformly region on sample unit sphere (center is the origin) where the other sphere is visible
pub fn sample_where_another_sphere_is_visible(
    rad_light: f32,
    pos_light_center: &[f32; 3],
    unirand: &[f32; 2],
) -> ([f32; 3], f32) {
    use crate::mat3_col_major::Mat3ColMajor;
    use crate::vec3::Vec3;
    let sin_theta_max_squared = rad_light * rad_light / pos_light_center.squared_norm();
    assert!((0f32..=1f32).contains(&sin_theta_max_squared));
    let cos_theta_max = (1f32 - sin_theta_max_squared).max(0.0).sqrt();
    let cos_theta = 1f32 - unirand[0] * (1f32 - cos_theta_max);
    assert!(cos_theta >= 0.0);
    assert!(
        cos_theta >= cos_theta_max,
        "{} {}",
        cos_theta,
        cos_theta_max
    );
    assert!(cos_theta <= 1f32);
    let sin_theta = (1f32 - cos_theta * cos_theta).max(0f32).sqrt();
    let phi = 2f32 * std::f32::consts::PI * unirand[1];
    // sample unit sphere assuming that the other light is in the z-axis direction
    let dir_lcl = [sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta];
    let mat3 = crate::mat3_col_major::transform_lcl2world_given_local_z(pos_light_center);
    let dir_world = mat3.mult_vec(&dir_lcl);
    let pdf = 1f32 / (2f32 * std::f32::consts::PI * (1f32 - cos_theta_max));
    (dir_world, pdf)
}
pub fn pdf_light_sample(light_center: &[f32; 3], light_rad: f32) -> f32 {
    use crate::vec3::Vec3;
    let sin_theta_max_squared = light_rad * light_rad / light_center.squared_norm();
    assert!(sin_theta_max_squared > 0f32 && sin_theta_max_squared < 1f32);
    let cos_theta_max = (1f32 - sin_theta_max_squared).max(0.0).sqrt();
    1f32 / (2f32 * std::f32::consts::PI * (1f32 - cos_theta_max))
}
