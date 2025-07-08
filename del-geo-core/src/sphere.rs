//! Methods for sphere.
//! The sphere is represented as `(radius: Real, center: &[Real;3])`

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

pub fn area<T>(r: T) -> T
where
    T: num_traits::Float + num_traits::FloatConst,
{
    r * r * T::from(4).unwrap() * T::PI()
}

/// <https://corysimon.github.io/articles/uniformdistn-on-sphere/>
pub fn sample_surface_uniform<T>(rnd: &[T; 2]) -> [T; 3]
where
    T: num_traits::Float + num_traits::FloatConst,
{
    let one = T::one();
    let two = one + one;
    let phi = (one - two * rnd[0]).acos();
    let theta = two * T::PI() * rnd[1];
    [theta.cos() * phi.sin(), theta.sin() * phi.sin(), phi.cos()]
}

/// uniformly region on sample unit sphere (center is the origin) where the other sphere is visible
pub fn sample_where_another_sphere_is_visible<T>(
    rad_light: T,
    pos_light_center: &[T; 3],
    unirand: &[T; 2],
) -> ([T; 3], T)
where
    T: num_traits::Float + num_traits::FloatConst + std::fmt::Display,
{
    use crate::mat3_col_major::Mat3ColMajor;
    use crate::vec3::Vec3;
    let one = T::one();
    let zero = T::zero();
    let two = one + one;
    let sin_theta_max_squared = rad_light * rad_light / pos_light_center.squared_norm();
    assert!((zero..=one).contains(&sin_theta_max_squared));
    let cos_theta_max = (one - sin_theta_max_squared).max(zero).sqrt();
    let cos_theta = one - unirand[0] * (one - cos_theta_max);
    assert!(cos_theta >= zero);
    assert!(
        cos_theta >= cos_theta_max,
        "{} {}",
        cos_theta,
        cos_theta_max
    );
    assert!(cos_theta <= one);
    let sin_theta = (one - cos_theta * cos_theta).max(zero).sqrt();
    let phi = two * T::PI() * unirand[1];
    // sample unit sphere assuming that the other light is in the z-axis direction
    let dir_lcl = [sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta];
    let mat3 = crate::mat3_col_major::transform_lcl2world_given_local_z(pos_light_center);
    let dir_world = mat3.mult_vec(&dir_lcl);
    let pdf = one / (two * T::PI() * (one - cos_theta_max));
    (dir_world, pdf)
}
pub fn pdf_light_sample<Real>(light_center: &[Real; 3], light_rad: Real) -> Real
where
    Real: num_traits::Float + num_traits::FloatConst,
{
    use crate::vec3::Vec3;
    let zero = Real::zero();
    let one = Real::one();
    let two = one + one;
    let sin_theta_max_squared = light_rad * light_rad / light_center.squared_norm();
    assert!(sin_theta_max_squared > zero && sin_theta_max_squared < one);
    let cos_theta_max = (one - sin_theta_max_squared).max(zero).sqrt();
    one / (two * Real::PI() * (one - cos_theta_max))
}
