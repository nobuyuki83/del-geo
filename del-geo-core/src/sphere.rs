use crate::vec3::Vec3;
pub fn intersection_ray<T>(rad: T, center: &[T; 3], ray_src: &[T; 3], ray_dir: &[T; 3]) -> Option<T>
where
    T: num_traits::Float + Copy + std::ops::MulAssign,
{
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
