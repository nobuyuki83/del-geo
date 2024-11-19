/// ray_dir is not always a unit vector
pub fn intersection_ray<T>(
    center: &nalgebra::Vector3<T>,
    rad: T,
    ray_src: &nalgebra::Vector3<T>,
    ray_dir: &nalgebra::Vector3<T>,
) -> Option<T>
where
    T: nalgebra::RealField + Copy,
{
    // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    let op = ray_src - center;
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
