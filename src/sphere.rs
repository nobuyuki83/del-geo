pub fn intersection_ray(
    center: &nalgebra::Vector3::<f32>,
    rad: f32,
    ray_src: &nalgebra::Vector3::<f32>,
    ray_dir: &nalgebra::Vector3::<f32>)
    -> Option<(nalgebra::Vector3<f32>, nalgebra::Vector3<f32>, f32)>
{
    let depth0 = (center - ray_src).dot(ray_dir);
    if depth0 < 0f32 { return None; }
    let sqdist = (ray_src + depth0 * ray_dir - center).norm_squared();
    if rad * rad - sqdist < 0f32 { return None; }
    let depth1 = depth0 - (rad * rad - sqdist).sqrt();
    if depth1 < 0f32 { return None; }
    let hit_pos = ray_src + depth1 * ray_dir;
    let hit_normal = (hit_pos - center).normalize();
    Some((hit_pos, hit_normal, depth1))
}