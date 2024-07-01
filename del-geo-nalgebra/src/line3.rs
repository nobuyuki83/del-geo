//! methods for 3D line

pub fn nearest_to_line3<T>(
    line_org_a: &nalgebra::Vector3<T>,
    line_dir_a: &nalgebra::Vector3<T>,
    line_org_b: &nalgebra::Vector3<T>,
    line_dir_b: &nalgebra::Vector3<T>,
) -> (T, T, T, nalgebra::Vector3<T>, nalgebra::Vector3<T>)
where
    T: nalgebra::RealField + Copy,
{
    let xaa = line_dir_a.dot(line_dir_a);
    let xab = line_dir_b.dot(line_dir_a);
    let xbb = line_dir_b.dot(line_dir_b);
    let scale = xaa * xbb - xab * xab;
    let xac = line_dir_a.dot(&(line_org_b - line_org_a));
    let xbc = line_dir_b.dot(&(line_org_b - line_org_a));
    let da = xbb * xac - xab * xbc;
    let db = xab * xac - xaa * xbc;
    let scaled_neraest_a = line_org_a.scale(scale) + line_dir_a.scale(da);
    let scaled_nearest_b = line_org_b.scale(scale) + line_dir_b.scale(db);
    (scale, da, db, scaled_neraest_a, scaled_nearest_b)
}
