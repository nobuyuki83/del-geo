//! methods for unit 3D vector

pub fn map_to_unit2_octahedron(dir: &[f32; 3]) -> [f32; 2] {
    let n = dir[0].abs() + dir[1].abs() + dir[2].abs();
    let oct = [dir[0] / n, dir[1] / n];
    let oct = if dir[2] < 0. {
        [
            (1. - oct[1].abs()) * oct[0].signum(),
            (1. - oct[0].abs()) * oct[1].signum(),
        ]
    } else {
        oct
    };
    [oct[0] * 0.5 + 0.5, oct[1] * 0.5 + 0.5]
}

/// <https://github.com/mmp/pbrt-v4/blob/1ae72cfa7344e79a7815a21ed3da746cdccee59b/src/pbrt/util/math.cpp#L317>
pub fn map_to_unit2_equal_area(d: &[f32; 3]) -> [f32; 2] {
    let x = d[0].abs();
    let y = d[1].abs();
    let z = d[2].abs();
    let r = (1. - z).sqrt();
    let phi = y.atan2(x);
    let phi = phi * std::f32::consts::FRAC_2_PI;
    let v = phi * r;
    let u = r - v;
    let (u, v) = if d[2] < 0. { (1. - v, 1. - u) } else { (u, v) };
    let u = u.copysign(d[0]);
    let v = v.copysign(-d[1]);
    [u * 0.5 + 0.5, v * 0.5 + 0.5]
}
