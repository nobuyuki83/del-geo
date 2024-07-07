//!

use crate::vec3::cross;

pub fn is_include_point<Real>(obb: &[Real; 12], p: &[Real; 3]) -> bool
where
    Real: num_traits::Float,
{
    let d = [p[0] - obb[0], p[1] - obb[1], p[2] - obb[2]];
    {
        let lx = obb[3] * obb[3] + obb[4] * obb[4] + obb[5] * obb[5];
        let dx = obb[3] * d[0] + obb[4] * d[1] + obb[5] * d[2];
        if dx.abs() > lx {
            return false;
        }
    }
    {
        let ly = obb[6] * obb[6] + obb[7] * obb[7] + obb[8] * obb[8];
        let dy = obb[6] * d[0] + obb[7] * d[1] + obb[8] * d[2];
        if dy.abs() > ly {
            return false;
        }
    }
    {
        let lz = obb[9] * obb[9] + obb[10] * obb[10] + obb[11] * obb[11];
        let dz = obb[9] * d[0] + obb[10] * d[1] + obb[11] * d[2];
        if dz.abs() > lz {
            return false;
        }
    }
    true
}

pub fn axis_size<Real>(obb: &[Real; 12]) -> ([[Real; 3]; 3], [Real; 3])
where
    Real: num_traits::Float,
{
    let lx = (obb[3] * obb[3] + obb[4] * obb[4] + obb[5] * obb[5]).sqrt();
    let ly = (obb[6] * obb[6] + obb[7] * obb[7] + obb[8] * obb[8]).sqrt();
    let lz = (obb[9] * obb[9] + obb[10] * obb[10] + obb[11] * obb[11]).sqrt();
    let lx_inv = Real::one() / lx;
    let ly_inv = Real::one() / ly;
    let lz_inv = Real::one() / lz;
    let axes = [
        [obb[3] * lx_inv, obb[4] * lx_inv, obb[5] * lx_inv],
        [obb[6] * ly_inv, obb[7] * ly_inv, obb[8] * ly_inv],
        [obb[9] * lz_inv, obb[10] * lz_inv, obb[11] * lz_inv],
    ];
    let sizes = [lx, ly, lz];
    (axes, sizes)
}

pub fn range_axis<Real>(obb: &[Real; 12], axis: &[Real; 3]) -> (Real, Real)
where
    Real: num_traits::Float,
{
    let c = obb[0] * axis[0] + obb[1] * axis[1] + obb[2] * axis[2];
    let x = (obb[3] * axis[0] + obb[4] * axis[1] + obb[5] * axis[2]).abs();
    let y = (obb[6] * axis[0] + obb[7] * axis[1] + obb[8] * axis[2]).abs();
    let z = (obb[9] * axis[0] + obb[10] * axis[1] + obb[11] * axis[2]).abs();
    (c - x - y - z, c + x + y + z)
}

pub fn distance_range<Real>(a: (Real, Real), b: (Real, Real)) -> Option<Real>
where
    Real: num_traits::Float,
{
    if a.0 > b.1 {
        return Some(a.0 - b.1);
    }
    if b.0 > a.1 {
        return Some(b.0 - a.1);
    }
    None
}

pub fn distance_to_obb3<Real>(obb_i: &[Real; 12], obb_j: &[Real; 12]) -> Option<Real>
where
    Real: num_traits::Float,
{
    let center_i: [Real; 3] = obb_i[0..3].try_into().unwrap();
    let axis_size_i = axis_size(obb_i);
    let center_j: [Real; 3] = obb_j[0..3].try_into().unwrap();
    let axis_size_j = axis_size(obb_j);
    let mut min_dist = Real::max_value();
    for i in 0..3 {
        let axis_i = axis_size_i.0[i];
        let lh_i = axis_size_i.1[i];
        let c_i = crate::vec3::dot(&axis_i, &center_i);
        let range_i = (c_i - lh_i, c_i + lh_i);
        let range_j = range_axis(obb_j, &axis_i);
        let Some(dist) = distance_range(range_i, range_j) else {
            return None;
        };
        if dist < min_dist {
            min_dist = dist;
        }
    }
    for j in 0..3 {
        let axis_j = axis_size_j.0[j];
        let lh_j = axis_size_j.1[j];
        let c_j = crate::vec3::dot(&axis_j, &center_j);
        let range_j = (c_j - lh_j, c_j + lh_j);
        let range_i = range_axis(obb_i, &axis_j);
        let Some(dist) = distance_range(range_i, range_j) else {
            return None;
        };
        if dist < min_dist {
            min_dist = dist;
        }
    }
    for i in 0..3 {
        let axis_i = axis_size_i.0[i];
        for j in 0..3 {
            let axis_j = axis_size_j.0[j];
            let axis = crate::vec3::cross(&axis_i, &axis_j);
            let range_i = range_axis(obb_i, &axis);
            let range_j = range_axis(obb_j, &axis);
            let Some(dist) = distance_range(range_i, range_j) else {
                return None;
            };
            if dist < min_dist {
                min_dist = dist;
            }
        }
    }
    Some(min_dist)
}
