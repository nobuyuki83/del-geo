//! 3D Oriented Bounding Box (OBB)

pub fn is_include_point<Real>(obb: &[Real; 12], p: &[Real; 3], eps: Real) -> bool
where
    Real: num_traits::Float,
{
    let s = Real::one() + eps;
    let d = [p[0] - obb[0], p[1] - obb[1], p[2] - obb[2]];
    {
        let lx = obb[3] * obb[3] + obb[4] * obb[4] + obb[5] * obb[5];
        let dx = obb[3] * d[0] + obb[4] * d[1] + obb[5] * d[2];
        if dx.abs() > lx * s {
            return false;
        }
    }
    {
        let ly = obb[6] * obb[6] + obb[7] * obb[7] + obb[8] * obb[8];
        let dy = obb[6] * d[0] + obb[7] * d[1] + obb[8] * d[2];
        if dy.abs() > ly * s {
            return false;
        }
    }
    {
        let lz = obb[9] * obb[9] + obb[10] * obb[10] + obb[11] * obb[11];
        let dz = obb[9] * d[0] + obb[10] * d[1] + obb[11] * d[2];
        if dz.abs() > lz * s {
            return false;
        }
    }
    true
}

// return the normalized axes and the magnitude of each axis
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

// Projection of an OBB at axis, return (min,max)
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

// Find the distance of two ranges, return None if they are overlapped
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

// Use Separating Axis Theorem to check if two OBBs are intersected
pub fn distance_to_obb3<Real>(obb_i: &[Real; 12], obb_j: &[Real; 12]) -> Real
where
    Real: num_traits::Float,
{
    let center_i: [Real; 3] = obb_i[0..3].try_into().unwrap();
    let axis_size_i = axis_size(obb_i);
    let center_j: [Real; 3] = obb_j[0..3].try_into().unwrap();
    let axis_size_j = axis_size(obb_j);
    let mut max_dist = Real::zero();

    for i in 0..3 {
        let axis_i = axis_size_i.0[i];
        let lh_i = axis_size_i.1[i];
        let c_i = crate::vec3::dot(&axis_i, &center_i);
        let range_i = (c_i - lh_i, c_i + lh_i);
        let range_j = range_axis(obb_j, &axis_i);
        let Some(dist) = distance_range(range_i, range_j) else {
            continue;
        };
        max_dist = max_dist.max(dist);
    }
    for j in 0..3 {
        let axis_j = axis_size_j.0[j];
        let lh_j = axis_size_j.1[j];
        let c_j = crate::vec3::dot(&axis_j, &center_j);
        let range_j = (c_j - lh_j, c_j + lh_j);
        let range_i = range_axis(obb_i, &axis_j);
        let Some(dist) = distance_range(range_i, range_j) else {
            continue;
        };
        max_dist = max_dist.max(dist);
    }
    for i in 0..3 {
        let axis_i = axis_size_i.0[i];
        for j in 0..3 {
            let axis_j = axis_size_j.0[j];
            let axis = crate::vec3::cross(&axis_i, &axis_j);
            let range_i = range_axis(obb_i, &axis);
            let range_j = range_axis(obb_j, &axis);
            let Some(dist) = distance_range(range_i, range_j) else {
                continue;
            };
            max_dist = max_dist.max(dist);
        }
    }
    max_dist
}

pub fn is_intersect<Real>(obb_i: &[Real; 12], obb_j: &[Real; 12]) -> bool
where
    Real: num_traits::Float,
{
    // If all distances between OBBs are not bigger than zero, intersected
    distance_to_obb3(obb_i, obb_j) == Real::zero()
}

/// Get 8 corners of an obb3
/// The first four points defines the front face of the obb3, the last four defines back face.
pub fn corner_points<Real>(obb: &[Real; 12]) -> [[Real; 3]; 8]
where
    Real: num_traits::Float,
{
    // extend respect each axis
    let e_x = obb[3] + obb[6] + obb[9];
    let e_y = obb[4] + obb[7] + obb[10];
    let e_z = obb[5] + obb[8] + obb[11];
    [
        [obb[0] + e_x, obb[1] + e_y, obb[2] + e_z],
        [obb[0] + e_x, obb[1] - e_y, obb[2] + e_z],
        [obb[0] - e_x, obb[1] + e_y, obb[2] + e_z],
        [obb[0] - e_x, obb[1] - e_y, obb[2] + e_z],
        [obb[0] + e_x, obb[1] + e_y, obb[2] - e_z],
        [obb[0] + e_x, obb[1] - e_y, obb[2] - e_z],
        [obb[0] - e_x, obb[1] + e_y, obb[2] - e_z],
        [obb[0] - e_x, obb[1] - e_y, obb[2] - e_z],
    ]
}

fn nearest_point3<Real>(obb: &[Real; 12], p: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    use crate::{mat3_array_of_array::inverse, mat3_array_of_array::matmul};

    let tran = [
        [obb[3], obb[6], obb[9]],
        [obb[4], obb[7], obb[10]],
        [obb[5], obb[8], obb[11]],
    ];

    // TODO: check not inversable
    let inv_tran = inverse(&tran);

    let zero = Real::zero();
    // TODO: use matxi mul vector
    let p_mat = [[p[0], zero, zero], [p[1], zero, zero], [p[2], zero, zero]];
    // transfrom point to obb's column space
    let mut c_coord = matmul(&inv_tran, &p_mat);

    let one = Real::one();
    c_coord[0][0] = c_coord[0][0].clamp(-one, one);
    c_coord[1][0] = c_coord[1][0].clamp(-one, one);
    c_coord[2][0] = c_coord[2][0].clamp(-one, one);

    // convert to original system
    let o_coord = matmul(&tran, &c_coord);
    [o_coord[0][0], o_coord[1][0], o_coord[2][0]]
}

fn convex_sets_distance<Real, const M: usize, const N: usize>(
    s1: &[Real; M],
    s2: &[Real; N],
) -> Real
where
    Real: num_traits::Float,
{
    todo!()
}

#[test]
fn test_nearest_point_obb3() {
    let halfsqrt2 = 1.41421356 / 2.;
    let obb: [f64; 12] = [
        0.0, 0.0, 0.0, halfsqrt2, halfsqrt2, 0.0, -halfsqrt2, halfsqrt2, 0.0, 0.0, 0.0, 1.0,
    ];
    let p: [f64; 3] = [2., 0., 0.];
    let nearset = nearest_point3(&obb, &p);
    assert_eq!(nearset[0], halfsqrt2 * 2.);

    let obb: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let p: [f64; 3] = [2., 0., 0.];
    let nearset = nearest_point3(&obb, &p);
    assert_eq!(nearset[0], 1.);

    let obb: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 1.0];
    let p: [f64; 3] = [2., 5., 0.];
    let nearset = nearest_point3(&obb, &p);
    assert_eq!(nearset[1], 2.);
}

#[test]
fn test_convex_sets_distance() {}

#[test]
fn test_distance() {
    // Orthogonal vectors
    let random_i: [f64; 3] = [1.0, 3.0, 6.1];
    let random_j: [f64; 3] = [1.0, 7.0, 4.1];
    let i_basic = crate::vec3::normalized(&random_i);
    let mut j_basic = crate::vec3::normalized(&random_j);
    let k_basic = crate::vec3::cross(&i_basic, &j_basic);
    j_basic = crate::vec3::cross(&k_basic, &i_basic);

    // Seperated
    let obb1: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let obb2: [f64; 12] = [3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

    assert!(is_intersect(&obb1, &obb2) == false);

    let obb1: [f64; 12] = [
        0.0, 0.0, 0.0, i_basic[0], i_basic[1], i_basic[2], j_basic[0], j_basic[1], j_basic[2],
        k_basic[0], k_basic[1], k_basic[2],
    ];
    let obb2: [f64; 12] = [3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    assert!(is_intersect(&obb1, &obb2) == false);

    // Partially intersected
    let obb1: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let obb2: [f64; 12] = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    assert!(is_intersect(&obb1, &obb2));
    let obb1: [f64; 12] = [
        0.0, 0.0, 0.0, i_basic[0], i_basic[1], i_basic[2], j_basic[0], j_basic[1], j_basic[2],
        k_basic[0], k_basic[1], k_basic[2],
    ];
    let obb2: [f64; 12] = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    assert!(is_intersect(&obb1, &obb2));

    // Constains
    let obb1: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let obb2: [f64; 12] = [0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5];
    assert!(is_intersect(&obb1, &obb2));
    let obb1: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let obb2: [f64; 12] = [
        0.0,
        0.0,
        0.0,
        i_basic[0] * 0.5,
        i_basic[1] * 0.5,
        i_basic[2] * 0.5,
        j_basic[0] * 0.5,
        j_basic[1] * 0.5,
        j_basic[2] * 0.5,
        k_basic[0] * 0.5,
        k_basic[1] * 0.5,
        k_basic[2] * 0.5,
    ];

    assert!(is_intersect(&obb1, &obb2));
}
