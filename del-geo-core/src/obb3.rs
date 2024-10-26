//! 3D Oriented Bounding Box (OBB)

use rand::distributions::{Distribution, Standard};

pub fn from_random<RAND, Real>(reng: &mut RAND) -> [Real; 12]
where
    RAND: rand::Rng,
    Real: num_traits::Float,
    Standard: Distribution<Real>,
{
    let one = Real::one();
    let aabb_m1p1 = [-one, -one, -one, one, one, one];
    let cntr = crate::aabb3::sample(&aabb_m1p1, reng);
    let u = crate::aabb3::sample(&aabb_m1p1, reng);
    let v = crate::aabb3::sample(&aabb_m1p1, reng);
    let v = crate::vec3::orthogonalize(&u, &v);
    let w = crate::aabb3::sample(&aabb_m1p1, reng);
    let w = crate::vec3::orthogonalize(&u, &w);
    let w = crate::vec3::orthogonalize(&v, &w);
    [
        cntr[0], cntr[1], cntr[2], u[0], u[1], u[2], v[0], v[1], v[2], w[0], w[1], w[2],
    ]
}

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

/// return the normalized axes and the magnitude of each axis
pub fn unit_axes_and_half_edge_lengths<Real>(obb: &[Real; 12]) -> ([[Real; 3]; 3], [Real; 3])
where
    Real: num_traits::Float,
{
    let l0 = (obb[3] * obb[3] + obb[4] * obb[4] + obb[5] * obb[5]).sqrt();
    let l1 = (obb[6] * obb[6] + obb[7] * obb[7] + obb[8] * obb[8]).sqrt();
    let l2 = (obb[9] * obb[9] + obb[10] * obb[10] + obb[11] * obb[11]).sqrt();
    let l0_inv = Real::one() / l0;
    let l1_inv = Real::one() / l1;
    let l2_inv = Real::one() / l2;
    let axes = [
        [obb[3] * l0_inv, obb[4] * l0_inv, obb[5] * l0_inv],
        [obb[6] * l1_inv, obb[7] * l1_inv, obb[8] * l1_inv],
        [obb[9] * l2_inv, obb[10] * l2_inv, obb[11] * l2_inv],
    ];
    let sizes = [l0, l1, l2];
    (axes, sizes)
}

/// Projection of an OBB at axis, return (min,max)
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

/// Find the distance of two ranges, return None if they are overlapped
fn distance_between_two_ranges<Real>(a: (Real, Real), b: (Real, Real)) -> Option<Real>
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

pub fn distance_to_obb3<Real>(obb_i: &[Real; 12], obb_j: &[Real; 12]) -> Real
where
    Real: num_traits::Float,
{
    let center_i: [Real; 3] = obb_i[0..3].try_into().unwrap();
    let axis_size_i = unit_axes_and_half_edge_lengths(obb_i);
    let center_j: [Real; 3] = obb_j[0..3].try_into().unwrap();
    let axis_size_j = unit_axes_and_half_edge_lengths(obb_j);
    let mut max_dist = Real::zero();

    for i in 0..3 {
        let axis_i = axis_size_i.0[i];
        let lh_i = axis_size_i.1[i];
        let c_i = crate::vec3::dot(&axis_i, &center_i);
        let range_i = (c_i - lh_i, c_i + lh_i);
        let range_j = range_axis(obb_j, &axis_i);
        let Some(dist) = distance_between_two_ranges(range_i, range_j) else {
            continue;
        };
        if dist > max_dist {
            max_dist = dist;
        }
    }
    for j in 0..3 {
        let axis_j = axis_size_j.0[j];
        let lh_j = axis_size_j.1[j];
        let c_j = crate::vec3::dot(&axis_j, &center_j);
        let range_j = (c_j - lh_j, c_j + lh_j);
        let range_i = range_axis(obb_i, &axis_j);
        let Some(dist) = distance_between_two_ranges(range_i, range_j) else {
            continue;
        };
        if dist > max_dist {
            max_dist = dist;
        }
    }
    for i in 0..3 {
        let axis_i = axis_size_i.0[i];
        for j in 0..3 {
            let axis_j = axis_size_j.0[j];
            let axis = crate::vec3::cross(&axis_i, &axis_j);
            let range_i = range_axis(obb_i, &axis);
            let range_j = range_axis(obb_j, &axis);
            let Some(dist) = distance_between_two_ranges(range_i, range_j) else {
                continue;
            };
            if dist > max_dist {
                max_dist = dist;
            }
        }
    }
    max_dist
}

pub fn nearest_to_point3<Real>(obb: &[Real; 12], p: &[Real; 3]) -> [Real; 3]
where Real: num_traits::Float
{
    if is_include_point(&obb, p, Real::zero()) {
        return *p;
    }
    let (axes, hlen) = unit_axes_and_half_edge_lengths(obb);
    let d = [p[0] - obb[0], p[1] - obb[1], p[2] - obb[2]];
    let t0 = crate::vec3::dot(&axes[0], &d).clamp(-hlen[0], hlen[0]);
    let t1 = crate::vec3::dot(&axes[1], &d).clamp(-hlen[1], hlen[1]);
    let t2 = crate::vec3::dot(&axes[2], &d).clamp(-hlen[2], hlen[2]);
    [
        obb[0] + t0 * axes[0][0] + t1 * axes[1][0] + t2 * axes[2][0],
        obb[1] + t0 * axes[0][1] + t1 * axes[1][1] + t2 * axes[2][1],
        obb[2] + t0 * axes[0][2] + t1 * axes[1][2] + t2 * axes[2][2],
    ]
}

#[test]
fn test_nearest_to_point3() {
    use rand::SeedableRng;
    let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    for _itr in 0..100 {
        let obb = from_random::<_, f64>(&mut reng);
        let p = crate::aabb3::sample(&[-1., -1., -1., 1., 1., 1.], &mut reng);
        let p_near = nearest_to_point3(&obb, &p);
        for _iter in 0..10 {
            let eps = 1.0e-2;
            let dp = crate::aabb3::sample(&[-eps, -eps, -eps, eps, eps, eps], &mut reng);
            let q = [
                p_near[0] + dp[0],
                p_near[1] + dp[1],
                p_near[2] + dp[2],
            ];
            let q = nearest_to_point3(&obb, &q);
            let len0 = crate::edge3::length(&p, &p_near);
            let len1 = crate::edge3::length(&p, &q);
            dbg!(len0, len1);
            assert!(len0<=len1);
        }
    }
}

/// Use Separating Axis Theorem (SAT) to check if two OBBs are intersected
pub fn is_intersect_to_obb3<Real>(obb_i: &[Real; 12], obb_j: &[Real; 12]) -> bool
where
    Real: num_traits::Float,
{
    let axes = {
        let (axes_i, _) = unit_axes_and_half_edge_lengths(&obb_i);
        let (axes_j, _) = unit_axes_and_half_edge_lengths(&obb_j);
        [
            axes_i[0], axes_i[1], axes_i[2], axes_j[0], axes_j[1], axes_j[2],
        ]
    };
    for axis in axes.iter() {
        let range_i = range_axis(obb_i, axis);
        let range_j = range_axis(obb_j, axis);
        if distance_between_two_ranges(range_i, range_j).is_some() {
            return false;
        }
    }
    true
}

#[test]
fn test_is_intersect_to_obb3() {
    use rand::SeedableRng;
    let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    for _iter in 0..100 {
        let obb_i = from_random::<_, f64>(&mut reng);
        let obb_j = from_random(&mut reng);
        let res0 = is_intersect_to_obb3(&obb_i, &obb_j);
        let p0 = arrayref::array_ref![obb_i, 0, 3]; // center
        let p1 = nearest_to_point3(&obb_j, p0);
        let p2 = nearest_to_point3(&obb_i, &p1);
        let p3 = nearest_to_point3(&obb_j, &p2);
        let p4 = nearest_to_point3(&obb_i, &p3);
        let p5 = nearest_to_point3(&obb_j, &p4);
        let p6 = nearest_to_point3(&obb_i, &p5);
        let len45 = crate::edge3::length(&p4, &p5);
        let len56 = crate::edge3::length(&p5, &p6);
        assert!(len56 <= len45);
        if len56 > 0. && len56 < len45 {
            continue;
        } // still converging
        let res1 = len56 < 0.0001;
        assert_eq!(res0, res1, "{} {}", len45, len56);
    }
}

#[test]
fn test2_is_intersect_to_obb3() {
    {
        // Seperated
        let obb1: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let obb2: [f64; 12] = [3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        assert_eq!(is_intersect_to_obb3(&obb1, &obb2), false);
    }
    {
        // Orthogonal vectors
        let random_i: [f64; 3] = [1.0, 3.0, 6.1];
        let random_j: [f64; 3] = [1.0, 7.0, 4.1];
        let i_basic = crate::vec3::normalized(&random_i);
        let mut j_basic = crate::vec3::normalized(&random_j);
        let k_basic = crate::vec3::cross(&i_basic, &j_basic);
        j_basic = crate::vec3::cross(&k_basic, &i_basic);
        let obb1: [f64; 12] = [
            0.0, 0.0, 0.0, i_basic[0], i_basic[1], i_basic[2], j_basic[0], j_basic[1], j_basic[2],
            k_basic[0], k_basic[1], k_basic[2],
        ];
        let obb2: [f64; 12] = [3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        assert_eq!(is_intersect_to_obb3(&obb1, &obb2), false);
        //
        let obb1: [f64; 12] = [
            0.0, 0.0, 0.0, i_basic[0], i_basic[1], i_basic[2], j_basic[0], j_basic[1], j_basic[2],
            k_basic[0], k_basic[1], k_basic[2],
        ];
        let obb2: [f64; 12] = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        assert_eq!(is_intersect_to_obb3(&obb1, &obb2), true);
        //
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
        assert_eq!(is_intersect_to_obb3(&obb1, &obb2), true);
    }
    {
        // Partially intersected
        let obb1: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let obb2: [f64; 12] = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        assert_eq!(is_intersect_to_obb3(&obb1, &obb2), true);
    }
    {
        // Constains
        let obb1: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let obb2: [f64; 12] = [0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5];
        assert_eq!(is_intersect_to_obb3(&obb1, &obb2), true);
    }
}

/// Get 8 corners of an obb3
/// The first four points defines the front face of the obb3, the last four defines back face.
pub fn corner_points<Real>(obb: &[Real; 12]) -> [[Real; 3]; 8]
where
    Real: num_traits::Float,
{
    [
        [
            obb[0] - obb[3] - obb[6] - obb[9],
            obb[1] - obb[4] - obb[7] - obb[10],
            obb[2] - obb[5] - obb[8] - obb[11],
        ],
        [
            obb[0] + obb[3] - obb[6] - obb[9],
            obb[1] + obb[4] - obb[7] - obb[10],
            obb[2] + obb[5] - obb[8] - obb[11],
        ],
        [
            obb[0] + obb[3] + obb[6] - obb[9],
            obb[1] + obb[4] + obb[7] - obb[10],
            obb[2] + obb[5] + obb[8] - obb[11],
        ],
        [
            obb[0] - obb[3] + obb[6] - obb[9],
            obb[1] - obb[4] + obb[7] - obb[10],
            obb[2] - obb[5] + obb[8] - obb[11],
        ],
        [
            obb[0] - obb[3] - obb[6] + obb[9],
            obb[1] - obb[4] - obb[7] + obb[10],
            obb[2] - obb[5] - obb[8] + obb[11],
        ],
        [
            obb[0] + obb[3] - obb[6] + obb[9],
            obb[1] + obb[4] - obb[7] + obb[10],
            obb[2] + obb[5] - obb[8] + obb[11],
        ],
        [
            obb[0] + obb[3] + obb[6] + obb[9],
            obb[1] + obb[4] + obb[7] + obb[10],
            obb[2] + obb[5] + obb[8] + obb[11],
        ],
        [
            obb[0] - obb[3] + obb[6] + obb[9],
            obb[1] - obb[4] + obb[7] + obb[10],
            obb[2] - obb[5] + obb[8] + obb[11],
        ],
    ]
}

/*
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
 */

/*
fn convex_sets_distance<Real, const M: usize, const N: usize>(
    s1: &[Real; M],
    s2: &[Real; N],
) -> Real
where
    Real: num_traits::Float,
{
    todo!()
}
 */

/*
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
 */

// #[test]
// fn test_convex_sets_distance() {}
