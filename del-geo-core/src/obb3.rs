#![allow(unused_imports)]
//! 3D Oriented Bounding Box (OBB)

use crate::vec3::Vec3;
use crate::{aabb3::AABB3Trait, edge3::Edge3Trait};
//use rand::distributions::{Distribution, Standard};
/// trait for 3D Oriented Bounding Box (OBB)
pub trait OBB3Trait<T> {
    fn is_include_point(&self, p: &[T; 3], eps: T) -> bool;
    fn unit_axes_and_half_edge_lengths(&self) -> ([[T; 3]; 3], [T; 3]);
    fn nearest_to_point3(&self, p: &[T; 3]) -> [T; 3];
}
impl<Real> OBB3Trait<Real> for [Real; 12]
where
    Real: num_traits::Float,
{
    fn is_include_point(&self, p: &[Real; 3], eps: Real) -> bool {
        is_include_point(self, p, eps)
    }
    fn unit_axes_and_half_edge_lengths(&self) -> ([[Real; 3]; 3], [Real; 3]) {
        unit_axes_and_half_edge_lengths(self)
    }
    fn nearest_to_point3(&self, p: &[Real; 3]) -> [Real; 3] {
        nearest_to_point3(self, p)
    }
}

pub fn from_random<RAND, Real>(reng: &mut RAND) -> [Real; 12]
where
    RAND: rand::Rng,
    Real: num_traits::Float,
    rand::distr::StandardUniform: rand::distr::Distribution<Real>,
{
    let one = Real::one();
    let aabb_m1p1 = [-one, -one, -one, one, one, one];
    let cntr = aabb_m1p1.sample(reng);
    let u = aabb_m1p1.sample(reng);
    let v = aabb_m1p1.sample(reng);
    let v = u.orthogonalize(&v);
    let w = aabb_m1p1.sample(reng);
    let w = u.orthogonalize(&w);
    let w = v.orthogonalize(&w);
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
fn range_axis<Real, const N: usize>(ps: &[[Real; 3]; N], axis: &[Real; 3]) -> (Real, Real)
where
    Real: num_traits::Float,
{
    let min0 = ps
        .iter()
        .map(|v| v.dot(axis))
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let max0 = ps
        .iter()
        .map(|v| v.dot(axis))
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    (min0, max0)
}

/*
pub fn distance_to_obb3<Real>(obb_i: &[Real; 12], obb_j: &[Real; 12]) -> Real
where
    Real: num_traits::Float,
{
    let center_i: [Real; 3] = obb_i[0..3].try_into().unwrap();
    let axis_size_i = unit_axes_and_half_edge_lengths(obb_i);
    let center_j: [Real; 3] = obb_j[0..3].try_into().unwrap();
    let axis_size_j = unit_axes_and_half_edge_lengths(obb_j);
    let mut max_dist = Real::zero();
    let cp_i = corner_points(obb_i);
    let cp_j = corner_points(obb_j);

    for i in 0..3 {
        let axis_i = axis_size_i.0[i];
        let lh_i = axis_size_i.1[i];
        let c_i = crate::vec3::dot(&axis_i, &center_i);
        let range_i = (c_i - lh_i, c_i + lh_i);
        let range_j = range_axis(&cp_j, &axis_i);
        let Some(dist) = distance_between_two_ranges(range_i, range_j) else {
            continue;
        };
        max_dist = max_dist.max(dist);
    }
    for j in 0..3 {
        let axis_j = axis_size_j.0[j];
        let lh_j = axis_size_j.1[j];
        let c_j = crate::vec3::dot(&axis_j, &center_j);
        let range_i = range_axis(&cp_i, &axis_j);
        let range_j = (c_j - lh_j, c_j + lh_j);
        let Some(dist) = distance_between_two_ranges(range_i, range_j) else {
            continue;
        };
        max_dist = max_dist.max(dist);
    }
    for i in 0..3 {
        let axis_i = axis_size_i.0[i];
        for j in 0..3 {
            let axis_j = axis_size_j.0[j];
            let axis = crate::vec3::cross(&axis_i, &axis_j);
            let range_i = range_axis(&cp_i, &axis);
            let range_j = range_axis(&cp_j, &axis);
            let Some(dist) = distance_between_two_ranges(range_i, range_j) else {
                continue;
            };
            max_dist = max_dist.max(dist);
        }
    }
    max_dist
}
 */

pub fn nearest_to_point3<Real>(obb: &[Real; 12], p: &[Real; 3]) -> [Real; 3]
where
    Real: num_traits::Float,
{
    if obb.is_include_point(p, Real::zero()) {
        return *p;
    }
    let (axes, hlen) = obb.unit_axes_and_half_edge_lengths();
    let d = p.sub(obb[..3].try_into().unwrap());
    let [t0, t1, t2] = std::array::from_fn::<_, 3, _>(|i| axes[i].dot(&d).clamp(-hlen[i], hlen[i]));
    axes[0]
        .scale(t0)
        .add(&axes[1].scale(t1))
        .add(&axes[2].scale(t2))
        .add(&obb[..3].try_into().unwrap())
}

#[test]
fn test_nearest_to_point3() {
    use rand::SeedableRng;
    let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    for _itr in 0..1000 {
        let obb = from_random::<_, f64>(&mut reng);
        let p = [-1., -1., -1., 1., 1., 1.].sample(&mut reng);
        let p_near = obb.nearest_to_point3(&p);
        for _iter in 0..100 {
            let eps = 1.0e-4;
            let dp = [-eps, -eps, -eps, eps, eps, eps].sample(&mut reng);
            let q = p_near.add(&dp);
            let q = obb.nearest_to_point3(&q);
            let len0 = p.length(&p_near);
            let len1 = p.length(&q);
            assert!(len0 <= len1);
        }
    }
}

/// Use Separating Axis Theorem (SAT) to check if two OBBs are intersected
pub fn is_intersect_to_obb3<Real>(obb_i: &[Real; 12], obb_j: &[Real; 12]) -> bool
where
    Real: num_traits::Float + std::fmt::Debug,
{
    let axes = {
        let (axes_i, _) = obb_i.unit_axes_and_half_edge_lengths();
        let (axes_j, _) = obb_j.unit_axes_and_half_edge_lengths();
        [
            axes_i[0],
            axes_i[1],
            axes_i[2],
            axes_j[0],
            axes_j[1],
            axes_j[2],
            axes_i[0].cross(&axes_j[0]),
            axes_i[0].cross(&axes_j[1]),
            axes_i[0].cross(&axes_j[2]),
            axes_i[1].cross(&axes_j[0]),
            axes_i[1].cross(&axes_j[1]),
            axes_i[1].cross(&axes_j[2]),
            axes_i[2].cross(&axes_j[0]),
            axes_i[2].cross(&axes_j[1]),
            axes_i[2].cross(&axes_j[2]),
        ]
    };
    let corner_i = corner_points(obb_i);
    let corner_j = corner_points(obb_j);
    for axis in axes.iter() {
        let range_i = range_axis(&corner_i, axis);
        let range_j = range_axis(&corner_j, axis);
        if crate::range::distance_to_range(range_i, range_j).is_some() {
            return false;
        }
    }
    true
}

#[test]
fn test_is_intersect_to_obb3() {
    use rand::SeedableRng;
    let mut reng = rand_chacha::ChaChaRng::seed_from_u64(0u64);
    for _iter in 0..1000 {
        let obb_i = from_random::<_, f64>(&mut reng);
        let obb_j = from_random(&mut reng);
        let p0 = obb_i[..3].try_into().unwrap(); // center
        let p1 = obb_j.nearest_to_point3(p0);
        let p2 = obb_i.nearest_to_point3(&p1);
        let p3 = obb_j.nearest_to_point3(&p2);
        let p4 = obb_i.nearest_to_point3(&p3);
        let p5 = obb_j.nearest_to_point3(&p4);
        let p6 = obb_i.nearest_to_point3(&p5);
        let len45 = p4.length(&p5);
        let len56 = p5.length(&p6);
        assert!(len56 <= len45);
        if len56 > 0. && len56 < len45 * 0.9999 {
            continue;
        } // still converging
        {
            // test intersect
            let res0 = is_intersect_to_obb3(&obb_i, &obb_j);
            let res1 = len56 < 0.0001;
            if res0 != res1 {
                let (mut tri2vtx_i, mut vtx2xyz_i) = del_msh_cpu::trimesh3_primitive::obb3(&obb_i);
                let (tri2vtx_j, vtx2xyz_j) = del_msh_cpu::trimesh3_primitive::obb3(&obb_j);
                del_msh_cpu::uniform_mesh::merge(
                    &mut tri2vtx_i,
                    &mut vtx2xyz_i,
                    &tri2vtx_j,
                    &vtx2xyz_j,
                    3,
                );
                // output mesh to visualize the failure case
                let _ = del_msh_cpu::io_obj::save_tri2vtx_vtx2xyz(
                    "../../target/fail_obb3.obj",
                    &tri2vtx_i,
                    &vtx2xyz_i,
                    3,
                );
            }
            assert_eq!(res0, res1, "{} {}", len45, len56);
        }
    }
}

#[test]
fn test2_is_intersect_to_obb3() {
    use std::f64::consts::PI;

    for i in 0..2 {
        let obb_i = [0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 1.];

        // special obb need to check additional axes other than six orientation axes
        // this test will fail if not checking addtional axes
        let d = 2.0 + 0.01 * (i as f64); // small distance to the obbs make them seperated
        let quad = crate::quaternion::around_axis(&[1., 0., -1.], -PI * 0.5f64);
        let rot_mat = crate::quaternion::to_mat3_col_major(&quad);
        let u = crate::mat3_col_major::mult_vec(&rot_mat, &[1., 0., 0.]);
        let v = crate::mat3_col_major::mult_vec(&rot_mat, &[0., 1., 0.]);
        let w = crate::mat3_col_major::mult_vec(&rot_mat, &[0., 0., 1.]);
        let obb_j = [
            d, 0., -d, u[0], u[1], u[2], v[0], v[1], v[2], w[0], w[1], w[2],
        ];

        let p0 = obb_i[..3].try_into().unwrap(); // center
        let p1 = obb_j.nearest_to_point3(p0);
        let p2 = obb_i.nearest_to_point3(&p1);
        let p3 = obb_j.nearest_to_point3(&p2);
        let p4 = obb_i.nearest_to_point3(&p3);
        let p5 = obb_j.nearest_to_point3(&p4);
        let p6 = obb_i.nearest_to_point3(&p5);

        let len45 = p4.length(&p5);
        let len56 = p5.length(&p6);
        assert!(len56 <= len45);

        let res1 = len56 < 0.0001; // this will be false since not intersected
        let res0 = is_intersect_to_obb3(&obb_i, &obb_j);
        assert_eq!(res0, res1, "{} {}", len45, len56);
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
