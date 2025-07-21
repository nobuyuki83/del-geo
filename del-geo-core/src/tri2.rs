//! methods for 2D triangle

use crate::mat2_col_major::add_four;
use crate::vec2::Vec2;

pub fn area<T>(p0: &[T; 2], p1: &[T; 2], p2: &[T; 2]) -> T
where
    T: num_traits::Float,
{
    let half = T::one() / (T::one() + T::one());
    half * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
}

/// # Return
/// `(dldp0: [T;2], dldp1: [T;2], dldp2: [T;2])`
pub fn dldw_area<T>(p0: &[T; 2], p1: &[T; 2], p2: &[T; 2], dldarea: T) -> ([T; 2], [T; 2], [T; 2])
where
    T: num_traits::Float,
{
    let two = T::one() + T::one();
    let half = T::one() / two;
    let dareadp0x2 = [p1[1] - p2[1], p2[0] - p1[0]];
    let dareadp1x2 = [p2[1] - p0[1], p0[0] - p2[0]];
    let dareadp2x2 = [p0[1] - p1[1], p1[0] - p0[0]];
    (
        dareadp0x2.scale(dldarea * half),
        dareadp1x2.scale(dldarea * half),
        dareadp2x2.scale(dldarea * half),
    )
}

#[test]
fn test_dldw_area() {
    let p0s = [[0.1, -0.2], [1.3, 0.2], [0.6, 0.45]];
    let dldarea = 1.3f32;
    let area0 = area(&p0s[0], &p0s[1], &p0s[2]);
    let l0 = dldarea * area0;
    let dl = dldw_area(&p0s[0], &p0s[1], &p0s[2], dldarea);
    let dl = [dl.0, dl.1, dl.2];
    let eps = 1.0e-4;
    for inode in 0..3 {
        for idim in 0..2 {
            let p1s = {
                let mut p1s = p0s;
                p1s[inode][idim] += eps;
                p1s
            };
            let area1 = area(&p1s[0], &p1s[1], &p1s[2]);
            let l1 = dldarea * area1;
            let val_num = (l1 - l0) / eps;
            let val_ana = dl[inode][idim];
            let diff_abs = (val_num - val_ana).abs();
            assert!(diff_abs < 5.0e-4, "{}", diff_abs);
        }
    }
}

pub fn is_inside<Real>(
    p0: &[Real; 2],
    p1: &[Real; 2],
    p2: &[Real; 2],
    q: &[Real; 2],
    sign: Real,
) -> Option<(Real, Real)>
where
    Real: num_traits::Float,
{
    let a0 = area(q, p1, p2) * sign;
    if a0 < Real::zero() {
        return None;
    }
    let a1 = area(q, p2, p0) * sign;
    if a1 < Real::zero() {
        return None;
    }
    let a2 = area(q, p0, p1) * sign;
    if a2 < Real::zero() {
        return None;
    }
    let sum_area_inv = Real::one() / (a0 + a1 + a2);
    Some((a0 * sum_area_inv, a1 * sum_area_inv))
}

/// shape function's derivative in x- and y-direction and the constant term
///
/// # Example
/// ```
/// let (dldx, c) = del_geo_core::tri2::dldx(&[0., 1.], &[1., 2.], &[3., 4.]);
/// let a = dldx[0][2]; // x-differentiation of L^2
/// ```
pub fn dldx<T>(p0: &[T; 2], p1: &[T; 2], p2: &[T; 2]) -> ([[T; 3]; 2], [T; 3])
where
    T: num_traits::Float,
{
    use crate::vec3::Vec3;
    let half = T::one() / (T::one() + T::one());
    let a0 = area(p0, p1, p2);
    let tmp1: T = half / a0;
    (
        [
            [p1[1] - p2[1], p2[1] - p0[1], p0[1] - p1[1]].scale(tmp1),
            [p2[0] - p1[0], p0[0] - p2[0], p1[0] - p0[0]].scale(tmp1),
        ],
        [
            p1[0] * p2[1] - p2[0] * p1[1],
            p2[0] * p0[1] - p0[0] * p2[1],
            p0[0] * p1[1] - p1[0] * p0[1],
        ]
        .scale(tmp1),
    )
}

pub fn barycentric_coords<Real>(
    p0: &[Real; 2],
    p1: &[Real; 2],
    p2: &[Real; 2],
    q: &[Real; 2],
) -> Option<(Real, Real, Real)>
where
    Real: num_traits::Float,
{
    let a0 = area(q, p1, p2);
    let a1 = area(q, p2, p0);
    let a2 = area(q, p0, p1);
    if (a0 + a1 + a2).is_zero() {
        return None;
    }
    let sum_area_inv = Real::one() / (a0 + a1 + a2);
    Some((a0 * sum_area_inv, a1 * sum_area_inv, a2 * sum_area_inv))
}

// -------------------------------------------
#[derive(Debug, Clone, Copy)]
pub struct Tri2<'a, Real> {
    pub p0: &'a [Real; 2],
    pub p1: &'a [Real; 2],
    pub p2: &'a [Real; 2],
}

impl<Real> Tri2<'_, Real>
where
    Real: num_traits::Float,
{
    pub fn is_inside(&self, q: &[Real; 2], sign: Real) -> Option<(Real, Real)> {
        let Tri2 { p0, p1, p2 } = self;
        is_inside(p0, p1, p2, q, sign)
    }

    pub fn area(&self) -> Real {
        let Tri2 { p0, p1, p2 } = self;
        area(p0, p1, p2)
    }
}

pub fn circumcenter<T>(p0: &[T; 2], p1: &[T; 2], p2: &[T; 2]) -> [T; 2]
where
    T: num_traits::Float + Copy + std::fmt::Debug,
{
    use crate::vec2::Vec2;
    let a0 = p1.sub(p2).squared_norm();
    let a1 = p2.sub(p0).squared_norm();
    let a2 = p0.sub(p1).squared_norm();

    let b0: T = a0 * (a1 + a2 - a0);
    let b1: T = a1 * (a0 + a2 - a1);
    let b2: T = a2 * (a0 + a1 - a2);
    let sum = T::one() / (b0 + b1 + b2);

    let c0 = b0 * sum;
    let c1 = b1 * sum;
    let c2 = b2 * sum;
    [
        p0[0] * c0 + p1[0] * c1 + p2[0] * c2,
        p0[1] * c0 + p1[1] * c1 + p2[1] * c2,
    ]
}

#[test]
fn test_circumcenter() {
    let p0 = [[1.3f64, 2.1], [3.2, 2.1], [1.5, 2.5]];
    let cc0 = circumcenter(&p0[0], &p0[1], &p0[2]);
    let d0 = cc0.sub(&p0[0]).squared_norm();
    let d1 = cc0.sub(&p0[1]).squared_norm();
    let d2 = cc0.sub(&p0[2]).squared_norm();
    assert!((d0 - d1).abs() < d0 * 1.0e-10);
    assert!((d0 - d2).abs() < d0 * 1.0e-10);
}

pub fn wdw_circumcenter<T>(p0: &[T; 2], p1: &[T; 2], p2: &[T; 2]) -> ([T; 2], [[T; 4]; 3])
where
    T: num_traits::Float + Copy + std::fmt::Debug + num_traits::Float,
{
    use crate::vec2::Vec2;
    use crate::vec2::XY;
    let p0 = XY { p: *p0 };
    let p1 = XY { p: *p1 };
    let p2 = XY { p: *p2 };
    let a0 = (p1 - p2).p.squared_norm();
    let a1 = (p2 - p0).p.squared_norm();
    let a2 = (p0 - p1).p.squared_norm();
    //
    let b0: T = a0 * (a1 + a2 - a0);
    let b1: T = a1 * (a2 + a0 - a1);
    let b2: T = a2 * (a0 + a1 - a2);
    //
    let sum = b0 + b1 + b2;
    let sum_inv = T::one() / sum;
    //
    let c0 = b0 * sum_inv;
    let c1 = b1 * sum_inv;
    let c2 = b2 * sum_inv;
    let cc = crate::vec2::add_three(&p0.p.scale(c0), &p1.p.scale(c1), &p2.p.scale(c2));
    // -----------------
    let two = T::one() + T::one();
    let db0 = [
        (p0 - p2 + p0 - p1).p.scale(two * a0),
        crate::vec2::add(
            &((p1 - p0 + p2 - p1).p.scale(two * a0)),
            &(p1 - p2).p.scale(two * (a1 + a2 - a0)),
        ),
        crate::vec2::add(
            &((p2 - p0 + p1 - p2).p.scale(two * a0)),
            &(p2 - p1).p.scale(two * (a1 + a2 - a0)),
        ),
    ];
    let db1 = [
        crate::vec2::add(
            &(p0 - p1 + p2 - p0).p.scale(two * a1),
            &(p0 - p2).p.scale(two * (a2 + a0 - a1)),
        ),
        (p1 - p0 + p1 - p2).p.scale(two * a1),
        crate::vec2::add(
            &(p2 - p1 + p0 - p2).p.scale(two * a1),
            &(p2 - p0).p.scale(two * (a2 + a0 - a1)),
        ),
    ];
    let db2 = [
        crate::vec2::add(
            &(p0 - p2 + p1 - p0).p.scale(two * a2),
            &(p0 - p1).p.scale(two * (a0 + a1 - a2)),
        ),
        crate::vec2::add(
            &(p1 - p2 + p0 - p1).p.scale(two * a2),
            &(p1 - p0).p.scale(two * (a0 + a1 - a2)),
        ),
        (p2 - p1 + p2 - p0).p.scale(two * a2),
    ];
    let tmp = -T::one() / (sum * sum);
    let dsum_inv = [
        crate::vec2::add_three(&db0[0], &db1[0], &db2[0]).scale(tmp),
        crate::vec2::add_three(&db0[1], &db1[1], &db2[1]).scale(tmp),
        crate::vec2::add_three(&db0[2], &db1[2], &db2[2]).scale(tmp),
    ];
    //
    use crate::mat2_col_major::{Mat2ColMajor, from_identity, from_outer_product};
    use crate::vec2;
    let dcc = [
        add_four(
            &from_identity().scale(c0),
            &from_outer_product(
                &p0.p,
                &vec2::add(&db0[0].scale(sum_inv), &dsum_inv[0].scale(b0)),
            ),
            &from_outer_product(
                &p1.p,
                &vec2::add(&db1[0].scale(sum_inv), &dsum_inv[0].scale(b1)),
            ),
            &from_outer_product(
                &p2.p,
                &vec2::add(&db2[0].scale(sum_inv), &dsum_inv[0].scale(b2)),
            ),
        ),
        add_four(
            &from_identity().scale(c1),
            &from_outer_product(
                &p0.p,
                &vec2::add(&db0[1].scale(sum_inv), &dsum_inv[1].scale(b0)),
            ),
            &from_outer_product(
                &p1.p,
                &vec2::add(&db1[1].scale(sum_inv), &dsum_inv[1].scale(b1)),
            ),
            &from_outer_product(
                &p2.p,
                &vec2::add(&db2[1].scale(sum_inv), &dsum_inv[1].scale(b2)),
            ),
        ),
        add_four(
            &from_identity().scale(c2),
            &from_outer_product(
                &p0.p,
                &vec2::add(&db0[2].scale(sum_inv), &dsum_inv[2].scale(b0)),
            ),
            &from_outer_product(
                &p1.p,
                &vec2::add(&db1[2].scale(sum_inv), &dsum_inv[2].scale(b1)),
            ),
            &from_outer_product(
                &p2.p,
                &vec2::add(&db2[2].scale(sum_inv), &dsum_inv[2].scale(b2)),
            ),
        ),
    ];

    (cc, dcc)
}

#[test]
fn test_dw_circumcenter() {
    let p0 = [[0.1, 0.2], [1.3, 0.2], [0.3, 1.5]];
    let (cc0, dcc0) = wdw_circumcenter(&p0[0], &p0[1], &p0[2]);
    let eps = 1.0e-4;
    for i_node in 0..3 {
        for i_dim in 0..2 {
            let p1 = {
                let mut p1 = p0;
                p1[i_node][i_dim] += eps;
                p1
            };
            let (cc1, _dcc1) = wdw_circumcenter(&p1[0], &p1[1], &p1[2]);
            let dcc_num = cc1.sub(&cc0).scale(1. / eps);
            let mut b = [0f64; 2];
            b[i_dim] = 1.0;
            let dcc_ana = crate::mat2_col_major::mult_vec(&dcc0[i_node], &b);
            let diff = dcc_num.sub(&dcc_ana).norm();
            /*
            println!(
                "{}, {} --> {:?}, {:?}, {:?}",
                i_node, i_dim, dcc_num, dcc_ana, diff
            );
             */
            assert!(diff < 1.0e-4);
        }
    }
}
