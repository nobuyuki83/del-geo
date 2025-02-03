//! methods for 2D edge (line segment)

use num_traits::AsPrimitive;

pub fn length<T>(ps: &[T; 2], pe: &[T; 2]) -> T
where
    T: num_traits::Float,
{
    use crate::vec2::Vec2;
    let [dx, dy] = ps.sub(pe);
    (dx * dx + dy * dy).sqrt()
}

pub fn squared_length<T>(ps: &[T; 2], pe: &[T; 2]) -> T
where
    T: num_traits::Float,
{
    use crate::vec2::Vec2;
    let [dx, dy] = ps.sub(pe);
    dx * dx + dy * dy
}

pub fn unit_edge_vector<T>(ps: &[T; 2], pe: &[T; 2]) -> [T; 2]
where
    T: num_traits::Float,
{
    use crate::vec2::Vec2;
    let [dx, dy] = pe.sub(ps);
    let linv: T = T::one() / (dx * dx + dy * dy).sqrt();
    [dx * linv, dy * linv]
}

/// `ratio==0` should output `p0`
pub fn position_from_ratio<T>(p0: &[T; 2], p1: &[T; 2], ratio: T) -> [T; 2]
where
    T: num_traits::Float,
{
    let one = T::one();
    [
        (one - ratio) * p0[0] + ratio * p1[0],
        (one - ratio) * p0[1] + ratio * p1[1],
    ]
}

pub fn culling_intersection<T>(
    po_s0: &[T; 2],
    po_e0: &[T; 2],
    po_s1: &[T; 2],
    po_e1: &[T; 2],
) -> bool
where
    T: num_traits::Float + 'static + Copy,
    f64: AsPrimitive<T>,
{
    let min0x = po_s0[0].min(po_e0[0]);
    let max0x = po_s0[0].max(po_e0[0]);
    let min1x = po_s1[0].min(po_e1[0]);
    let max1x = po_s1[0].max(po_e1[0]);
    let min0y = po_s0[1].min(po_e0[1]);
    let max0y = po_s0[1].max(po_e0[1]);
    let min1y = po_s1[1].min(po_e1[1]);
    let max1y = po_s1[1].max(po_e1[1]);
    let len =
        ((max0x - min0x) + (max0y - min0y) + (max1x - min1x) + (max1y - min1y)) * 0.0001_f64.as_();
    if max1x + len < min0x {
        return false;
    }
    if max0x + len < min1x {
        return false;
    }
    if max1y + len < min0y {
        return false;
    }
    if max0y + len < min1y {
        return false;
    }
    true
}

///
/// * Return
///
/// Some((r0,r1)): ratio of edge
///
/// r0 * s0 + (1-r0) * e0 == r1 * s1 + (1-r1) * e1
pub fn intersection_edge2<T>(s0: &[T; 2], e0: &[T; 2], s1: &[T; 2], e1: &[T; 2]) -> Option<(T, T)>
where
    T: num_traits::Float,
{
    let area1 = crate::tri2::area(s0, e0, s1);
    let area2 = crate::tri2::area(s0, e0, e1);
    let area3 = crate::tri2::area(s1, e1, s0);
    let area4 = crate::tri2::area(s1, e1, e0);
    if area1 * area2 > T::zero() {
        return None;
    }
    if area3 * area4 > T::zero() {
        return None;
    }
    let r1 = area1 / (area1 - area2);
    let r0 = area3 / (area3 - area4);
    Some((r0, r1))
}

#[test]
fn test_intersection_edge2() {
    let s0 = [0f32, 0.];
    let e0 = [1., 0.];
    let s1 = [0., -0.1];
    let e1 = [0.2, 0.1];
    let Some((r0, r1)) = intersection_edge2(&s0, &e0, &s1, &e1) else {
        panic!()
    };
    use crate::vec2::Vec2;
    let p0 = crate::vec2::axpy(r0, &e0.sub(&s0), &s0);
    let p1 = crate::vec2::axpy(r1, &e1.sub(&s1), &s1);
    assert!(length(&p0, &[0.1, 0.0]) < 1.0e-5);
    assert!(length(&p1, &[0.1, 0.0]) < 1.0e-5);
    assert!((r0 - 0.1).abs() < 1e-5f32);
    assert!((r1 - 0.5).abs() < 1e-5f32);
}

/// # return
///
///  `(dls0: [T;2], dle0: [T;2], dlds1: [T;2], dlde1: [T;2])`
pub fn dldw_intersection_edge2<T>(
    s0: &[T; 2],
    e0: &[T; 2],
    s1: &[T; 2],
    e1: &[T; 2],
    dldr0: T,
    dldr1: T,
) -> ([T; 2], [T; 2], [T; 2], [T; 2])
where
    T: num_traits::Float + std::fmt::Debug,
{
    let one = T::one();
    let a1 = crate::tri2::area(s0, e0, s1);
    let a2 = crate::tri2::area(s0, e0, e1);
    let a3 = crate::tri2::area(s1, e1, s0);
    let a4 = crate::tri2::area(s1, e1, e0);
    //let r1 = a1 / (a1 - a2);
    //let r0 = a3 / (a3 - a4);
    let dlda1 = dldr1 * (one / (a1 - a2) - a1 / (a1 - a2).powi(2));
    let dlda2 = dldr1 * (a1 / (a1 - a2).powi(2));
    let dlda3 = dldr0 * (one / (a3 - a4) - a3 / (a3 - a4).powi(2));
    let dlda4 = dldr0 * (a3 / (a3 - a4).powi(2));
    let (dlds0_1, dlde0_1, dlds1_1) = crate::tri2::dldw_area(s0, e0, s1, dlda1);
    let (dlds0_2, dlde0_2, dlde1_2) = crate::tri2::dldw_area(s0, e0, e1, dlda2);
    let (dlds1_3, dlde1_3, dlds0_3) = crate::tri2::dldw_area(s1, e1, s0, dlda3);
    let (dlds1_4, dlde1_4, dlde0_4) = crate::tri2::dldw_area(s1, e1, e0, dlda4);
    use crate::vec2::Vec2;
    (
        dlds0_1.add(&dlds0_2).add(&dlds0_3),
        dlde0_1.add(&dlde0_2).add(&dlde0_4),
        dlds1_1.add(&dlds1_3).add(&dlds1_4),
        dlde1_2.add(&dlde1_3).add(&dlde1_4),
    )
}

#[test]
fn test_dldw_intersection_edge2() {
    let dldr0 = 0.3f64;
    let dldr1 = 0.65;
    let p0 = [[0.01, 0.03], [1.02, 0.05], [0.03, -0.102], [0.203, 0.105]];
    let (r0_0, r1_0) = intersection_edge2(&p0[0], &p0[1], &p0[2], &p0[3]).unwrap();
    let l0 = dldr0 * r0_0 + dldr1 * r1_0;
    let dl = dldw_intersection_edge2(&p0[0], &p0[1], &p0[2], &p0[3], dldr0, dldr1);
    let dl = [dl.0, dl.1, dl.2, dl.3];
    let eps = 1.0e-6;
    for inode in 0..4 {
        for idim in 0..2 {
            let p1 = {
                let mut p1 = p0;
                p1[inode][idim] += eps;
                p1
            };
            let (r0_1, r1_1) = intersection_edge2(&p1[0], &p1[1], &p1[2], &p1[3]).unwrap();
            let l1 = dldr0 * r0_1 + dldr1 * r1_1;
            let diff_num = (l1 - l0) / eps;
            let diff_ana = dl[inode][idim];
            let diff_abs = (diff_num - diff_ana).abs();
            assert!(diff_abs < 2.0e-5, "{}", diff_abs);
        }
    }
}

pub fn winding_number<T>(ps: &[T; 2], pe: &[T; 2], po: &[T; 2]) -> T
where
    T: num_traits::Float + num_traits::FloatConst,
{
    use crate::vec2::Vec2;
    let half = T::one() / (T::one() + T::one());
    let p0 = ps.sub(po);
    let p1 = pe.sub(po);
    let y: T = p1[1] * p0[0] - p1[0] * p0[1];
    let x: T = p0[0] * p1[0] + p0[1] * p1[1];
    y.atan2(x) * T::FRAC_1_PI() * half
}

/// Find the nearest point on a line segment to the origin(0,0)
/// Returns (k,v), where k is the coeffcient between `[0,1]`, v is the point
pub fn nearest_origin<T>(ps: &[T; 2], pe: &[T; 2]) -> (T, [T; 2])
where
    T: num_traits::Float + 'static,
    f64: AsPrimitive<T>,
{
    use crate::vec2::Vec2;
    let d = pe.sub(ps);
    let a = d.squared_norm();
    if a.is_zero() {
        return (
            0.5f64.as_(),
            std::array::from_fn(|i| (ps[i] + pe[i]) * 0.5f64.as_()),
        );
    }
    let b = d.dot(ps);
    let r0 = (-b / a).clamp(0f64.as_(), 1f64.as_());
    (
        r0,
        std::array::from_fn(|i| (-r0 + 1f64.as_()) * ps[i] + r0 * pe[i]),
    )
}

#[test]
fn test_nearest_origin() {
    let (_r, pm) = nearest_origin(&[-0.1, 1.0], &[1.0, 1.0]);
    use crate::vec2::Vec2;
    assert!(pm.sub(&[0., 1.]).norm() < 1.0e-5);
}

/// Find the nearest point on a line segment to point p
/// Returns (k,v), where k is the coeffcient, v is the point
pub fn nearest_point2<T>(
    s: &[T; 2], // source
    e: &[T; 2], // end
    p: &[T; 2],
) -> (T, [T; 2])
where
    T: num_traits::Float + 'static,
    f64: AsPrimitive<T>,
{
    use crate::vec2::Vec2;
    let (r, p0) = nearest_origin(&s.sub(p), &e.sub(p));
    (r, [p0[0] + p[0], p0[1] + p[1]])
}

#[test]
fn test_nearest_point2() {
    let (_r, pm) = nearest_point2(&[-0.1, 1.0], &[1.0, 1.0], &[0.0, 0.3]);
    use crate::vec2::Vec2;
    assert!(pm.sub(&[0., 1.]).norm() < 1.0e-5);
}

pub fn intersection_length_against_aabb2(ps: &[f32; 2], pe: &[f32; 2], aabb2: &[f32; 4]) -> f32 {
    // 0 min, 1 max
    let edge_range_x = [ps[0].min(pe[0]), ps[0].max(pe[0])];
    let aabb_range_x = [aabb2[0], aabb2[2]];
    let Some(dx) = crate::range::intersection_length(&aabb_range_x, &edge_range_x) else {
        return 0f32;
    };
    //
    let edge_range_y = [ps[1].min(pe[1]), ps[1].max(pe[1])];
    let aabb_range_y = [aabb2[1], aabb2[3]];
    let Some(dy) = crate::range::intersection_length(&aabb_range_y, &edge_range_y) else {
        return 0f32;
    };
    //
    (dx * dx + dy * dy).sqrt()
}

#[test]
fn test_intersection_length_against_aabb2() {
    use std::f32::consts::SQRT_2;
    // inside AABB
    let ps = [1.0, 1.0];
    let pe = [2.0, 2.0];
    let aabb = [0.0, 0.0, 3.0, 3.0];
    let length = intersection_length_against_aabb2(&ps, &pe, &aabb);
    assert!((length - SQRT_2).abs() < 1e-6);

    // outside AABB
    let ps = [0.0, 0.0];
    let pe = [0.5, 0.5];
    let aabb = [2.0, 2.0, 3.0, 3.0];
    let length = intersection_length_against_aabb2(&ps, &pe, &aabb);
    assert_eq!(length, 0.0);

    // outside AABB
    let ps = [2.0, 0.0];
    let pe = [3.0, 0.5];
    let aabb = [2.0, 2.0, 3.0, 3.0];
    let length = intersection_length_against_aabb2(&ps, &pe, &aabb);
    assert_eq!(length, 0.0);

    // partially inside AABB
    let ps = [0.0, 0.0];
    let pe = [2.0, 2.0];
    let aabb = [1.0, 1.0, 3.0, 3.0];
    let length = intersection_length_against_aabb2(&ps, &pe, &aabb);
    assert!((length - SQRT_2).abs() < 1e-6);

    // segment with coordinates in reverse order
    let ps = [3.0, 3.0];
    let pe = [1.0, 1.0];
    let aabb = [0.0, 0.0, 2.0, 2.0];
    let length = intersection_length_against_aabb2(&ps, &pe, &aabb);
    assert!((length - SQRT_2).abs() < 1e-6);

    // zero-length segment
    let ps = [1.0, 1.0];
    let pe = [1.0, 1.0];
    let aabb = [0.0, 0.0, 2.0, 2.0];
    let length = intersection_length_against_aabb2(&ps, &pe, &aabb);
    assert_eq!(length, 0.0);

    // segment crossing aabb diagonally
    let ps = [0.0, 0.0];
    let pe = [3.0, 3.0];
    let aabb = [1.0, 1.0, 2.0, 2.0];
    let length = intersection_length_against_aabb2(&ps, &pe, &aabb);
    assert!((length - SQRT_2).abs() < 1e-6);
}

pub fn overlapping_pixels_dda<Real>(
    (img_width, img_height): (usize, usize),
    p0: &[Real; 2],
    p1: &[Real; 2],
) -> Vec<usize>
where
    Real: num_traits::Float + 'static + Copy + AsPrimitive<usize> + std::fmt::Debug,
    usize: AsPrimitive<Real>,
{
    let width_f: Real = img_width.as_();
    let height_f: Real = img_height.as_();
    let zero = Real::zero();
    let dx = p1[0] - p0[0];
    let dy = p1[1] - p0[1];
    let step = if dx.abs() > dy.abs() {
        dx.abs()
    } else {
        dy.abs()
    };
    let slope_y = dy / step;
    let slope_x = dx / step;
    let mut x = p0[0];
    let mut y = p0[1];
    let mut res = vec![];
    while (x - p0[0]).abs() <= (p1[0] - p0[0]).abs() && (y - p0[1]).abs() <= (p1[1] - p0[1]).abs() {
        if x >= zero && x < width_f && y >= zero && y < height_f {
            let ix: usize = x.as_();
            let iy: usize = y.as_();
            res.push(iy * img_width + ix);
        }
        x = x + slope_x;
        y = y + slope_y;
    }
    {
        let ix: usize = x.as_();
        let iy: usize = y.as_();
        res.push(iy * img_width + ix);
    }
    res
}
