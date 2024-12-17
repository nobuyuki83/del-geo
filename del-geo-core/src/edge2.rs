//! methods for 2D edge (line segment)

use crate::aabb;
use num_traits::AsPrimitive;

pub fn length<T>(ps: &[T; 2], pe: &[T; 2]) -> T
where
    T: num_traits::Float,
{
    let dx = ps[0] - pe[0];
    let dy = ps[1] - pe[1];
    (dx * dx + dy * dy).sqrt()
}

pub fn unit_edge_vector<T>(ps: &[T; 2], pe: &[T; 2]) -> [T; 2]
where
    T: num_traits::Float,
{
    let dx = pe[0] - ps[0];
    let dy = pe[1] - ps[1];
    let linv: T = T::one() / (dx * dx + dy * dy).sqrt();
    [dx * linv, dy * linv]
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
    let s0 = [0., 0.];
    let e0 = [1., 0.];
    let s1 = [0., -0.1];
    let e1 = [0.2, 0.1];
    let Some((r0, r1)) = intersection_edge2(&s0, &e0, &s1, &e1) else {
        panic!()
    };
    let p0 = crate::vec2::axpy(r0, &crate::vec2::sub(&e0, &s0), &s0);
    let p1 = crate::vec2::axpy(r1, &crate::vec2::sub(&e1, &s1), &s1);
    assert!(length(&p0, &[0.1, 0.0]) < 1.0e-5);
    assert!(length(&p1, &[0.1, 0.0]) < 1.0e-5);
}

pub fn winding_number<T>(ps: &[T; 2], pe: &[T; 2], po: &[T; 2]) -> T
where
    T: num_traits::Float + num_traits::FloatConst,
{
    let half = T::one() / (T::one() + T::one());
    let p0 = crate::vec2::sub(ps, po);
    let p1 = crate::vec2::sub(pe, po);
    let y: T = p1[1] * p0[0] - p1[0] * p0[1];
    let x: T = p0[0] * p1[0] + p0[1] * p1[1];
    y.atan2(x) * T::FRAC_1_PI() * half
}

/// Find the nearest point on a line segment to the origin(0,0)
/// Returns (k,v), where k is the coeffcient between [0,1], v is the point
pub fn nearest_origin<T>(ps: &[T; 2], pe: &[T; 2]) -> (T, [T; 2])
where
    T: num_traits::Float + 'static,
    f64: AsPrimitive<T>,
{
    let d = crate::vec2::sub(pe, ps);
    let a = crate::vec2::squared_length(&d);
    if a.is_zero() {
        return (
            0.5f64.as_(),
            std::array::from_fn(|i| (ps[i] + pe[i]) * 0.5f64.as_()),
        );
    }
    let b = crate::vec2::dot(&d, ps);
    let r0 = (-b / a).clamp(0f64.as_(), 1f64.as_());
    (
        r0,
        std::array::from_fn(|i| (-r0 + 1f64.as_()) * ps[i] + r0 * pe[i]),
    )
}

#[test]
fn test_nearest_origin() {
    let (_r, pm) = nearest_origin(&[-0.1, 1.0], &[1.0, 1.0]);
    assert!(crate::vec2::length(&crate::vec2::sub(&pm, &[0., 1.])) < 1.0e-5);
}

/// Find the nearest point on a line segment to point p
/// Returns (k,v), where k is the coeffcient, v is the point
pub fn nearest_point2(
    s: &[f32; 2], // source
    e: &[f32; 2], // end
    p: &[f32; 2],
) -> (f32, [f32; 2]) {
    use crate::vec2;
    let (r, p0) = nearest_origin(&vec2::sub(s, p), &vec2::sub(e, p));
    (r, [p0[0] + p[0], p0[1] + p[1]])
}

#[test]
fn test_nearest_point2() {
    let (_r, pm) = nearest_point2(&[-0.1, 1.0], &[1.0, 1.0], &[0.0, 0.3]);
    assert!(crate::vec2::length(&crate::vec2::sub(&pm, &[0., 1.])) < 1.0e-5);
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

pub fn area_left_side_of_edge2(aabb2: &[f32; 4], ps: &[f32; 2], pe: &[f32; 2]) -> f32 {
    let Some((tmin, tmax)) = aabb::intersections_against_line(aabb2, ps, &crate::vec2::sub(pe, ps))
    else {
        return 0f32;
    };
    // dbg!(tmin, tmax);
    let qs = crate::vec2::axpy(tmin, &crate::vec2::sub(pe, ps), ps);
    let qe = crate::vec2::axpy(tmax, &crate::vec2::sub(pe, ps), ps);
    // dbg!(qe, qs);
    if (qs[0] - aabb2[0]).abs() < f32::EPSILON {
        // left in
        assert!(aabb2[1] <= qs[1] && qs[1] < aabb2[3]);
        if (qe[1] - aabb2[1]).abs() < f32::EPSILON {
            // bottom out
            let l_s_y = qs[1] - aabb2[1];
            let l_e_x = qe[0] - aabb2[0];
            let a = (aabb2[3] - aabb2[1]) * (aabb2[2] - aabb2[0]) - l_s_y * l_e_x * 0.5;
            return a;
        } else if (qe[0] - aabb2[2]).abs() < f32::EPSILON {
            // right out
            assert!(aabb2[1] <= qe[1] && qe[1] < aabb2[3]);
            let l_s_y = aabb2[3] - qs[1];
            let l_e_y = aabb2[3] - qe[1];
            let a = (aabb2[2] - aabb2[0]) * (l_s_y + l_e_y) * 0.5;
            return a;
        } else if (qe[1] - aabb2[3]).abs() < f32::EPSILON {
            // top out
            let l_s_y = aabb2[3] - qs[1];
            let l_e_x = qe[0] - aabb2[0];
            return l_s_y * l_e_x * 0.5;
        }
    } else if (qs[1] - aabb2[1]).abs() < f32::EPSILON {
        // bottom in
        assert!(aabb2[0] <= qs[0] && qs[0] < aabb2[2]);
        if (qe[0] - aabb2[0]).abs() < f32::EPSILON {
            // left out
            let l_s_x = qs[0] - aabb2[0];
            let l_e_y = qe[1] - aabb2[1];
            let a = l_s_x * l_e_y * 0.5;
            return a;
        }
        if (qe[1] - aabb2[3]).abs() < f32::EPSILON {
            // top out
            assert!(aabb2[0] <= qe[0] && qe[0] < aabb2[2]);
            let l_s_x = qs[0] - aabb2[0];
            let l_e_x = qe[0] - aabb2[0];
            let a = (aabb2[3] - aabb2[1]) * (l_s_x + l_e_x) * 0.5;
            return a;
        }
        if (qe[0] - aabb2[2]).abs() < f32::EPSILON {
            // right out
            assert!(aabb2[1] <= qe[1] && qe[1] < aabb2[3]);
            let l_s_x = aabb2[2] - qs[0];
            let l_e_y = qe[1] - aabb2[1];
            let a = (aabb2[3] - aabb2[1]) * (aabb2[2] - aabb2[0]) - (l_s_x * l_e_y) * 0.5;
            return a;
        }
    } else if (qs[0] - aabb2[2]).abs() < f32::EPSILON {
        // right in
        assert!(aabb2[1] <= qs[1] && qs[1] < aabb2[3]);
        if (qe[0] - aabb2[0]).abs() < f32::EPSILON {
            // left out
            assert!(aabb2[1] <= qe[1] && qe[1] < aabb2[3]);
            let l_s_y = qs[1] - aabb2[1];
            let l_e_y = qe[1] - aabb2[1];
            let a = (aabb2[2] - aabb2[0]) * (l_s_y + l_e_y) * 0.5;
            return a;
        }
    } else if (qs[1] - aabb2[3]).abs() < f32::EPSILON && (qe[1] - aabb2[1]).abs() < f32::EPSILON {
        // up to down
        assert!(aabb2[0] <= qs[0] && qs[0] < aabb2[2]);
        assert!(aabb2[0] <= qe[0] && qe[0] < aabb2[2]);
        let l_s_x = aabb2[2] - qs[0];
        let l_e_x = aabb2[2] - qe[0];
        let a = (aabb2[3] - aabb2[1]) * (l_s_x + l_e_x) * 0.5;
        return a;
    }
    0f32
}

#[test]
fn test_area_of_aabb2_left_side_of_edge2() {
    let aabb2 = [0f32, 0.0, 1.0, 1.0];
    {
        let ps = [0.1f32, -0.5];
        let pe = [0.7f32, 1.5];
        let a = area_left_side_of_edge2(&aabb2, &ps, &pe); // down to up
        assert!((a - 0.4).abs() < 1.0e-7);
        let a = area_left_side_of_edge2(&aabb2, &pe, &ps); // up to down
        assert!((a - 0.6).abs() < 1.0e-7);
    }
    {
        let ps = [-0.5, 0.1];
        let pe = [1.5, 0.7f32];
        let a = area_left_side_of_edge2(&aabb2, &ps, &pe); // left to right
        assert!((a - 0.6).abs() < 1.0e-7);
        let a = area_left_side_of_edge2(&aabb2, &pe, &ps); // right to left
        assert!((a - 0.4).abs() < 1.0e-7);
    }
    {
        // left in
        let ps = [-0.2, 0.5];
        let a = area_left_side_of_edge2(&aabb2, &ps, &[0.5, -0.2]); // bottom out
        assert!((a - 0.955).abs() < 1.0e-7); // 0.955=1.0-0.3*0.3*0.5
        let a = area_left_side_of_edge2(&aabb2, &ps, &[0.5, 1.2]); // top out
        assert!((a - 0.045).abs() < 1.0e-7); // // 0.045=0.3*0.3*0.5
        let a = area_left_side_of_edge2(&aabb2, &ps, &[1.2, 0.7]); // right out
        assert!((a - 0.4).abs() < 1.0e-7);
    }
    {
        // bottom in
        let ps = [0.5f32, -0.2];
        let a = area_left_side_of_edge2(&aabb2, &ps, &[0.7, 1.2]); // top out
        assert!((a - 0.6).abs() < 1.0e-7);
        let a = area_left_side_of_edge2(&aabb2, &ps, &[-0.2, 0.5]); // left out
        assert!((a - 0.045).abs() < 1.0e-7);
        let a = area_left_side_of_edge2(&aabb2, &ps, &[1.2, 0.5]); // right out
        assert!((a - 0.955).abs() < 1.0e-7);
    }
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
    res
}
