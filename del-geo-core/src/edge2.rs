//! methods for 2D edge (line segment)

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

pub fn intersection_edge2<T>(
    po_s0: &[T; 2],
    po_e0: &[T; 2],
    po_s1: &[T; 2],
    po_e1: &[T; 2],
) -> Option<(T, T)>
where
    T: num_traits::Float + 'static + Copy,
    f64: AsPrimitive<T>,
{
    let area1 = crate::tri2::area(po_s0, po_e0, po_s1);
    let area2 = crate::tri2::area(po_s0, po_e0, po_e1);
    let area3 = crate::tri2::area(po_s1, po_e1, po_s0);
    let area4 = crate::tri2::area(po_s1, po_e1, po_e0);
    if area1 * area2 > 0_f64.as_() {
        return None;
    }
    if area3 * area4 > 0_f64.as_() {
        return None;
    }
    let r1 = area1 / (area1 - area2);
    let r0 = area3 / (area3 - area4);
    Some((r0, r1))
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
    let seg_range_x = [ps[0].min(pe[0]), ps[0].max(pe[0])];
    let aabb_range_x = [aabb2[0], aabb2[2]];

    // no intersction
    if seg_range_x[0] > aabb_range_x[1] || seg_range_x[1] < aabb_range_x[0] {
        return 0f32;
    }

    let seg_range_y = [ps[1].min(pe[1]), ps[1].max(pe[1])];
    let aabb_range_y = [aabb2[1], aabb2[3]];

    if seg_range_y[0] > aabb_range_y[1] || seg_range_y[1] < aabb_range_y[0] {
        return 0f32;
    }

    // fully inside aabb
    if seg_range_x[1] <= aabb_range_x[1]
        && seg_range_x[0] >= aabb_range_x[0]
        && seg_range_y[1] <= aabb_range_y[1]
        && seg_range_y[0] >= aabb_range_y[0]
    {
        return crate::vec2::length(&crate::vec2::sub(pe, ps));
    }

    let dx = range_intersection_length(&aabb_range_x, &seg_range_x);
    let dy = range_intersection_length(&aabb_range_y, &seg_range_y);
    (dx * dx + dy * dy).sqrt()
}

pub fn range_intersection_length(r1: &[f32; 2], r2: &[f32; 2]) -> f32 {
    // separeted
    if r1[1] <= r2[0] || r1[0] >= r2[1] {
        return 0f32;
    };

    // contains
    if r1[0] >= r2[0] && r1[1] <= r2[1] {
        return r1[1] - r1[0];
    } else if r1[0] <= r2[0] && r1[1] >= r2[1] {
        return r2[1] - r2[0];
    };

    if r1[1] > r2[1] {
        r2[1] - r1[0]
    } else if r2[1] > r1[1] {
        r1[1] - r2[0]
    } else {
        (r1[1] - r2[0]).min(r1[1] - r1[0])
    }
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

    // partially inside AABB
    let ps = [0.0, 0.0];
    let pe = [2.0, 2.0];
    let aabb = [1.0, 1.0, 3.0, 3.0];
    let length = intersection_length_against_aabb2(&ps, &pe, &aabb);
    assert!((length - SQRT_2).abs() < 1e-6);

    // touching border
    let ps = [1.0, 0.0];
    let pe = [1.0, 3.0];
    let aabb = [1.0, 1.0, 2.0, 2.0];
    let length = intersection_length_against_aabb2(&ps, &pe, &aabb);
    assert!((length - 1.0).abs() < 1e-6);

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
